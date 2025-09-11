#!/usr/bin/env python3
"""
Tree Species Summary by Ecoregions

Creates tree species summaries grouped by US_L3CODE and US_L4CODE using FIA plot data
and ecoregions vector data. Uses only the latest survey data per plot to avoid double counting.

Usage:
    python ecoregions_tree_summary.py --state WA --level L3
    python ecoregions_tree_summary.py --state WA --level L4
"""

import sqlite3
import pandas as pd
import geopandas as gpd
from pathlib import Path
import argparse
import sys
import rasterio
import rasterio.mask
from rasterio.features import shapes
import numpy as np

def get_species_mapping(db_path):
    """Get species code to common name mapping from FIA database"""
    conn = sqlite3.connect(db_path)
    
    species_query = """
    SELECT SPCD, COMMON_NAME, SCIENTIFIC_NAME 
    FROM REF_SPECIES 
    WHERE COMMON_NAME IS NOT NULL
    ORDER BY SPCD;
    """
    
    species_df = pd.read_sql_query(species_query, conn)
    conn.close()
    
    # Create mapping dictionary
    species_map = dict(zip(species_df['SPCD'], species_df['COMMON_NAME']))
    return species_map

def extract_latest_plot_data(db_path, state_code=53):
    """Extract latest survey data per plot for specified state"""
    conn = sqlite3.connect(db_path)
    
    # Get latest survey year per plot
    latest_survey_query = f"""
    WITH latest_plots AS (
        SELECT 
            STATECD, UNITCD, COUNTYCD, PLOT,
            MAX(INVYR) as latest_invyr
        FROM PLOT 
        WHERE STATECD = {state_code} 
        AND PLOT_STATUS_CD = 1  -- Forested plots only
        GROUP BY STATECD, UNITCD, COUNTYCD, PLOT
    )
    SELECT 
        p.CN as PLT_CN,
        p.STATECD, p.UNITCD, p.COUNTYCD, p.PLOT,
        p.INVYR, p.LAT, p.LON, p.ELEV
    FROM PLOT p
    INNER JOIN latest_plots lp ON (
        p.STATECD = lp.STATECD AND 
        p.UNITCD = lp.UNITCD AND 
        p.COUNTYCD = lp.COUNTYCD AND 
        p.PLOT = lp.PLOT AND 
        p.INVYR = lp.latest_invyr
    )
    WHERE p.PLOT_STATUS_CD = 1;
    """
    
    plots_df = pd.read_sql_query(latest_survey_query, conn)
    
    # Get tree data for these latest plots
    plt_cn_list = "','".join(plots_df['PLT_CN'].astype(str))
    
    trees_query = f"""
    SELECT 
        t.PLT_CN,
        t.SPCD,
        COUNT(*) as tree_count
    FROM TREE t
    WHERE t.PLT_CN IN ('{plt_cn_list}')
    AND t.STATUSCD = 1  -- Live trees only
    GROUP BY t.PLT_CN, t.SPCD;
    """
    
    trees_df = pd.read_sql_query(trees_query, conn)
    conn.close()
    
    return plots_df, trees_df

def load_ecoregions(state='WA'):
    """Load ecoregions shapefile for specified state"""
    eco_path = f"data/ecoregions/{state.lower()}_eco_l4/{state.lower()}_eco_l4.shp"
    
    if not Path(eco_path).exists():
        raise FileNotFoundError(f"Ecoregions file not found: {eco_path}")
    
    eco_gdf = gpd.read_file(eco_path)
    
    # Ensure we have the required columns
    required_cols = ['US_L3CODE', 'US_L3NAME', 'US_L4CODE', 'US_L4NAME']
    missing_cols = [col for col in required_cols if col not in eco_gdf.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns in ecoregions data: {missing_cols}")
    
    return eco_gdf

def spatial_join_plots_ecoregions(plots_df, eco_gdf):
    """Perform spatial join between plots and ecoregions"""
    # Convert plots to GeoDataFrame
    plots_gdf = gpd.GeoDataFrame(
        plots_df,
        geometry=gpd.points_from_xy(plots_df.LON, plots_df.LAT),
        crs='EPSG:4326'
    )
    
    # Ensure both have same CRS
    if eco_gdf.crs != plots_gdf.crs:
        eco_gdf = eco_gdf.to_crs(plots_gdf.crs)
    
    # Spatial join
    joined_gdf = gpd.sjoin(plots_gdf, eco_gdf, how='left', predicate='within')
    
    return joined_gdf

def calculate_tree_cover_by_ecoregion(eco_gdf, level='L3'):
    """Calculate average tree cover percentage for each ecoregion"""
    
    # Check for tree cover raster
    tcc_path = "data/raw/nlcd_tcc_conus_2021_v2021-4.tif"
    if not Path(tcc_path).exists():
        print(f"WARNING: Tree cover raster not found at {tcc_path}")
        return None
    
    tcc_results = []
    
    print(f"Calculating tree cover for {len(eco_gdf)} ecoregions...")
    
    with rasterio.open(tcc_path) as src:
        # Group by ecoregion code to handle duplicates
        for (code, name), group in eco_gdf.groupby([f'US_{level}CODE', f'US_{level}NAME']):
            try:
                # Merge all polygons for this ecoregion
                try:
                    merged_geom = group.geometry.union_all()
                except AttributeError:
                    # Fallback for older versions of geopandas
                    merged_geom = group.geometry.unary_union
                
                # Ensure it's in the same CRS as the raster
                if eco_gdf.crs != src.crs:
                    eco_single = gpd.GeoDataFrame([1], geometry=[merged_geom], crs=eco_gdf.crs)
                    eco_single = eco_single.to_crs(src.crs)
                    merged_geom = eco_single.geometry.iloc[0]
                
                # Mask the raster with the ecoregion polygon
                masked_data, masked_transform = rasterio.mask.mask(
                    src, [merged_geom], crop=True, nodata=src.nodata
                )
                
                # Calculate average tree cover (excluding nodata values)
                valid_data = masked_data[masked_data != src.nodata]
                if len(valid_data) > 0:
                    avg_tcc = float(np.mean(valid_data))
                    area_sqkm = group.geometry.to_crs('EPSG:3857').area.sum() / 1e6  # Convert to km²
                else:
                    avg_tcc = 0.0
                    area_sqkm = 0.0
                
                tcc_results.append({
                    f'US_{level}CODE': code,
                    f'US_{level}NAME': name,
                    'avg_tree_cover_pct': round(avg_tcc, 1),
                    'area_sqkm': round(area_sqkm, 1)
                })
                
            except Exception as e:
                print(f"Warning: Could not calculate tree cover for {name} ({code}): {e}")
                tcc_results.append({
                    f'US_{level}CODE': code,
                    f'US_{level}NAME': name,
                    'avg_tree_cover_pct': np.nan,
                    'area_sqkm': np.nan
                })
    
    return pd.DataFrame(tcc_results)

def create_species_summary(joined_plots, trees_df, species_map, eco_gdf, level='L3'):
    """Create tree species summary by ecoregion"""
    
    # Determine which columns to keep from joined_plots
    base_cols = ['PLT_CN', f'US_{level}CODE', f'US_{level}NAME']
    
    # If we're working with L4 data, also include L3 information
    if level == 'L4' and 'US_L3CODE' in joined_plots.columns:
        base_cols.extend(['US_L3CODE', 'US_L3NAME'])
    
    # Merge plots with tree data
    plot_trees = pd.merge(
        joined_plots[base_cols], 
        trees_df, 
        on='PLT_CN'
    )
    
    # Add species common names
    plot_trees['COMMON_NAME'] = plot_trees['SPCD'].map(species_map)
    
    # Remove unmapped species
    plot_trees = plot_trees.dropna(subset=['COMMON_NAME'])
    
    # Define groupby columns
    group_cols = [f'US_{level}CODE', f'US_{level}NAME', 'COMMON_NAME']
    
    # Add L3 columns if we're working with L4 data and they exist
    if level == 'L4' and 'US_L3CODE' in plot_trees.columns:
        group_cols = ['US_L3CODE', 'US_L3NAME'] + group_cols
    
    # Group by ecoregion and species
    summary = plot_trees.groupby(group_cols).agg({
        'tree_count': 'sum',
        'PLT_CN': 'nunique'  # unique plots per species in ecoregion
    }).reset_index()
    
    # Rename columns for clarity
    summary = summary.rename(columns={'PLT_CN': 'plots_with_species'})
    
    # Get total plots per ecoregion (using the primary level columns)
    total_plots_cols = [f'US_{level}CODE', f'US_{level}NAME']
    if level == 'L4' and 'US_L3CODE' in joined_plots.columns:
        total_plots_cols = ['US_L3CODE', 'US_L3NAME'] + total_plots_cols
    
    total_plots = joined_plots.groupby(total_plots_cols).size().reset_index(name='total_plots')
    
    # Merge summary with total plots
    merge_cols = [f'US_{level}CODE', f'US_{level}NAME']
    if level == 'L4' and 'US_L3CODE' in summary.columns:
        merge_cols = ['US_L3CODE', 'US_L3NAME'] + merge_cols
    
    summary = pd.merge(summary, total_plots, on=merge_cols)
    
    # Calculate tree cover by ecoregion
    tcc_summary = calculate_tree_cover_by_ecoregion(eco_gdf, level)
    if tcc_summary is not None:
        summary = pd.merge(summary, tcc_summary, on=[f'US_{level}CODE', f'US_{level}NAME'], how='left')
    
    return summary

def format_output(summary, level='L3'):
    """Format summary into requested output format"""
    
    output_lines = []
    
    # Group by ecoregion
    for (code, name), group in summary.groupby([f'US_{level}CODE', f'US_{level}NAME']):
        total_plots = group['total_plots'].iloc[0]
        
        # Get tree cover and area info (should be same for all rows in group)
        avg_tcc = group['avg_tree_cover_pct'].iloc[0] if 'avg_tree_cover_pct' in group.columns else None
        area_sqkm = group['area_sqkm'].iloc[0] if 'area_sqkm' in group.columns else None
        
        output_lines.append(f"US_{level}CODE {code}")
        output_lines.append(f"US_{level}NAME {name}")
        output_lines.append(f"Total Plots in US_{level}CODE = {code} area: {total_plots}")
        
        if avg_tcc is not None and not pd.isna(avg_tcc):
            output_lines.append(f"Average Tree Cover: {avg_tcc}%")
        if area_sqkm is not None and not pd.isna(area_sqkm):
            output_lines.append(f"Total Area: {area_sqkm} km²")
        
        output_lines.append("")
        
        # Sort species by tree count (descending)
        species_sorted = group.sort_values('tree_count', ascending=False)
        
        for _, row in species_sorted.iterrows():
            species_name = row['COMMON_NAME'].title()
            tree_count = int(row['tree_count'])
            plots_count = int(row['plots_with_species'])
            
            output_lines.append(f"{species_name} - {tree_count} trees (found in {plots_count} plots)")
        
        output_lines.append("\n" + "="*60 + "\n")
    
    return "\n".join(output_lines)

def format_csv_output(summary, level='L3'):
    """Format summary as CSV for analysis and charting"""
    
    # Create a flat structure for CSV
    csv_data = []
    
    for _, row in summary.iterrows():
        csv_record = {
            f'{level}_CODE': row[f'US_{level}CODE'],
            f'{level}_NAME': row[f'US_{level}NAME'],
            'species_common_name': row['COMMON_NAME'],
            'tree_count': int(row['tree_count']),
            'plots_with_species': int(row['plots_with_species']),
            'total_plots_in_ecoregion': int(row['total_plots'])
        }
        
        # Add L3 information if we're working with L4 data
        if level == 'L4':
            if 'US_L3CODE' in row and not pd.isna(row['US_L3CODE']):
                csv_record['L3_CODE'] = row['US_L3CODE']
            else:
                csv_record['L3_CODE'] = None
                
            if 'US_L3NAME' in row and not pd.isna(row['US_L3NAME']):
                csv_record['L3_NAME'] = row['US_L3NAME']
            else:
                csv_record['L3_NAME'] = None
        
        # Add tree cover and area if available
        if 'avg_tree_cover_pct' in row and not pd.isna(row['avg_tree_cover_pct']):
            csv_record['avg_tree_cover_pct'] = row['avg_tree_cover_pct']
        else:
            csv_record['avg_tree_cover_pct'] = None
            
        if 'area_sqkm' in row and not pd.isna(row['area_sqkm']):
            csv_record['area_sqkm'] = row['area_sqkm']
        else:
            csv_record['area_sqkm'] = None
        
        # Calculate derived metrics
        csv_record['pct_plots_with_species'] = round(
            (csv_record['plots_with_species'] / csv_record['total_plots_in_ecoregion']) * 100, 1
        )
        
        csv_data.append(csv_record)
    
    # Convert to DataFrame and sort by ecoregion, then by tree count
    csv_df = pd.DataFrame(csv_data)
    csv_df = csv_df.sort_values([f'{level}_CODE', 'tree_count'], ascending=[True, False])
    
    return csv_df

def main():
    parser = argparse.ArgumentParser(description='Generate tree species summaries by ecoregions')
    parser.add_argument('--state', default='WA', choices=['WA', 'OR'], 
                       help='State to analyze (default: WA)')
    parser.add_argument('--level', default='L3', choices=['L3', 'L4'],
                       help='Ecoregion level to summarize by (default: L3)')
    parser.add_argument('--output', type=str,
                       help='Output file path (default: print to stdout)')
    parser.add_argument('--format', default='txt', choices=['txt', 'csv'],
                       help='Output format: txt for formatted text, csv for spreadsheet (default: txt)')
    
    args = parser.parse_args()
    
    # Set up file paths
    db_path = f"data/raw/trees_SQLite_FIADB_{args.state}.db"
    state_code = 53 if args.state == 'WA' else 41  # Washington=53, Oregon=41
    
    if not Path(db_path).exists():
        print(f"ERROR: FIA database not found: {db_path}")
        sys.exit(1)
    
    print(f"Processing {args.state} tree species summary by US_{args.level}CODE...")
    
    # Step 1: Get species mapping
    print("Getting species mapping...")
    species_map = get_species_mapping(db_path)
    
    # Step 2: Extract latest plot data
    print("Extracting latest survey data per plot...")
    plots_df, trees_df = extract_latest_plot_data(db_path, state_code)
    print(f"Found {len(plots_df)} plots with latest survey data")
    print(f"Found {len(trees_df)} plot-species combinations")
    
    # Step 3: Load ecoregions
    print("Loading ecoregions data...")
    eco_gdf = load_ecoregions(args.state)
    print(f"Loaded {len(eco_gdf)} ecoregion polygons")
    
    # Step 4: Spatial join
    print("Performing spatial join...")
    joined_plots = spatial_join_plots_ecoregions(plots_df, eco_gdf)
    
    # Check for plots without ecoregion assignment
    unassigned = joined_plots[f'US_{args.level}CODE'].isna().sum()
    if unassigned > 0:
        print(f"WARNING: {unassigned} plots could not be assigned to ecoregions")
    
    # Step 5: Create summary
    print("Creating species summary...")
    summary = create_species_summary(joined_plots, trees_df, species_map, eco_gdf, args.level)
    
    # Step 6: Format and output results
    if args.format == 'csv':
        csv_output = format_csv_output(summary, args.level)
        
        if args.output:
            csv_output.to_csv(args.output, index=False)
            print(f"CSV summary written to: {args.output}")
            print(f"Columns: {list(csv_output.columns)}")
            print(f"Records: {len(csv_output)} species-ecoregion combinations")
        else:
            print("\n" + "="*60)
            print(f"TREE SPECIES CSV DATA - US_{args.level}CODE - {args.state}")
            print("="*60 + "\n")
            print(csv_output.to_string(index=False))
    else:
        # Original text format
        output_text = format_output(summary, args.level)
        
        if args.output:
            with open(args.output, 'w') as f:
                f.write(output_text)
            print(f"Summary written to: {args.output}")
        else:
            print("\n" + "="*60)
            print(f"TREE SPECIES SUMMARY BY US_{args.level}CODE - {args.state}")
            print("="*60 + "\n")
            print(output_text)

if __name__ == "__main__":
    main()