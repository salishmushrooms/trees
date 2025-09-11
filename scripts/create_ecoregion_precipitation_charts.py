#!/usr/bin/env python3
"""
Create monthly precipitation normal histograms for Washington L3 ecoregions.

This script:
1. Loads Washington ecoregion shapefile (L4 boundaries)
2. Extracts monthly PRISM precipitation normals for each L3 region
3. Creates histograms showing precipitation distribution by month
4. Outputs charts to outputs/ecoregions/precipitation/

Usage:
    python create_ecoregion_precipitation_charts.py --l3code 1  # Test with one region
    python create_ecoregion_precipitation_charts.py --all      # Process all regions
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.mask import mask
import matplotlib.pyplot as plt
from pathlib import Path

def load_ecoregions(shapefile_path):
    """Load Washington ecoregions shapefile."""
    print(f"Loading ecoregions from {shapefile_path}")
    gdf = gpd.read_file(shapefile_path)
    print(f"Loaded {len(gdf)} ecoregions")
    print(f"Available columns: {gdf.columns.tolist()}")
    print(f"L3 codes available: {sorted(gdf['US_L3CODE'].unique())}")
    return gdf

def extract_precipitation_for_region(region_geom, region_crs, prism_raster_path):
    """Extract precipitation values for a specific ecoregion geometry."""
    try:
        with rasterio.open(prism_raster_path) as src:
            print(f"  Raster CRS: {src.crs}")
            print(f"  Region CRS: {region_crs}")
            
            # Reproject region geometry to match raster CRS if needed
            if str(region_crs) != str(src.crs):
                from shapely.geometry import shape
                import pyproj
                from pyproj import Transformer
                
                # Create transformer
                transformer = Transformer.from_crs(region_crs, src.crs, always_xy=True)
                
                # Transform geometry
                from shapely.ops import transform
                region_geom_transformed = transform(transformer.transform, region_geom)
                print(f"  Transformed geometry bounds: {region_geom_transformed.bounds}")
            else:
                region_geom_transformed = region_geom
            
            # Check if geometry is valid and not empty
            if not region_geom_transformed.is_valid or region_geom_transformed.is_empty:
                print(f"  Invalid or empty geometry after transformation")
                return np.array([])
            
            # Mask raster with region geometry
            out_image, out_transform = mask(src, [region_geom_transformed], crop=True, nodata=src.nodata)
            out_image = out_image[0]  # Get first band
            
            # Filter out nodata values
            valid_data = out_image[out_image != src.nodata]
            if len(valid_data) == 0:
                return np.array([])
            
            return valid_data
    except Exception as e:
        print(f"Error processing {prism_raster_path}: {e}")
        return np.array([])

def create_precipitation_chart(monthly_data, region_name, l3_code, output_dir):
    """Create bar chart showing monthly precipitation normals for a region."""
    months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 
              'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    
    # Calculate monthly statistics and convert from mm to inches
    # 1 mm = 0.0393701 inches
    mm_to_inches = 0.0393701
    
    monthly_stats = []
    for month_idx, month_name in enumerate(months, 1):
        if month_idx in monthly_data:
            data = monthly_data[month_idx]
            if len(data) > 0:
                monthly_stats.append({
                    'month': month_name,
                    'mean': np.mean(data) * mm_to_inches,
                    'median': np.median(data) * mm_to_inches,
                    'min': np.min(data) * mm_to_inches,
                    'max': np.max(data) * mm_to_inches,
                    'std': np.std(data) * mm_to_inches
                })
    
    if not monthly_stats:
        print(f"No valid data for region {region_name}")
        return None
    
    # Create single bar chart
    fig, ax = plt.subplots(1, 1, figsize=(12, 6))
    
    df_stats = pd.DataFrame(monthly_stats)
    ax.bar(df_stats['month'], df_stats['mean'], color='skyblue', alpha=0.7, edgecolor='navy')
    ax.set_title(f'{region_name} - Monthly Precipitation Normals (1991-2020)', fontsize=14, fontweight='bold')
    ax.set_ylabel('Precipitation (inches)', fontsize=12)
    ax.grid(True, alpha=0.3)
    
    # Add value labels on bars (show 1 decimal place for inches)
    for i, v in enumerate(df_stats['mean']):
        ax.text(i, v + max(df_stats['mean']) * 0.01, f'{v:.1f}', 
                ha='center', va='bottom', fontweight='bold')
    
    plt.tight_layout()
    
    # Save the plot with _inches suffix
    output_file = output_dir / f"L3_{l3_code}_{region_name.replace(' ', '_').replace('/', '_')}_precipitation_inches.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved chart: {output_file}")
    plt.close()
    
    # Save individual statistics to CSV with _inches suffix
    stats_file = output_dir / f"L3_{l3_code}_{region_name.replace(' ', '_').replace('/', '_')}_precipitation_stats_inches.csv"
    df_stats.to_csv(stats_file, index=False)
    print(f"Saved statistics: {stats_file}")
    
    # Return stats with region info for comparison CSV
    df_stats['l3_code'] = l3_code
    df_stats['region_name'] = region_name
    return df_stats

def main():
    parser = argparse.ArgumentParser(description='Create ecoregion precipitation charts')
    parser.add_argument('--l3code', type=int, help='Process specific L3 code (for testing)')
    parser.add_argument('--all', action='store_true', help='Process all L3 regions')
    args = parser.parse_args()
    
    if not args.l3code and not args.all:
        parser.error('Must specify either --l3code or --all')
    
    # Paths
    base_dir = Path(__file__).parent.parent
    shapefile_path = base_dir / "data/ecoregions/wa_eco_l4/wa_eco_l4.shp"
    prism_dir = base_dir / "data/raw/PRISM_ppt_30yr_normal_800mM4_all_bil"
    output_dir = base_dir / "outputs/ecoregions/precipitation"
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load ecoregions
    ecoregions = load_ecoregions(shapefile_path)
    
    # Determine which L3 codes to process
    if args.l3code:
        l3_codes = [str(args.l3code)]  # Convert to string to match shapefile
        print(f"Processing L3 code: {args.l3code}")
    else:
        l3_codes = sorted(ecoregions['US_L3CODE'].unique())
        print(f"Processing all L3 codes: {l3_codes}")
    
    # Initialize list to collect all region statistics
    all_region_stats = []
    
    # Process each L3 region
    for l3_code in l3_codes:
        print(f"\n=== Processing L3 Code: {l3_code} ===")
        
        # Filter ecoregions for this L3 code
        l3_regions = ecoregions[ecoregions['US_L3CODE'] == l3_code]
        
        if len(l3_regions) == 0:
            print(f"No regions found for L3 code {l3_code}")
            continue
        
        # Get region name from first record
        region_name = l3_regions.iloc[0]['US_L3NAME']
        print(f"Region name: {region_name}")
        print(f"Number of L4 subregions: {len(l3_regions)}")
        
        # Dissolve L4 regions to create single L3 geometry
        l3_dissolved = l3_regions.dissolve(by='US_L3CODE').iloc[0]
        region_geom = l3_dissolved.geometry
        region_crs = l3_regions.crs
        
        print(f"Region bounds: {region_geom.bounds}")
        print(f"Region CRS: {region_crs}")
        
        # Extract monthly precipitation data
        monthly_data = {}
        
        for month in range(1, 13):
            month_str = f"{month:02d}"
            raster_file = prism_dir / f"PRISM_ppt_30yr_normal_800mM4_{month_str}_bil.bil"
            
            if not raster_file.exists():
                print(f"Warning: Raster file not found: {raster_file}")
                continue
            
            print(f"Extracting precipitation data for month {month_str}...")
            precip_data = extract_precipitation_for_region(region_geom, region_crs, raster_file)
            
            if len(precip_data) > 0:
                monthly_data[month] = precip_data
                print(f"  Month {month_str}: {len(precip_data)} pixels, "
                      f"mean={np.mean(precip_data):.1f}mm, "
                      f"range={np.min(precip_data):.1f}-{np.max(precip_data):.1f}mm")
            else:
                print(f"  Month {month_str}: No valid data")
        
        # Create chart and collect stats
        if monthly_data:
            region_stats = create_precipitation_chart(monthly_data, region_name, l3_code, output_dir)
            if region_stats is not None:
                all_region_stats.append(region_stats)
        else:
            print(f"No precipitation data available for L3 code {l3_code}")
    
    # Create comparison CSV with all regions
    if all_region_stats:
        print(f"\n=== Creating Comparison CSV ===")
        
        # Combine all region statistics
        comparison_df = pd.concat(all_region_stats, ignore_index=True)
        
        # Pivot to create a matrix with regions as rows and months as columns
        comparison_pivot = comparison_df.pivot(index=['l3_code', 'region_name'], 
                                              columns='month', 
                                              values='mean')
        
        # Reorder columns to be in chronological order
        month_order = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 
                       'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
        comparison_pivot = comparison_pivot.reindex(columns=month_order)
        
        # Save comparison CSV with _inches suffix
        comparison_file = output_dir / "all_ecoregions_precipitation_comparison_inches.csv"
        comparison_pivot.to_csv(comparison_file)
        print(f"Saved comparison CSV: {comparison_file}")
        
        # Also create a detailed version with all statistics
        detailed_file = output_dir / "all_ecoregions_precipitation_detailed_inches.csv"
        comparison_df.to_csv(detailed_file, index=False)
        print(f"Saved detailed CSV: {detailed_file}")
    else:
        print("No regional statistics were collected")

if __name__ == "__main__":
    main()