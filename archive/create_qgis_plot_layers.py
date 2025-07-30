#!/usr/bin/env python3
"""
Create QGIS-Compatible Species Carbon Layer

This script creates a specialized layer showing total above ground carbon
for each major tree species by plot in Washington state. Perfect for
visualizing species abundance and doing species-level calculations in QGIS.

KEY FIELD MEANINGS:
- *_CARBON_AG: Above ground carbon stored in each species (metric tons)
- *_BASAL_AREA: Cross-sectional area of tree trunks at breast height (sq ft)
- *_AVG_DBH: Average diameter at breast height for species (inches)
- SHANNON_DIVERSITY: Species diversity index (higher = more diverse)
- Only includes the most recent survey per unique plot location
"""

import sqlite3
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import numpy as np
from pathlib import Path

# Database connection
DB_PATH = "data/raw/trees_SQLite_FIADB_WA.db"
OUTPUT_DIR = Path("outputs")
OUTPUT_DIR.mkdir(exist_ok=True)

def create_species_carbon_layer():
    """
    Create species-level carbon layer showing total above ground carbon by species per plot
    """
    print("Creating species-level carbon layer for QGIS...")
    print("="*60)
    
    # First, get the most common species in Washington to focus on major species
    species_query = """
    SELECT 
        t.SPCD,
        s.COMMON_NAME,
        s.SCIENTIFIC_NAME,
        COUNT(DISTINCT t.PLT_CN) as PLOT_COUNT,
        COUNT(*) as TREE_COUNT,
        SUM(t.CARBON_AG) as TOTAL_CARBON
    FROM TREE t
    JOIN REF_SPECIES s ON t.SPCD = s.SPCD
    JOIN PLOT p ON t.PLT_CN = p.CN
    WHERE p.PLOT_STATUS_CD = 1
        AND t.STATUSCD = 1  -- Live trees only
        AND t.CARBON_AG IS NOT NULL
        AND t.CARBON_AG > 0
    GROUP BY t.SPCD, s.COMMON_NAME, s.SCIENTIFIC_NAME
    HAVING PLOT_COUNT >= 25  -- Species found in at least 25 plots
    ORDER BY TOTAL_CARBON DESC
    LIMIT 15;  -- Top 15 species by total carbon
    """
    
    conn = sqlite3.connect(DB_PATH)
    species_df = pd.read_sql_query(species_query, conn)
    
    print(f"Found {len(species_df)} major species with significant carbon contribution:")
    for _, row in species_df.iterrows():
        print(f"  {row['COMMON_NAME']}: {row['PLOT_COUNT']} plots, {row['TOTAL_CARBON']:.1f} total carbon")
    
    # Get the species codes for our major species
    major_species_codes = species_df['SPCD'].tolist()
    
    # Main query to get species carbon data by plot
    main_query = """
    WITH latest_plots AS (
        -- Get the most recent survey year for each unique plot location
        SELECT 
            STATECD, UNITCD, COUNTYCD, PLOT,
            MAX(INVYR) as LATEST_INVYR
        FROM PLOT 
        WHERE PLOT_STATUS_CD = 1
            AND LAT IS NOT NULL 
            AND LON IS NOT NULL
        GROUP BY STATECD, UNITCD, COUNTYCD, PLOT
    )
    SELECT 
        -- Plot Information
        p.CN as PLOT_CN,
        p.STATECD,
        p.UNITCD,
        p.COUNTYCD,
        p.PLOT as PLOT_NUM,
        p.INVYR as INVENTORY_YEAR,
        p.LAT,
        p.LON,
        p.ELEV as ELEVATION_FT,
        
        -- Plot Context
        c.FORTYPCD as FOREST_TYPE_CD,
        ft.MEANING as FOREST_TYPE,
        c.OWNCD as OWNERSHIP_CD,
        c.SLOPE as SLOPE_PCT,
        c.ASPECT as ASPECT_DEG,
        c.LIVE_CANOPY_CVR_PCT as CANOPY_COVER_PCT,
        
        -- Total Plot Metrics
        COALESCE(plot_totals.TOTAL_TREES, 0) as TOTAL_TREES,
        COALESCE(plot_totals.TOTAL_CARBON_AG, 0) as TOTAL_CARBON_AG,
        COALESCE(plot_totals.TOTAL_BASAL_AREA, 0) as TOTAL_BASAL_AREA_SQFT,
        COALESCE(plot_totals.SPECIES_COUNT, 0) as SPECIES_COUNT,
        
        -- Douglas-fir (202)
        COALESCE(df.TOTAL_CARBON_AG, 0) as DOUGLAS_FIR_CARBON_AG,
        COALESCE(df.TREE_COUNT, 0) as DOUGLAS_FIR_COUNT,
        COALESCE(df.BASAL_AREA_SQFT, 0) as DOUGLAS_FIR_BASAL_AREA,
        COALESCE(df.AVG_DIAMETER, 0) as DOUGLAS_FIR_AVG_DBH,
        
        -- Western Hemlock (263)
        COALESCE(wh.TOTAL_CARBON_AG, 0) as WESTERN_HEMLOCK_CARBON_AG,
        COALESCE(wh.TREE_COUNT, 0) as WESTERN_HEMLOCK_COUNT,
        COALESCE(wh.BASAL_AREA_SQFT, 0) as WESTERN_HEMLOCK_BASAL_AREA,
        COALESCE(wh.AVG_DIAMETER, 0) as WESTERN_HEMLOCK_AVG_DBH,
        
        -- Western Red Cedar (242)
        COALESCE(wrc.TOTAL_CARBON_AG, 0) as WESTERN_RED_CEDAR_CARBON_AG,
        COALESCE(wrc.TREE_COUNT, 0) as WESTERN_RED_CEDAR_COUNT,
        COALESCE(wrc.BASAL_AREA_SQFT, 0) as WESTERN_RED_CEDAR_BASAL_AREA,
        COALESCE(wrc.AVG_DIAMETER, 0) as WESTERN_RED_CEDAR_AVG_DBH,
        
        -- Pacific Silver Fir (17)
        COALESCE(psf.TOTAL_CARBON_AG, 0) as PACIFIC_SILVER_FIR_CARBON_AG,
        COALESCE(psf.TREE_COUNT, 0) as PACIFIC_SILVER_FIR_COUNT,
        COALESCE(psf.BASAL_AREA_SQFT, 0) as PACIFIC_SILVER_FIR_BASAL_AREA,
        COALESCE(psf.AVG_DIAMETER, 0) as PACIFIC_SILVER_FIR_AVG_DBH,
        
        -- Ponderosa Pine (122)
        COALESCE(pp.TOTAL_CARBON_AG, 0) as PONDEROSA_PINE_CARBON_AG,
        COALESCE(pp.TREE_COUNT, 0) as PONDEROSA_PINE_COUNT,
        COALESCE(pp.BASAL_AREA_SQFT, 0) as PONDEROSA_PINE_BASAL_AREA,
        COALESCE(pp.AVG_DIAMETER, 0) as PONDEROSA_PINE_AVG_DBH,
        
        -- Noble Fir (11)
        COALESCE(nf.TOTAL_CARBON_AG, 0) as NOBLE_FIR_CARBON_AG,
        COALESCE(nf.TREE_COUNT, 0) as NOBLE_FIR_COUNT,
        COALESCE(nf.BASAL_AREA_SQFT, 0) as NOBLE_FIR_BASAL_AREA,
        COALESCE(nf.AVG_DIAMETER, 0) as NOBLE_FIR_AVG_DBH,
        
        -- Grand Fir (15)
        COALESCE(gf.TOTAL_CARBON_AG, 0) as GRAND_FIR_CARBON_AG,
        COALESCE(gf.TREE_COUNT, 0) as GRAND_FIR_COUNT,
        COALESCE(gf.BASAL_AREA_SQFT, 0) as GRAND_FIR_BASAL_AREA,
        COALESCE(gf.AVG_DIAMETER, 0) as GRAND_FIR_AVG_DBH,
        
        -- Sitka Spruce (98)
        COALESCE(ss.TOTAL_CARBON_AG, 0) as SITKA_SPRUCE_CARBON_AG,
        COALESCE(ss.TREE_COUNT, 0) as SITKA_SPRUCE_COUNT,
        COALESCE(ss.BASAL_AREA_SQFT, 0) as SITKA_SPRUCE_BASAL_AREA,
        COALESCE(ss.AVG_DIAMETER, 0) as SITKA_SPRUCE_AVG_DBH,
        
        -- Lodgepole Pine (108)
        COALESCE(lp.TOTAL_CARBON_AG, 0) as LODGEPOLE_PINE_CARBON_AG,
        COALESCE(lp.TREE_COUNT, 0) as LODGEPOLE_PINE_COUNT,
        COALESCE(lp.BASAL_AREA_SQFT, 0) as LODGEPOLE_PINE_BASAL_AREA,
        COALESCE(lp.AVG_DIAMETER, 0) as LODGEPOLE_PINE_AVG_DBH,
        
        -- Western Larch (73)
        COALESCE(wl.TOTAL_CARBON_AG, 0) as WESTERN_LARCH_CARBON_AG,
        COALESCE(wl.TREE_COUNT, 0) as WESTERN_LARCH_COUNT,
        COALESCE(wl.BASAL_AREA_SQFT, 0) as WESTERN_LARCH_BASAL_AREA,
        COALESCE(wl.AVG_DIAMETER, 0) as WESTERN_LARCH_AVG_DBH,
        
        -- Engelmann Spruce (93)
        COALESCE(es.TOTAL_CARBON_AG, 0) as ENGELMANN_SPRUCE_CARBON_AG,
        COALESCE(es.TREE_COUNT, 0) as ENGELMANN_SPRUCE_COUNT,
        COALESCE(es.BASAL_AREA_SQFT, 0) as ENGELMANN_SPRUCE_BASAL_AREA,
        COALESCE(es.AVG_DIAMETER, 0) as ENGELMANN_SPRUCE_AVG_DBH,
        
        -- Subalpine Fir (19)
        COALESCE(sf.TOTAL_CARBON_AG, 0) as SUBALPINE_FIR_CARBON_AG,
        COALESCE(sf.TREE_COUNT, 0) as SUBALPINE_FIR_COUNT,
        COALESCE(sf.BASAL_AREA_SQFT, 0) as SUBALPINE_FIR_BASAL_AREA,
        COALESCE(sf.AVG_DIAMETER, 0) as SUBALPINE_FIR_AVG_DBH,
        
        -- Mountain Hemlock (264)
        COALESCE(mh.TOTAL_CARBON_AG, 0) as MOUNTAIN_HEMLOCK_CARBON_AG,
        COALESCE(mh.TREE_COUNT, 0) as MOUNTAIN_HEMLOCK_COUNT,
        COALESCE(mh.BASAL_AREA_SQFT, 0) as MOUNTAIN_HEMLOCK_BASAL_AREA,
        COALESCE(mh.AVG_DIAMETER, 0) as MOUNTAIN_HEMLOCK_AVG_DBH,
        
        -- Red Alder (351) - Important hardwood
        COALESCE(ra.TOTAL_CARBON_AG, 0) as RED_ALDER_CARBON_AG,
        COALESCE(ra.TREE_COUNT, 0) as RED_ALDER_COUNT,
        COALESCE(ra.BASAL_AREA_SQFT, 0) as RED_ALDER_BASAL_AREA,
        COALESCE(ra.AVG_DIAMETER, 0) as RED_ALDER_AVG_DBH,
        
        -- Bigleaf Maple (312)
        COALESCE(bm.TOTAL_CARBON_AG, 0) as BIGLEAF_MAPLE_CARBON_AG,
        COALESCE(bm.TREE_COUNT, 0) as BIGLEAF_MAPLE_COUNT,
        COALESCE(bm.BASAL_AREA_SQFT, 0) as BIGLEAF_MAPLE_BASAL_AREA,
        COALESCE(bm.AVG_DIAMETER, 0) as BIGLEAF_MAPLE_AVG_DBH
        
    FROM PLOT p
    INNER JOIN latest_plots lp ON (
        p.STATECD = lp.STATECD AND 
        p.UNITCD = lp.UNITCD AND 
        p.COUNTYCD = lp.COUNTYCD AND 
        p.PLOT = lp.PLOT AND 
        p.INVYR = lp.LATEST_INVYR
    )
    LEFT JOIN COND c ON p.CN = c.PLT_CN AND c.CONDID = 1
    LEFT JOIN REF_FOREST_TYPE ft ON c.FORTYPCD = ft.VALUE
    
    -- Get total plot metrics
    LEFT JOIN (
        SELECT 
            PLT_CN,
            COUNT(*) as TOTAL_TREES,
            SUM(CARBON_AG) as TOTAL_CARBON_AG,
            SUM(DIA * DIA * 0.005454) as TOTAL_BASAL_AREA,
            COUNT(DISTINCT SPCD) as SPECIES_COUNT
        FROM TREE
        WHERE STATUSCD = 1 AND CARBON_AG IS NOT NULL AND CARBON_AG > 0
        GROUP BY PLT_CN
    ) plot_totals ON p.CN = plot_totals.PLT_CN
    
    -- Join individual species data
    LEFT JOIN (
        SELECT 
            PLT_CN,
            COUNT(*) as TREE_COUNT,
            SUM(CARBON_AG) as TOTAL_CARBON_AG,
            AVG(DIA) as AVG_DIAMETER,
            SUM(DIA * DIA * 0.005454) as BASAL_AREA_SQFT
        FROM TREE
        WHERE STATUSCD = 1 AND CARBON_AG > 0 AND SPCD = 202
        GROUP BY PLT_CN
    ) df ON p.CN = df.PLT_CN
    
    LEFT JOIN (
        SELECT 
            PLT_CN,
            COUNT(*) as TREE_COUNT,
            SUM(CARBON_AG) as TOTAL_CARBON_AG,
            AVG(DIA) as AVG_DIAMETER,
            SUM(DIA * DIA * 0.005454) as BASAL_AREA_SQFT
        FROM TREE
        WHERE STATUSCD = 1 AND CARBON_AG > 0 AND SPCD = 263
        GROUP BY PLT_CN
    ) wh ON p.CN = wh.PLT_CN
    
    LEFT JOIN (
        SELECT 
            PLT_CN,
            COUNT(*) as TREE_COUNT,
            SUM(CARBON_AG) as TOTAL_CARBON_AG,
            AVG(DIA) as AVG_DIAMETER,
            SUM(DIA * DIA * 0.005454) as BASAL_AREA_SQFT
        FROM TREE
        WHERE STATUSCD = 1 AND CARBON_AG > 0 AND SPCD = 242
        GROUP BY PLT_CN
    ) wrc ON p.CN = wrc.PLT_CN
    
    LEFT JOIN (
        SELECT 
            PLT_CN,
            COUNT(*) as TREE_COUNT,
            SUM(CARBON_AG) as TOTAL_CARBON_AG,
            AVG(DIA) as AVG_DIAMETER,
            SUM(DIA * DIA * 0.005454) as BASAL_AREA_SQFT
        FROM TREE
        WHERE STATUSCD = 1 AND CARBON_AG > 0 AND SPCD = 17
        GROUP BY PLT_CN
    ) psf ON p.CN = psf.PLT_CN
    
    LEFT JOIN (
        SELECT 
            PLT_CN,
            COUNT(*) as TREE_COUNT,
            SUM(CARBON_AG) as TOTAL_CARBON_AG,
            AVG(DIA) as AVG_DIAMETER,
            SUM(DIA * DIA * 0.005454) as BASAL_AREA_SQFT
        FROM TREE
        WHERE STATUSCD = 1 AND CARBON_AG > 0 AND SPCD = 122
        GROUP BY PLT_CN
    ) pp ON p.CN = pp.PLT_CN
    
    LEFT JOIN (
        SELECT 
            PLT_CN,
            COUNT(*) as TREE_COUNT,
            SUM(CARBON_AG) as TOTAL_CARBON_AG,
            AVG(DIA) as AVG_DIAMETER,
            SUM(DIA * DIA * 0.005454) as BASAL_AREA_SQFT
        FROM TREE
        WHERE STATUSCD = 1 AND CARBON_AG > 0 AND SPCD = 11
        GROUP BY PLT_CN
    ) nf ON p.CN = nf.PLT_CN
    
    LEFT JOIN (
        SELECT 
            PLT_CN,
            COUNT(*) as TREE_COUNT,
            SUM(CARBON_AG) as TOTAL_CARBON_AG,
            AVG(DIA) as AVG_DIAMETER,
            SUM(DIA * DIA * 0.005454) as BASAL_AREA_SQFT
        FROM TREE
        WHERE STATUSCD = 1 AND CARBON_AG > 0 AND SPCD = 15
        GROUP BY PLT_CN
    ) gf ON p.CN = gf.PLT_CN
    
    LEFT JOIN (
        SELECT 
            PLT_CN,
            COUNT(*) as TREE_COUNT,
            SUM(CARBON_AG) as TOTAL_CARBON_AG,
            AVG(DIA) as AVG_DIAMETER,
            SUM(DIA * DIA * 0.005454) as BASAL_AREA_SQFT
        FROM TREE
        WHERE STATUSCD = 1 AND CARBON_AG > 0 AND SPCD = 98
        GROUP BY PLT_CN
    ) ss ON p.CN = ss.PLT_CN
    
    LEFT JOIN (
        SELECT 
            PLT_CN,
            COUNT(*) as TREE_COUNT,
            SUM(CARBON_AG) as TOTAL_CARBON_AG,
            AVG(DIA) as AVG_DIAMETER,
            SUM(DIA * DIA * 0.005454) as BASAL_AREA_SQFT
        FROM TREE
        WHERE STATUSCD = 1 AND CARBON_AG > 0 AND SPCD = 108
        GROUP BY PLT_CN
    ) lp ON p.CN = lp.PLT_CN
    
    LEFT JOIN (
        SELECT 
            PLT_CN,
            COUNT(*) as TREE_COUNT,
            SUM(CARBON_AG) as TOTAL_CARBON_AG,
            AVG(DIA) as AVG_DIAMETER,
            SUM(DIA * DIA * 0.005454) as BASAL_AREA_SQFT
        FROM TREE
        WHERE STATUSCD = 1 AND CARBON_AG > 0 AND SPCD = 73
        GROUP BY PLT_CN
    ) wl ON p.CN = wl.PLT_CN
    
    LEFT JOIN (
        SELECT 
            PLT_CN,
            COUNT(*) as TREE_COUNT,
            SUM(CARBON_AG) as TOTAL_CARBON_AG,
            AVG(DIA) as AVG_DIAMETER,
            SUM(DIA * DIA * 0.005454) as BASAL_AREA_SQFT
        FROM TREE
        WHERE STATUSCD = 1 AND CARBON_AG > 0 AND SPCD = 93
        GROUP BY PLT_CN
    ) es ON p.CN = es.PLT_CN
    
    LEFT JOIN (
        SELECT 
            PLT_CN,
            COUNT(*) as TREE_COUNT,
            SUM(CARBON_AG) as TOTAL_CARBON_AG,
            AVG(DIA) as AVG_DIAMETER,
            SUM(DIA * DIA * 0.005454) as BASAL_AREA_SQFT
        FROM TREE
        WHERE STATUSCD = 1 AND CARBON_AG > 0 AND SPCD = 19
        GROUP BY PLT_CN
    ) sf ON p.CN = sf.PLT_CN
    
    LEFT JOIN (
        SELECT 
            PLT_CN,
            COUNT(*) as TREE_COUNT,
            SUM(CARBON_AG) as TOTAL_CARBON_AG,
            AVG(DIA) as AVG_DIAMETER,
            SUM(DIA * DIA * 0.005454) as BASAL_AREA_SQFT
        FROM TREE
        WHERE STATUSCD = 1 AND CARBON_AG > 0 AND SPCD = 264
        GROUP BY PLT_CN
    ) mh ON p.CN = mh.PLT_CN
    
    LEFT JOIN (
        SELECT 
            PLT_CN,
            COUNT(*) as TREE_COUNT,
            SUM(CARBON_AG) as TOTAL_CARBON_AG,
            AVG(DIA) as AVG_DIAMETER,
            SUM(DIA * DIA * 0.005454) as BASAL_AREA_SQFT
        FROM TREE
        WHERE STATUSCD = 1 AND CARBON_AG > 0 AND SPCD = 351
        GROUP BY PLT_CN
    ) ra ON p.CN = ra.PLT_CN
    
    LEFT JOIN (
        SELECT 
            PLT_CN,
            COUNT(*) as TREE_COUNT,
            SUM(CARBON_AG) as TOTAL_CARBON_AG,
            AVG(DIA) as AVG_DIAMETER,
            SUM(DIA * DIA * 0.005454) as BASAL_AREA_SQFT
        FROM TREE
        WHERE STATUSCD = 1 AND CARBON_AG > 0 AND SPCD = 312
        GROUP BY PLT_CN
    ) bm ON p.CN = bm.PLT_CN
    
    WHERE p.LAT IS NOT NULL 
        AND p.LON IS NOT NULL
        AND p.PLOT_STATUS_CD = 1
        AND plot_totals.TOTAL_TREES > 0  -- Only plots with trees
    ORDER BY plot_totals.TOTAL_CARBON_AG DESC;
    """
    
    print("Executing species carbon query (this may take a few minutes)...")
    df = pd.read_sql_query(main_query, conn)
    
    # Get some statistics about the duplicate removal
    duplicate_check_query = """
    SELECT 
        COUNT(*) as total_plot_records,
        COUNT(DISTINCT STATECD || '-' || UNITCD || '-' || COUNTYCD || '-' || PLOT) as unique_plot_locations,
        COUNT(*) - COUNT(DISTINCT STATECD || '-' || UNITCD || '-' || COUNTYCD || '-' || PLOT) as duplicates_removed
    FROM PLOT 
    WHERE PLOT_STATUS_CD = 1 AND LAT IS NOT NULL AND LON IS NOT NULL;
    """
    dup_stats = pd.read_sql_query(duplicate_check_query, conn)
    
    conn.close()
    
    print(f"Extracted species carbon data for {len(df):,} plots")
    print(f"  📊 Total plot records in database: {dup_stats.iloc[0]['total_plot_records']:,}")
    print(f"  📍 Unique plot locations: {dup_stats.iloc[0]['unique_plot_locations']:,}")
    print(f"  🔄 Duplicate surveys removed: {dup_stats.iloc[0]['duplicates_removed']:,}")
    print(f"  ✅ Keeping only latest survey per plot location")
    
    # Create GeoDataFrame
    geometry = [Point(xy) for xy in zip(df['LON'], df['LAT'])]
    gdf = gpd.GeoDataFrame(df, geometry=geometry, crs='EPSG:4326')
    
    # Add derived metrics for visualization
    gdf['COORDINATES_ACCURACY'] = 'FUZZED_0.5_TO_1.0_MILE'
    
    # Calculate diversity metrics
    species_carbon_cols = [col for col in gdf.columns if col.endswith('_CARBON_AG') and col != 'TOTAL_CARBON_AG']
    
    # Species richness (number of species with carbon > 0)
    gdf['SPECIES_RICHNESS'] = (gdf[species_carbon_cols] > 0).sum(axis=1)
    
    # Dominant species (species with highest carbon)
    def get_dominant_species(row):
        species_values = {}
        for col in species_carbon_cols:
            species_name = col.replace('_CARBON_AG', '').replace('_', ' ').title()
            species_values[species_name] = row[col]
        if species_values:
            return max(species_values.items(), key=lambda x: x[1])[0]
        return 'None'
    
    gdf['DOMINANT_SPECIES'] = gdf.apply(get_dominant_species, axis=1)
    
    # Calculate diversity index (Shannon-Weaver)
    def calculate_shannon_diversity(row):
        values = [row[col] for col in species_carbon_cols if row[col] > 0]
        if len(values) <= 1:
            return 0
        total = sum(values)
        if total == 0:
            return 0
        proportions = [v/total for v in values]
        return -sum(p * np.log(p) for p in proportions if p > 0)
    
    gdf['SHANNON_DIVERSITY'] = gdf.apply(calculate_shannon_diversity, axis=1)
    
    # Calculate relative abundance percentages for major species
    for col in species_carbon_cols:
        pct_col = col.replace('_CARBON_AG', '_PERCENT')
        gdf[pct_col] = (gdf[col] / gdf['TOTAL_CARBON_AG'] * 100).fillna(0).round(1)
    
    # Add ownership type labels
    gdf['OWNERSHIP_TYPE'] = gdf['OWNERSHIP_CD'].map({
        11: 'National Forest Service',
        21: 'Other Federal',
        31: 'State and Local Government',
        32: 'Private Corporation',
        46: 'Private Individual'
    }).fillna('Unknown')
    
    # Export to GeoJSON
    output_file = OUTPUT_DIR / "fia_species_carbon_wa_complete.geojson"
    gdf.to_file(output_file, driver='GeoJSON')
    print(f"Species carbon layer saved to: {output_file}")
    
    # Create a summary report
    print(f"\n" + "="*60)
    print(f"SPECIES CARBON LAYER SUMMARY")
    print(f"="*60)
    print(f"  Total plots with species data: {len(gdf):,}")
    print(f"  Average species richness per plot: {gdf['SPECIES_RICHNESS'].mean():.1f}")
    print(f"  Average Shannon diversity: {gdf['SHANNON_DIVERSITY'].mean():.2f}")
    print(f"  Total carbon range: {gdf['TOTAL_CARBON_AG'].min():.1f} to {gdf['TOTAL_CARBON_AG'].max():.1f}")
    
    print(f"\nTop dominant species:")
    dominant_counts = gdf['DOMINANT_SPECIES'].value_counts().head(10)
    for species, count in dominant_counts.items():
        pct = 100 * count / len(gdf)
        print(f"  {species}: {count:,} plots ({pct:.1f}%)")
    
    print(f"\nSpecies carbon statistics (top 5 species):")
    species_stats = []
    for col in species_carbon_cols:
        species_name = col.replace('_CARBON_AG', '').replace('_', ' ').title()
        plot_count = (gdf[col] > 0).sum()
        if plot_count > 0:
            mean_carbon = gdf[col][gdf[col] > 0].mean()
            max_carbon = gdf[col].max()
            total_carbon = gdf[col].sum()
            species_stats.append({
                'species': species_name,
                'plots': plot_count,
                'mean_carbon': mean_carbon,
                'max_carbon': max_carbon,
                'total_carbon': total_carbon
            })
    
    species_stats.sort(key=lambda x: x['total_carbon'], reverse=True)
    for stat in species_stats[:8]:
        print(f"  {stat['species']}:")
        print(f"    Present in {stat['plots']} plots")
        print(f"    Mean carbon (when present): {stat['mean_carbon']:.1f}")
        print(f"    Max carbon: {stat['max_carbon']:.1f}")
        print(f"    Total carbon: {stat['total_carbon']:.1f}")
    
    print(f"\nQGIS USAGE INSTRUCTIONS:")
    print(f"  1. Load 'fia_species_carbon_wa_complete.geojson' in QGIS")
    print(f"  2. Style by dominant species using 'DOMINANT_SPECIES' field")
    print(f"  3. Use graduated symbols for individual species carbon (e.g., 'DOUGLAS_FIR_CARBON_AG')")
    print(f"  4. Filter by 'OWNERSHIP_TYPE' or elevation for targeted analysis")
    print(f"  5. Use field calculator to compute species ratios or diversity metrics")
    print(f"  6. Create choropleth maps using percentage fields (e.g., 'DOUGLAS_FIR_PERCENT')")
    
    print(f"\nKEY FIELDS FOR ANALYSIS:")
    print(f"  • Species carbon: *_CARBON_AG (total carbon per species)")
    print(f"  • Species count: *_COUNT (number of trees per species)")
    print(f"  • Species percentage: *_PERCENT (relative abundance %)")
    print(f"  • DOMINANT_SPECIES: Species with highest carbon")
    print(f"  • SPECIES_RICHNESS: Number of species present")
    print(f"  • SHANNON_DIVERSITY: Diversity index")
    
    return gdf

def main():
    """
    Main execution function - creates species carbon layer
    """
    try:
        # Create species carbon layer
        species_gdf = create_species_carbon_layer()
        
        print(f"\n✅ SUCCESS! Created species carbon layer for QGIS")
        print(f"\nFile created:")
        print(f"  📊 outputs/fia_species_carbon_wa_complete.geojson")
        
        print(f"\nNext steps:")
        print(f"  1. Open QGIS and load the GeoJSON file")
        print(f"  2. Explore species distribution and abundance patterns")
        print(f"  3. Use for species composition analysis and habitat modeling")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main() 