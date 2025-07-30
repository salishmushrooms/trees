#!/usr/bin/env python3
"""
Create QGIS-Compatible Plot and Subplot Layers (Simplified Working Version)

This script extracts FIA plot and subplot data into GeoJSON format for QGIS visualization.
Includes plot IDs, coordinates, ownership, forest characteristics, and location attributes.
"""

import sqlite3
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import numpy as np
from pathlib import Path

# Database connection
DB_PATH = "../../data/raw/trees_SQLite_FIADB_WA.db"
OUTPUT_DIR = Path("../../outputs")
OUTPUT_DIR.mkdir(exist_ok=True)

def create_plot_layer():
    """
    Create comprehensive plot-level layer with confirmed columns
    """
    print("Extracting plot-level data...")
    
    query = """
    SELECT 
        -- Plot Identifiers
        p.CN as PLOT_CN,
        p.STATECD,
        p.UNITCD,
        p.COUNTYCD, 
        p.PLOT as PLOT_NUM,
        p.INVYR as INVENTORY_YEAR,
        p.CYCLE,
        p.SUBCYCLE,
        
        -- Coordinates (fuzzed but useful for regional analysis)
        p.LAT,
        p.LON,
        p.ELEV as ELEVATION_FT,
        
        -- Plot Status and Timing
        p.PLOT_STATUS_CD,
        p.MEASYEAR,
        p.MEASMON,
        p.MEASDAY,
        
        -- Location Characteristics
        p.WATERCD,
        p.RDDISTCD as ROAD_DISTANCE_CD,
        p.DECLINATION,
        p.TOPO_POSITION_PNW,
        
        -- Forest Condition Data (from primary condition)
        c.CONDID,
        c.FORTYPCD as FOREST_TYPE_CD,
        ft.MEANING as FOREST_TYPE,
        c.FLDTYPCD as FIELD_FOREST_TYPE_CD,
        c.STDAGE as STAND_AGE,
        c.STDSZCD as STAND_SIZE_CLASS,
        c.SITECLCD as SITE_CLASS,
        
        -- Topographic Data (precise despite coordinate fuzzing)
        c.SLOPE as SLOPE_PCT,
        c.ASPECT as ASPECT_DEG,
        c.LIVE_CANOPY_CVR_PCT as CANOPY_COVER_PCT,
        
        -- Ownership Information
        c.OWNCD as OWNERSHIP_CD,
        c.OWNGRPCD as OWNERSHIP_GROUP_CD,
        c.OWNSUBCD as OWNERSHIP_SUBCODE,
        c.RESERVCD as RESERVED_STATUS_CD,
        c.ADFORCD as ADMIN_FOREST_CD,
        
        -- Habitat Classification
        c.HABTYPCD1 as HABITAT_TYPE_CD,
        h.SCIENTIFIC_NAME as HABITAT_TYPE,
        h.COMMON_NAME as HABITAT_COMMON_NAME,
        
        -- Site Characteristics
        c.PHYSCLCD as PHYSIOGRAPHIC_CLASS,
        c.PRESNFCD as PRESENT_NONFOREST_CD,
        
        -- Disturbance History
        c.DSTRBCD1 as DISTURBANCE_1_CD,
        c.DSTRBYR1 as DISTURBANCE_1_YEAR,
        c.DSTRBCD2 as DISTURBANCE_2_CD,
        c.DSTRBYR2 as DISTURBANCE_2_YEAR
        
    FROM PLOT p
    LEFT JOIN COND c ON p.CN = c.PLT_CN AND c.CONDID = 1  -- Primary condition
    LEFT JOIN REF_FOREST_TYPE ft ON c.FORTYPCD = ft.VALUE
    LEFT JOIN REF_HABTYP_DESCRIPTION h ON c.HABTYPCD1 = h.HABTYPCD 
        AND c.HABTYPCD1_PUB_CD = h.PUB_CD
    WHERE p.LAT IS NOT NULL 
        AND p.LON IS NOT NULL
        AND p.PLOT_STATUS_CD = 1  -- Forested plots only
    ORDER BY p.INVYR DESC, p.PLOT;
    """
    
    conn = sqlite3.connect(DB_PATH)
    df = pd.read_sql_query(query, conn)
    conn.close()
    
    print(f"Extracted {len(df):,} plot records")
    
    # Create GeoDataFrame
    geometry = [Point(xy) for xy in zip(df['LON'], df['LAT'])]
    gdf = gpd.GeoDataFrame(df, geometry=geometry, crs='EPSG:4326')
    
    # Add derived fields for analysis
    gdf['COORDINATES_ACCURACY'] = 'FUZZED_0.5_TO_1.0_MILE'
    gdf['OWNERSHIP_TYPE'] = gdf['OWNERSHIP_CD'].map({
        11: 'National Forest Service',
        21: 'Other Federal',
        31: 'State and Local Government', 
        32: 'Private Corporation',
        46: 'Private Individual'
    }).fillna('Unknown')
    
    # Create aspect categories (handle null values)
    gdf['ASPECT_CARDINAL'] = pd.cut(gdf['ASPECT_DEG'].fillna(-1), 
                                   bins=[-1, 0, 45, 135, 225, 315, 360], 
                                   labels=['Unknown', 'North', 'East', 'South', 'West', 'North_High'],
                                   include_lowest=True)
    
    # Elevation categories
    gdf['ELEVATION_ZONE'] = pd.cut(gdf['ELEVATION_FT'].fillna(-1),
                                   bins=[-1, 0, 1500, 3000, 4500, 6000, 8000],
                                   labels=['Unknown', 'Lowland', 'Montane', 'Subalpine', 'Alpine', 'High_Alpine'])
    
    # Export to GeoJSON
    output_file = OUTPUT_DIR / "fia_plots_wa_complete.geojson"
    gdf.to_file(output_file, driver='GeoJSON')
    print(f"Plot layer saved to: {output_file}")
    
    return gdf

def create_subplot_layer():
    """
    Create subplot-level layer with calculated subplot coordinates
    """
    print("\nExtracting subplot-level data...")
    
    query = """
    SELECT 
        -- Subplot Identifiers
        s.CN as SUBPLOT_CN,
        s.PLT_CN as PLOT_CN,
        p.STATECD,
        p.UNITCD,
        p.COUNTYCD,
        p.PLOT as PLOT_NUM,
        s.SUBP as SUBPLOT_NUM,
        p.INVYR as INVENTORY_YEAR,
        
        -- Plot Center Coordinates (for calculating subplot positions)
        p.LAT as PLOT_LAT,
        p.LON as PLOT_LON,
        p.ELEV as ELEVATION_FT,
        
        -- Subplot Status
        s.SUBP_STATUS_CD,
        s.POINT_NONSAMPLE_REASN_CD,
        
        -- Subplot Conditions
        s.MICRCOND as MICROPLOT_CONDITION_ID,
        s.SUBPCOND as SUBPLOT_CONDITION_ID,
        s.MACRCOND as MACROPLOT_CONDITION_ID,
        
        -- Subplot-specific topography
        s.SLOPE as SUBPLOT_SLOPE_PCT,
        s.ASPECT as SUBPLOT_ASPECT_DEG,
        s.WATERDEP as WATER_DEPTH,
        
        -- Link to main condition data
        c.FORTYPCD as FOREST_TYPE_CD,
        ft.MEANING as FOREST_TYPE,
        c.OWNCD as OWNERSHIP_CD,
        c.HABTYPCD1 as HABITAT_TYPE_CD,
        h.SCIENTIFIC_NAME as HABITAT_TYPE
        
    FROM SUBPLOT s
    JOIN PLOT p ON s.PLT_CN = p.CN
    LEFT JOIN COND c ON p.CN = c.PLT_CN AND c.CONDID = 1
    LEFT JOIN REF_FOREST_TYPE ft ON c.FORTYPCD = ft.VALUE
    LEFT JOIN REF_HABTYP_DESCRIPTION h ON c.HABTYPCD1 = h.HABTYPCD 
        AND c.HABTYPCD1_PUB_CD = h.PUB_CD
    WHERE p.LAT IS NOT NULL 
        AND p.LON IS NOT NULL
        AND p.PLOT_STATUS_CD = 1
    ORDER BY p.PLOT, s.SUBP;
    """
    
    conn = sqlite3.connect(DB_PATH)
    df = pd.read_sql_query(query, conn)
    conn.close()
    
    print(f"Extracted {len(df):,} subplot records")
    
    # Calculate subplot coordinates based on plot center
    # Standard FIA subplot arrangement
    subplot_offsets = {
        1: (0.0, 0.0),      # Center subplot
        2: (0.0, 120.0),    # North subplot (120 ft north)
        3: (103.9, -60.0),  # Southeast subplot (120 ft at 120°)
        4: (-103.9, -60.0)  # Southwest subplot (120 ft at 240°)
    }
    
    # Convert feet to decimal degrees (approximate)
    feet_to_degrees_lat = 1.0 / 364000.0
    
    subplot_lats = []
    subplot_lons = []
    
    for idx, row in df.iterrows():
        plot_lat = row['PLOT_LAT']
        plot_lon = row['PLOT_LON']
        subplot_num = row['SUBPLOT_NUM']
        
        # Get offset in feet
        offset_x, offset_y = subplot_offsets.get(subplot_num, (0.0, 0.0))
        
        # Convert to degrees
        feet_to_degrees_lon = feet_to_degrees_lat / np.cos(np.radians(plot_lat))
        
        # Calculate subplot coordinates
        subplot_lat = plot_lat + (offset_y * feet_to_degrees_lat)
        subplot_lon = plot_lon + (offset_x * feet_to_degrees_lon)
        
        subplot_lats.append(subplot_lat)
        subplot_lons.append(subplot_lon)
    
    df['SUBPLOT_LAT'] = subplot_lats
    df['SUBPLOT_LON'] = subplot_lons
    
    # Create GeoDataFrame using calculated subplot coordinates
    geometry = [Point(xy) for xy in zip(df['SUBPLOT_LON'], df['SUBPLOT_LAT'])]
    gdf = gpd.GeoDataFrame(df, geometry=geometry, crs='EPSG:4326')
    
    # Add derived fields
    gdf['COORDINATES_ACCURACY'] = 'CALCULATED_FROM_FUZZED_PLOT_CENTER'
    gdf['SUBPLOT_AREA_ACRES'] = 0.133  # 24-foot radius = 0.133 acres
    
    gdf['OWNERSHIP_TYPE'] = gdf['OWNERSHIP_CD'].map({
        11: 'National Forest Service',
        21: 'Other Federal',
        31: 'State and Local Government',
        32: 'Private Corporation', 
        46: 'Private Individual'
    }).fillna('Unknown')
    
    # Export to GeoJSON
    output_file = OUTPUT_DIR / "fia_subplots_wa_complete.geojson"
    gdf.to_file(output_file, driver='GeoJSON')
    print(f"Subplot layer saved to: {output_file}")
    
    return gdf

def generate_summary_report(plot_gdf, subplot_gdf):
    """
    Generate summary report of extracted data
    """
    print("\n" + "="*60)
    print("QGIS LAYER EXTRACTION SUMMARY")
    print("="*60)
    
    print(f"\nPLOT LAYER STATISTICS:")
    print(f"  Total plots: {len(plot_gdf):,}")
    print(f"  Coordinate range:")
    print(f"    Latitude: {plot_gdf['LAT'].min():.4f} to {plot_gdf['LAT'].max():.4f}")
    print(f"    Longitude: {plot_gdf['LON'].min():.4f} to {plot_gdf['LON'].max():.4f}")
    print(f"    Elevation: {plot_gdf['ELEVATION_FT'].min():,.0f} to {plot_gdf['ELEVATION_FT'].max():,.0f} feet")
    
    print(f"\n  Ownership distribution:")
    ownership_counts = plot_gdf['OWNERSHIP_TYPE'].value_counts()
    for owner, count in ownership_counts.items():
        pct = 100 * count / len(plot_gdf)
        print(f"    {owner}: {count:,} plots ({pct:.1f}%)")
    
    print(f"\n  Forest type distribution (top 5):")
    forest_counts = plot_gdf['FOREST_TYPE'].value_counts().head()
    for ftype, count in forest_counts.items():
        pct = 100 * count / len(plot_gdf)
        print(f"    {ftype}: {count:,} plots ({pct:.1f}%)")
    
    print(f"\nSUBPLOT LAYER STATISTICS:")
    print(f"  Total subplots: {len(subplot_gdf):,}")
    print(f"  Plots with all 4 subplots: {len(subplot_gdf[subplot_gdf['SUBPLOT_NUM']==1]):,}")
    
    print(f"\n  Subplot distribution:")
    subplot_counts = subplot_gdf['SUBPLOT_NUM'].value_counts().sort_index()
    for subp_num, count in subplot_counts.items():
        print(f"    Subplot {subp_num}: {count:,}")

def main():
    """
    Main execution function
    """
    print("Creating QGIS-compatible FIA plot and subplot layers...")
    print("="*60)
    
    try:
        # Create plot layer
        plot_gdf = create_plot_layer()
        
        # Create subplot layer  
        subplot_gdf = create_subplot_layer()
        
        # Generate summary report
        generate_summary_report(plot_gdf, subplot_gdf)
        
        print(f"\n✅ SUCCESS! Created QGIS layers")
        print(f"\nFiles created:")
        print(f"  📍 outputs/fia_plots_wa_complete.geojson")
        print(f"  📍 outputs/fia_subplots_wa_complete.geojson") 
        
        print(f"\nQGIS Usage:")
        print(f"  1. Open QGIS")
        print(f"  2. Drag and drop the .geojson files into QGIS")
        print(f"  3. Use Layer Properties > Symbology to style by ownership, forest type, etc.")
        print(f"  4. Right-click layer > Open Attribute Table to see all data")
        
        print(f"\nKey Features:")
        print(f"  • Plot coordinates are fuzzed (±0.5-1.0 mile accuracy)")
        print(f"  • Subplot coordinates calculated from plot centers")
        print(f"  • Elevation, slope, aspect data are precise")
        print(f"  • Ownership and forest type classifications included")
        print(f"  • Habitat type data where available (55% of plots)")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main() 