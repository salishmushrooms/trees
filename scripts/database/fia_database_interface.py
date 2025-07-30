#!/usr/bin/env python3
"""
FIA Database Interface Module

Reusable functions for querying the FIA database with common patterns.
Use this module for consistent database access across different scripts.

Example usage:
    from scripts.database.fia_database_interface import *
    
    # Get Douglas-fir plots above 2000 feet
    douglas_fir_plots = query_plots_by_criteria(
        forest_types=[201], 
        elevation_range=(2000, 8000)
    )
    
    # Export to QGIS
    export_to_qgis(douglas_fir_plots, "douglas_fir_high_elevation")
"""

import sqlite3
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
from pathlib import Path

DB_PATH = "../../data/raw/trees_SQLite_FIADB_WA.db"

def get_database_connection():
    """Get database connection"""
    return sqlite3.connect(DB_PATH)

def query_plots_by_criteria(forest_types=None, ownership_codes=None, 
                           elevation_range=None, inventory_years=None,
                           slope_range=None, aspect_range=None):
    """
    Query plots by various criteria
    
    Args:
        forest_types: List of forest type codes (e.g., [201, 301])
        ownership_codes: List of ownership codes (e.g., [11, 46])
        elevation_range: Tuple of (min_elev, max_elev) in feet
        inventory_years: List of inventory years (e.g., [2017, 2018])
        slope_range: Tuple of (min_slope, max_slope) in percent
        aspect_range: Tuple of (min_aspect, max_aspect) in degrees
    
    Returns:
        GeoDataFrame with matching plots
    """
    conditions = ["p.PLOT_STATUS_CD = 1"]  # Forested plots
    
    if forest_types:
        conditions.append(f"c.FORTYPCD IN ({','.join(map(str, forest_types))})")
    
    if ownership_codes:
        conditions.append(f"c.OWNCD IN ({','.join(map(str, ownership_codes))})")
    
    if elevation_range:
        conditions.append(f"p.ELEV BETWEEN {elevation_range[0]} AND {elevation_range[1]}")
    
    if inventory_years:
        conditions.append(f"p.INVYR IN ({','.join(map(str, inventory_years))})")
    
    if slope_range:
        conditions.append(f"c.SLOPE BETWEEN {slope_range[0]} AND {slope_range[1]}")
    
    if aspect_range:
        conditions.append(f"c.ASPECT BETWEEN {aspect_range[0]} AND {aspect_range[1]}")
    
    where_clause = " AND ".join(conditions)
    
    query = f"""
    SELECT 
        p.CN as PLOT_CN, p.PLOT as PLOT_NUM, p.LAT, p.LON, p.ELEV,
        p.INVYR as INVENTORY_YEAR, c.FORTYPCD, ft.MEANING as FOREST_TYPE,
        c.OWNCD as OWNERSHIP_CD, c.SLOPE, c.ASPECT, c.HABTYPCD1,
        h.SCIENTIFIC_NAME as HABITAT_TYPE,
        c.LIVE_CANOPY_CVR_PCT as CANOPY_COVER
    FROM PLOT p
    LEFT JOIN COND c ON p.CN = c.PLT_CN AND c.CONDID = 1
    LEFT JOIN REF_FOREST_TYPE ft ON c.FORTYPCD = ft.VALUE
    LEFT JOIN REF_HABTYP_DESCRIPTION h ON c.HABTYPCD1 = h.HABTYPCD 
        AND c.HABTYPCD1_PUB_CD = h.PUB_CD
    WHERE {where_clause}
    ORDER BY p.PLOT;
    """
    
    conn = get_database_connection()
    df = pd.read_sql_query(query, conn)
    conn.close()
    
    # Convert to GeoDataFrame
    geometry = [Point(xy) for xy in zip(df['LON'], df['LAT'])]
    gdf = gpd.GeoDataFrame(df, geometry=geometry, crs='EPSG:4326')
    
    # Add ownership type labels
    gdf['OWNERSHIP_TYPE'] = gdf['OWNERSHIP_CD'].map({
        11: 'National Forest Service',
        21: 'Other Federal',
        31: 'State and Local Government',
        32: 'Private Corporation', 
        46: 'Private Individual'
    }).fillna('Unknown')
    
    return gdf

def query_trees_by_plot(plot_cns, species_codes=None, min_diameter=None):
    """
    Get tree data for specific plots
    
    Args:
        plot_cns: List of plot CN identifiers
        species_codes: Optional list of species codes to filter
        min_diameter: Minimum diameter at breast height (inches)
    
    Returns:
        DataFrame with tree data
    """
    if isinstance(plot_cns, str):
        plot_cns = [plot_cns]
    
    plot_cn_str = "','".join(plot_cns)
    conditions = [f"t.PLT_CN IN ('{plot_cn_str}')"]
    
    if species_codes:
        conditions.append(f"t.SPCD IN ({','.join(map(str, species_codes))})")
    
    if min_diameter:
        conditions.append(f"t.DIA >= {min_diameter}")
    
    where_clause = " AND ".join(conditions)
    
    query = f"""
    SELECT 
        t.PLT_CN, t.SUBP, t.SPCD, s.COMMON_NAME as SPECIES,
        t.DIA, t.HT, t.CARBON_AG, t.STATUSCD,
        p.PLOT as PLOT_NUM, p.LAT, p.LON
    FROM TREE t
    JOIN REF_SPECIES s ON t.SPCD = s.SPCD
    JOIN PLOT p ON t.PLT_CN = p.CN
    WHERE {where_clause}
    ORDER BY t.PLT_CN, t.SUBP, t.SPCD;
    """
    
    conn = get_database_connection()
    df = pd.read_sql_query(query, conn)
    conn.close()
    
    return df

def get_habitat_types_by_region(forest_type_codes=None, ownership_codes=None):
    """
    Get habitat type distributions by forest type
    
    Args:
        forest_type_codes: List of forest type codes to filter
        ownership_codes: List of ownership codes to filter
    
    Returns:
        DataFrame with habitat type summaries
    """
    conditions = ["c.HABTYPCD1 IS NOT NULL"]
    
    if forest_type_codes:
        conditions.append(f"c.FORTYPCD IN ({','.join(map(str, forest_type_codes))})")
    
    if ownership_codes:
        conditions.append(f"c.OWNCD IN ({','.join(map(str, ownership_codes))})")
    
    where_clause = " AND ".join(conditions)
    
    query = f"""
    SELECT 
        c.FORTYPCD, ft.MEANING as FOREST_TYPE,
        c.HABTYPCD1, h.SCIENTIFIC_NAME as HABITAT_TYPE,
        h.COMMON_NAME as HABITAT_COMMON_NAME,
        c.OWNCD as OWNERSHIP_CD,
        COUNT(*) as PLOT_COUNT,
        AVG(p.ELEV) as AVG_ELEVATION,
        AVG(c.SLOPE) as AVG_SLOPE
    FROM COND c
    JOIN PLOT p ON c.PLT_CN = p.CN
    JOIN REF_FOREST_TYPE ft ON c.FORTYPCD = ft.VALUE
    JOIN REF_HABTYP_DESCRIPTION h ON c.HABTYPCD1 = h.HABTYPCD 
        AND c.HABTYPCD1_PUB_CD = h.PUB_CD
    WHERE {where_clause}
    GROUP BY c.FORTYPCD, c.HABTYPCD1, c.OWNCD
    ORDER BY PLOT_COUNT DESC;
    """
    
    conn = get_database_connection()
    df = pd.read_sql_query(query, conn)
    conn.close()
    
    return df

def get_species_composition_by_plot(plot_cns):
    """
    Get detailed species composition for specific plots
    
    Args:
        plot_cns: List of plot CN identifiers
    
    Returns:
        DataFrame with species composition summary
    """
    if isinstance(plot_cns, str):
        plot_cns = [plot_cns]
    
    plot_cn_str = "','".join(plot_cns)
    
    query = f"""
    SELECT 
        t.PLT_CN, p.PLOT as PLOT_NUM, p.LAT, p.LON, p.ELEV,
        t.SPCD, s.COMMON_NAME as SPECIES,
        COUNT(*) as TREE_COUNT,
        AVG(t.DIA) as AVG_DIAMETER,
        SUM(t.CARBON_AG) as TOTAL_CARBON,
        MAX(t.DIA) as MAX_DIAMETER
    FROM TREE t
    JOIN REF_SPECIES s ON t.SPCD = s.SPCD
    JOIN PLOT p ON t.PLT_CN = p.CN
    WHERE t.PLT_CN IN ('{plot_cn_str}')
        AND t.STATUSCD = 1  -- Live trees only
    GROUP BY t.PLT_CN, t.SPCD
    ORDER BY t.PLT_CN, TOTAL_CARBON DESC;
    """
    
    conn = get_database_connection()
    df = pd.read_sql_query(query, conn)
    conn.close()
    
    return df

def export_to_qgis(gdf, filename, output_dir="outputs"):
    """
    Export GeoDataFrame to QGIS-compatible format
    
    Args:
        gdf: GeoDataFrame to export
        filename: Output filename (without extension)
        output_dir: Output directory
        
    Returns:
        Path to exported file
    """
    output_path = Path(output_dir) / f"{filename}.geojson"
    output_path.parent.mkdir(exist_ok=True)
    
    gdf.to_file(output_path, driver='GeoJSON')
    print(f"Exported {len(gdf):,} records to: {output_path}")
    
    return output_path

def get_forest_type_summary():
    """Get summary of all forest types in the database"""
    query = """
    SELECT 
        c.FORTYPCD, ft.MEANING as FOREST_TYPE,
        COUNT(*) as PLOT_COUNT,
        AVG(p.ELEV) as AVG_ELEVATION,
        MIN(p.ELEV) as MIN_ELEVATION,
        MAX(p.ELEV) as MAX_ELEVATION
    FROM COND c
    JOIN PLOT p ON c.PLT_CN = p.CN
    JOIN REF_FOREST_TYPE ft ON c.FORTYPCD = ft.VALUE
    WHERE p.PLOT_STATUS_CD = 1
    GROUP BY c.FORTYPCD, ft.MEANING
    ORDER BY PLOT_COUNT DESC;
    """
    
    conn = get_database_connection()
    df = pd.read_sql_query(query, conn)
    conn.close()
    
    return df

def get_ownership_summary():
    """Get summary of plot ownership distribution"""
    query = """
    SELECT 
        c.OWNCD as OWNERSHIP_CD,
        COUNT(*) as PLOT_COUNT,
        AVG(p.ELEV) as AVG_ELEVATION,
        COUNT(DISTINCT c.FORTYPCD) as FOREST_TYPES_COUNT
    FROM COND c
    JOIN PLOT p ON c.PLT_CN = p.CN
    WHERE p.PLOT_STATUS_CD = 1
    GROUP BY c.OWNCD
    ORDER BY PLOT_COUNT DESC;
    """
    
    conn = get_database_connection()
    df = pd.read_sql_query(query, conn)
    conn.close()
    
    # Add ownership type labels
    df['OWNERSHIP_TYPE'] = df['OWNERSHIP_CD'].map({
        11: 'National Forest Service',
        21: 'Other Federal',
        31: 'State and Local Government',
        32: 'Private Corporation', 
        46: 'Private Individual'
    }).fillna('Unknown')
    
    return df

# Example usage functions
def example_douglas_fir_plots():
    """Example: Get all Douglas-fir plots above 2000 feet"""
    return query_plots_by_criteria(
        forest_types=[201], 
        elevation_range=(2000, 8000)
    )

def example_national_forest_plots():
    """Example: Get National Forest plots from recent inventory"""
    return query_plots_by_criteria(
        ownership_codes=[11],
        inventory_years=[2017, 2018, 2019]
    )

def example_steep_slope_plots():
    """Example: Get plots on steep slopes (>30%) for erosion analysis"""
    return query_plots_by_criteria(
        slope_range=(30, 100)
    )

def example_north_facing_plots():
    """Example: Get north-facing plots for cool microclimate analysis"""
    return query_plots_by_criteria(
        aspect_range=(315, 45)  # North-facing aspects
    )

if __name__ == "__main__":
    # Test the interface
    print("Testing FIA Database Interface...")
    print("="*50)
    
    # Test forest type summary
    print("\nForest Type Summary (top 10):")
    forest_summary = get_forest_type_summary()
    print(forest_summary.head(10).to_string(index=False))
    
    # Test ownership summary
    print(f"\nOwnership Summary:")
    ownership_summary = get_ownership_summary()
    print(ownership_summary.to_string(index=False))
    
    # Test Douglas-fir query
    print(f"\nTesting plot queries...")
    df_plots = example_douglas_fir_plots()
    print(f"Douglas-fir plots above 2000 feet: {len(df_plots):,}")
    
    # Test National Forest query
    nf_plots = example_national_forest_plots()
    print(f"National Forest plots from recent years: {len(nf_plots):,}")
    
    # Test steep slopes
    steep_plots = example_steep_slope_plots()
    print(f"Plots on steep slopes (>30%): {len(steep_plots):,}")
    
    print(f"\n✅ Database interface working correctly!")
    print(f"\nUsage examples:")
    print(f"  from scripts.database.fia_database_interface import *")
    print(f"  plots = query_plots_by_criteria(forest_types=[201], elevation_range=(2000, 4000))")
    print(f"  export_to_qgis(plots, 'my_analysis')") 