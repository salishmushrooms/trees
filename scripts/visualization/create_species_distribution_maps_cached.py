#!/usr/bin/env python3
"""
CACHED Species Distribution Mapping - Optimized for Multiple Species

This optimized version caches expensive masks that are shared across species:
- Hydrography mask (water bodies) - universal for all species
- Elevation masks - cached per threshold as needed
- Tree canopy cover mask - cached for 20% threshold

Masks are automatically cached on first use and reused in subsequent runs.

Usage:
  python create_species_distribution_maps_cached.py  # Uses caching automatically
  python create_species_distribution_maps_cached.py --no-cache  # Skip caching
"""

import sqlite3
import pandas as pd
import geopandas as gpd
import numpy as np
from pathlib import Path
from shapely.geometry import Point
from shapely.ops import unary_union
from tqdm import tqdm
import warnings
import argparse
import json
from datetime import datetime
from shapely.geometry import mapping
warnings.filterwarnings('ignore')

# GeoJSON optimization functions for MBTiles
def round_coordinates(geom, precision=6):
    """Round coordinates to specified decimal places for smaller file sizes"""
    def round_coords(coords):
        return [round(float(coord), precision) for coord in coords]
    
    def round_nested_coords(coords_list):
        if not coords_list:
            return coords_list
        
        if isinstance(coords_list[0], (list, tuple)):
            return [round_nested_coords(sublist) for sublist in coords_list]
        else:
            return round_coords(coords_list)
    
    geom_dict = mapping(geom)
    coords = geom_dict['coordinates']
    
    if geom_dict['type'] == 'Point':
        geom_dict['coordinates'] = round_coords(coords)
    elif geom_dict['type'] in ['LineString', 'MultiPoint']:
        geom_dict['coordinates'] = round_nested_coords(coords)
    elif geom_dict['type'] in ['Polygon', 'MultiLineString']:
        geom_dict['coordinates'] = [round_nested_coords(ring) for ring in coords]
    elif geom_dict['type'] == 'MultiPolygon':
        geom_dict['coordinates'] = [[round_nested_coords(ring) for ring in polygon] for polygon in coords]
    
    return geom_dict

def get_essential_properties(properties):
    """Keep only essential properties for web mapping"""
    essential_props = {
        'species': properties.get('species'),
        'species_code': properties.get('species_code'),
        'total_plots': properties.get('total_plots')
    }
    # Remove None values
    return {k: v for k, v in essential_props.items() if v is not None}

def save_optimized_geojson(gdf, output_file, precision=6):
    """Save GeoDataFrame as optimized GeoJSON for MBTiles"""
    print(f"  💾 Saving optimized GeoJSON (precision: {precision} decimal places)...")
    
    # Create optimized features
    optimized_features = []
    for idx, row in gdf.iterrows():
        # Round coordinates
        optimized_geom = round_coordinates(row.geometry, precision)
        
        # Filter properties
        properties = get_essential_properties(dict(row.drop('geometry')))
        
        feature = {
            "type": "Feature",
            "geometry": optimized_geom,
            "properties": properties
        }
        optimized_features.append(feature)
    
    # Create optimized GeoJSON structure
    optimized_geojson = {
        "type": "FeatureCollection",
        "features": optimized_features
    }
    
    # Save with compact JSON format
    with open(output_file, 'w') as f:
        json.dump(optimized_geojson, f, separators=(',', ':'))
    
    return output_file

# Area loss diagnostic functions
def calculate_area_km2(geometry):
    """Calculate area in square kilometers from WGS84 geometry"""
    # Rough conversion: 1 sq degree ≈ 12,321 km² at 45°N latitude (Pacific Northwest)
    if hasattr(geometry, 'geoms'):
        area_degrees = sum(poly.area for poly in geometry.geoms)
    else:
        area_degrees = geometry.area
    
    # More accurate conversion for Pacific Northwest (around 45°N)
    area_km2 = area_degrees * 12321  # km² per sq degree at ~45°N
    return area_km2

def log_area_change(step_name, original_geom, modified_geom, original_area_km2=None):
    """Log area changes with detailed diagnostics"""
    if original_area_km2 is None:
        original_area_km2 = calculate_area_km2(original_geom)
    
    modified_area_km2 = calculate_area_km2(modified_geom)
    
    area_lost_km2 = original_area_km2 - modified_area_km2
    area_reduction_pct = (area_lost_km2 / original_area_km2 * 100) if original_area_km2 > 0 else 0
    
    # Count polygons
    original_count = len(original_geom.geoms) if hasattr(original_geom, 'geoms') else 1
    modified_count = len(modified_geom.geoms) if hasattr(modified_geom, 'geoms') else 1
    
    print(f"    📊 {step_name} IMPACT:")
    print(f"      🏞️  Area: {original_area_km2:,.0f} km² → {modified_area_km2:,.0f} km²")
    print(f"      📉 Lost: {area_lost_km2:,.0f} km² ({area_reduction_pct:.1f}% reduction)")
    print(f"      🔢 Polygons: {original_count:,} → {modified_count:,}")
    
    # Warning thresholds
    if area_reduction_pct > 70:
        print(f"      ⚠️  WARNING: Very high area loss ({area_reduction_pct:.1f}%)")
        print(f"      💡 Consider adjusting {step_name.lower()} parameters")
    elif area_reduction_pct > 50:
        print(f"      ⚠️  NOTICE: High area loss ({area_reduction_pct:.1f}%)")
        
    return modified_area_km2

def log_cumulative_loss(original_area_km2, current_area_km2, steps_completed, pause_on_high_loss=False):
    """Log cumulative area loss from original"""
    total_lost_km2 = original_area_km2 - current_area_km2
    total_reduction_pct = (total_lost_km2 / original_area_km2 * 100) if original_area_km2 > 0 else 0
    
    print(f"    📈 CUMULATIVE LOSS after {len(steps_completed)} steps:")
    print(f"      🏞️  Original: {original_area_km2:,.0f} km²")
    print(f"      🏞️  Current: {current_area_km2:,.0f} km²") 
    print(f"      📉 Total lost: {total_lost_km2:,.0f} km² ({total_reduction_pct:.1f}%)")
    print(f"      📋 Steps: {' → '.join(steps_completed)}")
    
    if total_reduction_pct > 80:
        print(f"      🚨 CRITICAL: Excessive habitat loss! Consider relaxing constraints.")
        if pause_on_high_loss:
            print(f"      ⏸️  Pausing for review - press Enter to continue or Ctrl+C to stop...")
            input()
    elif total_reduction_pct > 60:
        print(f"      ⚠️  WARNING: High habitat loss - review constraints carefully.")
        if pause_on_high_loss and total_reduction_pct > 70:
            print(f"      ⏸️  Pausing for review - press Enter to continue or Ctrl+C to stop...")
            input()
        
    return total_reduction_pct

# Try to import rasterio for elevation processing
try:
    import rasterio
    import rasterio.features
    from rasterio.mask import mask
    HAS_RASTERIO = True
except ImportError:
    HAS_RASTERIO = False
    print("⚠️  Rasterio not available - elevation masking will be skipped")

# Try to import scipy for morphological operations
try:
    from scipy import ndimage
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    print("⚠️  Scipy not available - tree cover processing will be less optimized")

# Configuration
# Multiple state databases for Pacific Northwest coverage
DB_PATHS = {
    'WA': "data/raw/trees_SQLite_FIADB_WA.db",
    'OR': "data/raw/trees_SQLite_FIADB_OR.db", 
    'ID': "data/raw/trees_SQLite_FIADB_ID.db"
}
OUTPUT_DIR = Path("outputs/species_distribution")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Cache directory for reusable masks
CACHE_DIR = Path("cache/species_masks")
CACHE_DIR.mkdir(parents=True, exist_ok=True)

# Species configuration
DOUGLAS_FIR_CODE = 202
DOUGLAS_FIR_NAME = "Douglas Fir"

PONDEROSA_PINE_CODE = 122
PONDEROSA_PINE_NAME = "Ponderosa Pine"

BLACK_COTTONWOOD_CODE = 747
BLACK_COTTONWOOD_NAME = "Black Cottonwood"

# Buffer configuration (in degrees, approximately)
# Multiple buffer size options for different use cases
BUFFER_OPTIONS = {
    'small': {'km': 1.0, 'desc': 'Local habitat (1km radius)'},
    'medium': {'km': 2.0, 'desc': 'Regional habitat (2km radius)'},  
    'large': {'km': 4.0, 'desc': 'Landscape connectivity (4km radius)'},
    'xlarge': {'km': 6.0, 'desc': 'Broad habitat corridors (6km radius)'},
    'custom': {'km': 8.0, 'desc': 'Custom size (8km radius)'}
}

# Default buffer size - can be overridden
DEFAULT_BUFFER_SIZE = 'medium'  # Changed from 8km to 2km default
BUFFER_RADIUS_KM = BUFFER_OPTIONS[DEFAULT_BUFFER_SIZE]['km']
BUFFER_RADIUS_DEG = BUFFER_RADIUS_KM / 111.0  # Convert km to degrees (approx)

# Workflow configuration
WORKFLOW_MODE = 'full'  # Options: 'buffers_only', 'staged', 'full'
EXPORT_INTERMEDIATE_FILES = True  # Export files for QGIS editing
SKIP_CONSTRAINTS_IF_EDITING = True  # Skip constraints when exporting for editing

# Elevation configuration
DEM_FILE = "/Users/JJC/morel-maps/QGIS/static-layers/DEM-western-us-combined.vrt"
DOUGLAS_FIR_MAX_ELEVATION_FT = 6000  # Douglas Fir upper elevation limit
PONDEROSA_PINE_MAX_ELEVATION_FT = 4000  # Ponderosa Pine upper elevation limit
MAX_ELEVATION_FT = PONDEROSA_PINE_MAX_ELEVATION_FT  # Current species elevation limit
APPLY_ELEVATION_MASK = True  # Set to False to skip elevation masking

# Tree canopy cover configuration
# Option 1: Use pre-created regional masks (RECOMMENDED for multiple species)
USE_REGIONAL_TCC_MASK = True  # Use pre-created regional masks instead of creating species-specific masks
REGIONAL_TCC_MASK_FILE = None  # Auto-detect based on parameters, or specify path

# Option 2: Create species-specific masks (legacy mode)
TCC_FILE = "data/raw/nlcd_tcc_conus_2021_v2021-4.tif"
MIN_TREE_COVER_PCT = 5  # Minimum tree canopy cover percentage (updated to 5% for current processing)
MIN_LOW_COVER_AREA_SQKM = 2.0  # Minimum area (sq km) for low-cover polygons to include 
TCC_PROCESSING_RESOLUTION = 30  # Process at native NLCD resolution (30m) for maximum detail
APPLY_TCC_MASK = True  # Set to False to skip tree canopy cover masking

# Legacy high-resolution TCC processing settings (for species-specific masks)
TCC_HIGH_RES_MODE = True  # Enable high-resolution processing for fine-grained detail
TCC_MORPHOLOGICAL_SMOOTHING = True  # Apply morphological operations to reduce noise while preserving detail

# Hydrography configuration - Pacific Northwest multi-state coverage
# Multiple potential sources - script will try them in order
HYDROGRAPHY_SOURCES = [
    # Option 1: NHD Washington State Water Bodies
    {
        'file': 'data/raw/hydrography/NHD_H_Washington_State_Shape/Shape/NHDWaterbody.shp',
        'type': 'vector',
        'description': 'NHD Washington State Water Bodies'
    },
    # Option 2: NHD Oregon State Water Bodies
    {
        'file': 'data/raw/hydrography/NHD_H_Oregon_State_Shape/Shape/NHDWaterbody.shp',
        'type': 'vector',
        'description': 'NHD Oregon State Water Bodies'
    },
    # Option 3: NHD Idaho State Water Bodies
    {
        'file': 'data/raw/hydrography/NHD_H_Idaho_State_Shape/Shape/NHDWaterbody.shp',
        'type': 'vector',
        'description': 'NHD Idaho State Water Bodies'
    },


    # Option 6: Simple built-in water polygons (DISABLED - too aggressive)
    # {
    #     'file': 'built-in',
    #     'type': 'built-in',
    #     'description': 'Simple built-in water polygons'
    # }
]
APPLY_HYDROGRAPHY_MASK = False  # Set to False to skip hydrography masking (tree canopy handles water exclusion)

# Geographic bounds for Pacific Northwest (WA, OR, ID)
# Includes US-Canada border clipping at 49°N and Oregon southern border at 42°N
PACIFIC_NW_BOUNDS = [-125.0, 42.0, -110.0, 49.0]  # W, S, E, N (WGS84)
US_CANADA_BORDER = 49.0  # Latitude to clip at northern border

# Simplification tolerance (for reducing polygon complexity)
SIMPLIFICATION_TOLERANCE = 0.002  # degrees (increased to reduce small artifacts)

# Cache management functions
def get_cache_filename(mask_type, **params):
    """Generate cache filename for a mask type and parameters"""
    if mask_type == "hydrography":
        return CACHE_DIR / "pnw_hydrography_mask.geojson"
    elif mask_type == "elevation":
        elevation_ft = params['elevation_ft']
        # Handle both old format (single number) and new format (min_max)
        if '_' in str(elevation_ft):
            return CACHE_DIR / f"pnw_elevation_mask_{elevation_ft}ft.geojson"
        else:
            return CACHE_DIR / f"pnw_elevation_mask_{elevation_ft}ft.geojson"
    elif mask_type == "tcc":
        tcc_pct = params['tcc_pct']
        resolution = params.get('resolution', TCC_PROCESSING_RESOLUTION)
        min_area = params.get('min_area_sqkm', MIN_LOW_COVER_AREA_SQKM)
        return CACHE_DIR / f"pnw_tcc_mask_{tcc_pct}pct_{resolution}m_{min_area}sqkm.geojson"
    else:
        raise ValueError(f"Unknown mask type: {mask_type}")

def load_cached_mask(mask_type, **params):
    """Load a cached mask if it exists"""
    cache_file = get_cache_filename(mask_type, **params)
    
    if cache_file.exists():
        print(f"    📦 Loading cached {mask_type} mask: {cache_file.name}")
        file_size_mb = cache_file.stat().st_size / (1024 * 1024)
        
        # Different size limits for different mask types
        size_limit_mb = {
            'elevation': 250,     # Elevation masks are large due to high-res DEM data
            'hydrography': 800,   # Hydrography masks can be very large
            'tcc': 300           # Tree cover masks can be substantial
        }.get(mask_type, 300)
        
        # Skip very large cache files that are likely to fail
        if file_size_mb > size_limit_mb:
            print(f"    ⚠️  Cache file too large ({file_size_mb:.1f} MB > {size_limit_mb} MB) - skipping for stability")
            print(f"    🗑️  Consider clearing cache: rm {cache_file}")
            return None
        
        try:
            # Try to load with error recovery and validation
            print(f"    🔄 Loading {file_size_mb:.1f} MB cache file...")
            gdf = gpd.read_file(cache_file)
            
            if len(gdf) == 0:
                print(f"    ⚠️  Cache file is empty")
                cache_file.unlink()
                print(f"    🗑️  Deleted empty cache file")
                return None
                
            geometry = gdf.iloc[0].geometry
            if geometry is None or geometry.is_empty:
                print(f"    ⚠️  Cache file has invalid geometry")
                cache_file.unlink()
                print(f"    🗑️  Deleted invalid cache file")
                return None
            
            # Validate geometry
            if not geometry.is_valid:
                print(f"    🔧 Cache geometry is invalid - attempting to fix...")
                from shapely.validation import make_valid
                geometry = make_valid(geometry)
                if not geometry.is_valid:
                    print(f"    ❌ Could not fix invalid geometry")
                    cache_file.unlink()
                    print(f"    🗑️  Deleted unfixable cache file")
                    return None
                print(f"    ✅ Fixed invalid geometry")
            
            print(f"    ✅ Loaded cached mask ({file_size_mb:.1f} MB)")
            return geometry
            
        except Exception as e:
            print(f"    ❌ Error loading cache ({file_size_mb:.1f} MB): {e}")
            print(f"    🔄 Will rebuild mask from scratch...")
            # Delete corrupted cache file
            try:
                cache_file.unlink()
                print(f"    🗑️  Deleted corrupted cache file")
            except:
                pass
            return None
    
    return None

def save_cached_mask(mask_geometry, mask_type, **params):
    """Save a mask to cache"""
    if mask_geometry is None:
        return
    
    cache_file = get_cache_filename(mask_type, **params)
    
    print(f"    💾 Saving {mask_type} mask to cache: {cache_file.name}")
    
    # Create metadata
    metadata = {
        'mask_type': mask_type,
        'created_date': datetime.now().isoformat(),
        'bounds': PACIFIC_NW_BOUNDS,
        **params
    }
    
    # Create GeoDataFrame
    gdf = gpd.GeoDataFrame([metadata], geometry=[mask_geometry], crs='EPSG:4326')
    
    try:
        # For large geometries, simplify slightly before caching to reduce file size
        if mask_type == 'elevation' and hasattr(mask_geometry, 'area'):
            # Very light simplification for elevation masks (they tend to be large)
            simplified_geom = mask_geometry.simplify(0.0005, preserve_topology=True)
            if simplified_geom.is_valid and not simplified_geom.is_empty:
                print(f"    🎯 Applied light simplification to reduce cache size")
                mask_geometry = simplified_geom
        
        gdf.loc[0, 'geometry'] = mask_geometry  # Update with potentially simplified geometry
        gdf.to_file(cache_file, driver='GeoJSON')
        file_size_mb = cache_file.stat().st_size / (1024 * 1024)
        print(f"    ✅ Cached {mask_type} mask ({file_size_mb:.1f} MB)")
        
        # Warn if cache file is getting large
        size_limit_mb = {
            'elevation': 250,
            'hydrography': 800,
            'tcc': 300
        }.get(mask_type, 300)
        
        if file_size_mb > size_limit_mb * 0.8:  # 80% of limit
            print(f"    ⚠️  Cache file is large ({file_size_mb:.1f} MB)")
            print(f"    💡 This is normal for {mask_type} masks covering the Pacific Northwest")
            
    except Exception as e:
        print(f"    ❌ Error saving cache: {e}")
        # Delete partial cache file if it exists
        if cache_file.exists():
            try:
                cache_file.unlink()
                print(f"    🗑️  Cleaned up partial cache file")
            except:
                pass

def calculate_species_elevation_limits(species_code, species_name):
    """Calculate species-specific elevation limits from FIA data"""
    print(f"📊 Calculating elevation limits for {species_name} from FIA data...")
    
    query = """
    WITH latest_plots AS (
        SELECT 
            STATECD, UNITCD, COUNTYCD, PLOT,
            MAX(INVYR) as LATEST_INVYR
        FROM PLOT 
        WHERE PLOT_STATUS_CD = 1
            AND LAT IS NOT NULL 
            AND LON IS NOT NULL
            AND ELEV IS NOT NULL
        GROUP BY STATECD, UNITCD, COUNTYCD, PLOT
    )
    SELECT 
        p.ELEV as ELEVATION_FT
    FROM PLOT p
    INNER JOIN latest_plots lp ON (
        p.STATECD = lp.STATECD AND 
        p.UNITCD = lp.UNITCD AND 
        p.COUNTYCD = lp.COUNTYCD AND 
        p.PLOT = lp.PLOT AND 
        p.INVYR = lp.LATEST_INVYR
    )
    INNER JOIN TREE t ON p.CN = t.PLT_CN
    WHERE t.SPCD = ?
        AND t.STATUSCD = 1  -- Live trees only
        AND p.ELEV IS NOT NULL
        AND p.ELEV > 0  -- Valid elevations only
    ORDER BY p.ELEV;
    """
    
    # Collect elevation data from all states
    all_elevations = []
    
    for state_code, db_path in DB_PATHS.items():
        if not Path(db_path).exists():
            print(f"  ⚠️  {state_code} database not found: {db_path}")
            continue
            
        try:
            conn = sqlite3.connect(db_path)
            df = pd.read_sql_query(query, conn, params=[species_code])
            conn.close()
            
            elevations = df['ELEVATION_FT'].dropna().tolist()
            all_elevations.extend(elevations)
            
            print(f"    {state_code}: {len(elevations):,} elevation records")
            
        except Exception as e:
            print(f"  ❌ Error loading {state_code} elevations: {e}")
            continue
    
    if not all_elevations:
        print(f"  ❌ No elevation data found - using default limits")
        return None, None
    
    # Calculate actual min/max from plot data
    import numpy as np
    all_elevations = np.array(all_elevations)
    
    percentile_5 = np.percentile(all_elevations, 5)
    percentile_95 = np.percentile(all_elevations, 95)
    min_elev = np.min(all_elevations)
    max_elev = np.max(all_elevations)
    
    print(f"  📈 Elevation Analysis:")
    print(f"    Range: {min_elev:.0f} - {max_elev:.0f} ft")
    print(f"    5th percentile: {percentile_5:.0f} ft")
    print(f"    95th percentile: {percentile_95:.0f} ft")
    
    # Use actual min/max from plot data, rounded to 500ft intervals
    min_limit = max(0, int((min_elev // 500) * 500))  # Round down, min 0
    max_limit = int(((max_elev // 500) + 1) * 500)   # Round up
    
    print(f"  🎯 Using actual data limits: {min_limit:,} - {max_limit:,} ft")
    print(f"    This captures 100% of plot observations (full species range)")
    
    return min_limit, max_limit

def get_species_plots(species_code, species_name):
    """Get all species plot locations from multiple Pacific Northwest state databases"""
    print(f"🌲 Extracting {species_name} plots from Pacific Northwest FIA databases...")
    
    query = """
    WITH latest_plots AS (
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
        p.LAT,
        p.LON,
        p.ELEV as ELEVATION_FT,
        COUNT(t.CN) as TREE_COUNT,
        SUM(t.CARBON_AG) as TOTAL_CARBON_AG,
        SUM(t.DIA * t.DIA * 0.005454) as BASAL_AREA_SQFT
    FROM PLOT p
    INNER JOIN latest_plots lp ON (
        p.STATECD = lp.STATECD AND 
        p.UNITCD = lp.UNITCD AND 
        p.COUNTYCD = lp.COUNTYCD AND 
        p.PLOT = lp.PLOT AND 
        p.INVYR = lp.LATEST_INVYR
    )
    INNER JOIN TREE t ON p.CN = t.PLT_CN
    WHERE t.SPCD = ?
        AND t.STATUSCD = 1  -- Live trees only
        AND p.LAT IS NOT NULL 
        AND p.LON IS NOT NULL
    GROUP BY p.LAT, p.LON, p.ELEV
    ORDER BY TOTAL_CARBON_AG DESC;
    """
    
    # Combine data from all available state databases
    all_plots = []
    state_counts = {}
    
    for state_code, db_path in DB_PATHS.items():
        if not Path(db_path).exists():
            print(f"  ⚠️  {state_code} database not found: {db_path}")
            continue
            
        try:
            print(f"  📁 Loading {state_code} data from: {db_path}")
            conn = sqlite3.connect(db_path)
            df = pd.read_sql_query(query, conn, params=[species_code])
            conn.close()
            
            # Add state information
            df['STATE'] = state_code
            all_plots.append(df)
            state_counts[state_code] = len(df)
            
            print(f"    📍 {state_code}: {len(df):,} plots")
            
        except Exception as e:
            print(f"  ❌ Error loading {state_code} database: {e}")
            continue
    
    if not all_plots:
        print(f"  ❌ No databases found or accessible")
        return pd.DataFrame()
    
    # Combine all state data
    combined_df = pd.concat(all_plots, ignore_index=True)
    
    # Apply Pacific Northwest geographic bounds clipping
    print(f"  🗺️  Applying Pacific Northwest bounds: {PACIFIC_NW_BOUNDS}")
    bounds_mask = (
        (combined_df['LON'] >= PACIFIC_NW_BOUNDS[0]) &  # West
        (combined_df['LAT'] >= PACIFIC_NW_BOUNDS[1]) &  # South
        (combined_df['LON'] <= PACIFIC_NW_BOUNDS[2]) &  # East
        (combined_df['LAT'] <= PACIFIC_NW_BOUNDS[3])    # North (US-Canada border)
    )
    
    combined_df = combined_df[bounds_mask]
    
    print(f"  📊 Pacific Northwest Summary:")
    for state_code in state_counts.keys():
        state_plots = len(combined_df[combined_df['STATE'] == state_code])
        print(f"    {state_code}: {state_plots:,} plots (after bounds clipping)")
    
    print(f"  📍 Total: {len(combined_df):,} plots with {species_name}")
    print(f"  🌿 Total trees: {combined_df['TREE_COUNT'].sum():,}")
    print(f"  💚 Total carbon: {combined_df['TOTAL_CARBON_AG'].sum():.1f}")
    
    return combined_df

def get_all_plots_with_species_info(target_species_codes=None):
    """Get ALL FIA plots with information about species presence for comparison analysis"""
    print(f"📊 Extracting ALL Pacific Northwest FIA plots for survey coverage analysis...")
    
    if target_species_codes is None:
        target_species_codes = [DOUGLAS_FIR_CODE, PONDEROSA_PINE_CODE]
    
    # Create species code placeholders for query
    species_placeholders = ",".join("?" * len(target_species_codes))
    
    query = f"""
    WITH latest_plots AS (
        SELECT 
            STATECD, UNITCD, COUNTYCD, PLOT,
            MAX(INVYR) as LATEST_INVYR
        FROM PLOT 
        WHERE PLOT_STATUS_CD = 1
            AND LAT IS NOT NULL 
            AND LON IS NOT NULL
        GROUP BY STATECD, UNITCD, COUNTYCD, PLOT
    ),
    plot_species_summary AS (
        SELECT 
            p.CN as PLOT_CN,
            p.LAT,
            p.LON,
            p.ELEV as ELEVATION_FT,
            p.INVYR,
            COUNT(t.CN) as TOTAL_TREES,
            COUNT(DISTINCT t.SPCD) as SPECIES_COUNT,
            SUM(CASE WHEN t.SPCD IN ({species_placeholders}) THEN 1 ELSE 0 END) as TARGET_SPECIES_TREES,
            GROUP_CONCAT(DISTINCT CASE WHEN t.SPCD IN ({species_placeholders}) THEN t.SPCD END) as TARGET_SPECIES_PRESENT,
            SUM(t.CARBON_AG) as TOTAL_CARBON_AG,
            SUM(t.DIA * t.DIA * 0.005454) as TOTAL_BASAL_AREA_SQFT
        FROM PLOT p
        INNER JOIN latest_plots lp ON (
            p.STATECD = lp.STATECD AND 
            p.UNITCD = lp.UNITCD AND 
            p.COUNTYCD = lp.COUNTYCD AND 
            p.PLOT = lp.PLOT AND 
            p.INVYR = lp.LATEST_INVYR
        )
        LEFT JOIN TREE t ON p.CN = t.PLT_CN AND t.STATUSCD = 1  -- Live trees only
        WHERE p.LAT IS NOT NULL 
            AND p.LON IS NOT NULL
        GROUP BY p.CN, p.LAT, p.LON, p.ELEV, p.INVYR
    )
    SELECT 
        LAT,
        LON,
        ELEVATION_FT,
        INVYR as SURVEY_YEAR,
        TOTAL_TREES,
        SPECIES_COUNT,
        TARGET_SPECIES_TREES,
        TARGET_SPECIES_PRESENT,
        TOTAL_CARBON_AG,
        TOTAL_BASAL_AREA_SQFT,
        CASE 
            WHEN TARGET_SPECIES_TREES > 0 THEN 1 
            ELSE 0 
        END as HAS_TARGET_SPECIES
    FROM plot_species_summary
    ORDER BY LAT DESC, LON ASC;
    """
    
    # Combine data from all available state databases
    all_plots = []
    state_counts = {}
    
    query_params = target_species_codes + target_species_codes  # Doubled for the two IN clauses
    
    for state_code, db_path in DB_PATHS.items():
        if not Path(db_path).exists():
            print(f"  ⚠️  {state_code} database not found: {db_path}")
            continue
            
        try:
            print(f"  📁 Loading {state_code} data from: {db_path}")
            conn = sqlite3.connect(db_path)
            df = pd.read_sql_query(query, conn, params=query_params)
            conn.close()
            
            # Add state information
            df['STATE'] = state_code
            all_plots.append(df)
            
            # Calculate summary stats
            total_plots = len(df)
            plots_with_target = len(df[df['HAS_TARGET_SPECIES'] == 1])
            plots_without_target = total_plots - plots_with_target
            
            state_counts[state_code] = {
                'total': total_plots,
                'with_target': plots_with_target,
                'without_target': plots_without_target
            }
            
            print(f"    📍 {state_code}: {total_plots:,} total plots")
            print(f"      🌲 {plots_with_target:,} with target species ({plots_with_target/total_plots*100:.1f}%)")
            print(f"      🚫 {plots_without_target:,} without target species ({plots_without_target/total_plots*100:.1f}%)")
            
        except Exception as e:
            print(f"  ❌ Error loading {state_code} database: {e}")
            continue
    
    if not all_plots:
        print(f"  ❌ No databases found or accessible")
        return pd.DataFrame()
    
    # Combine all state data
    combined_df = pd.concat(all_plots, ignore_index=True)
    
    # Apply Pacific Northwest geographic bounds clipping
    print(f"  🗺️  Applying Pacific Northwest bounds: {PACIFIC_NW_BOUNDS}")
    bounds_mask = (
        (combined_df['LON'] >= PACIFIC_NW_BOUNDS[0]) &  # West
        (combined_df['LAT'] >= PACIFIC_NW_BOUNDS[1]) &  # South
        (combined_df['LON'] <= PACIFIC_NW_BOUNDS[2]) &  # East
        (combined_df['LAT'] <= PACIFIC_NW_BOUNDS[3])    # North (US-Canada border)
    )
    
    combined_df = combined_df[bounds_mask]
    
    # Add species name information
    species_name_map = {
        DOUGLAS_FIR_CODE: 'Douglas Fir',
        PONDEROSA_PINE_CODE: 'Ponderosa Pine'
    }
    
    def decode_species_present(species_str):
        if pd.isna(species_str) or species_str == '':
            return 'None'
        species_codes = [int(x) for x in str(species_str).split(',') if x.strip().isdigit()]
        species_names = [species_name_map.get(code, f'Species_{code}') for code in species_codes]
        return ', '.join(species_names)
    
    combined_df['TARGET_SPECIES_NAMES'] = combined_df['TARGET_SPECIES_PRESENT'].apply(decode_species_present)
    
    print(f"\n  📊 Pacific Northwest TOTAL Summary:")
    total_plots = len(combined_df)
    total_with_target = len(combined_df[combined_df['HAS_TARGET_SPECIES'] == 1])
    total_without_target = total_plots - total_with_target
    
    for state_code in state_counts.keys():
        state_plots = len(combined_df[combined_df['STATE'] == state_code])
        state_with_target = len(combined_df[(combined_df['STATE'] == state_code) & (combined_df['HAS_TARGET_SPECIES'] == 1)])
        print(f"    {state_code}: {state_plots:,} plots, {state_with_target:,} with target species")
    
    print(f"  📍 TOTAL: {total_plots:,} plots surveyed")
    print(f"  🌲 {total_with_target:,} plots WITH target species ({total_with_target/total_plots*100:.1f}%)")
    print(f"  🚫 {total_without_target:,} plots WITHOUT target species ({total_without_target/total_plots*100:.1f}%)")
    print(f"  🗓️  Survey years: {combined_df['SURVEY_YEAR'].min()}-{combined_df['SURVEY_YEAR'].max()}")
    
    return combined_df

def export_all_plots_comparison(target_species_codes=None, suffix=""):
    """Export all plots for comparison with species-specific habitat maps"""
    print(f"📤 EXPORTING ALL PLOTS FOR SURVEY COVERAGE ANALYSIS")
    print("="*60)
    
    # Get all plots data
    all_plots_df = get_all_plots_with_species_info(target_species_codes)
    
    if len(all_plots_df) == 0:
        print(f"❌ No plot data found")
        return None
    
    # Create GeoDataFrame
    plot_points_gdf = gpd.GeoDataFrame(
        all_plots_df, 
        geometry=[Point(lon, lat) for lon, lat in zip(all_plots_df['LON'], all_plots_df['LAT'])],
        crs='EPSG:4326'
    )
    
    # Create output filename
    species_names = []
    if target_species_codes:
        for code in target_species_codes:
            if code == DOUGLAS_FIR_CODE:
                species_names.append("douglas_fir")
            elif code == PONDEROSA_PINE_CODE:
                species_names.append("ponderosa_pine")
            else:
                species_names.append(f"species_{code}")
    else:
        species_names = ["douglas_fir", "ponderosa_pine"]
    
    species_suffix = "_".join(species_names)
    output_file = OUTPUT_DIR / f"ALL_plots_pnw_{species_suffix}_comparison{suffix}.geojson"
    
    # Export the file
    plot_points_gdf.to_file(output_file)
    
    file_size_mb = output_file.stat().st_size / (1024 * 1024)
    
    print(f"\n✅ ALL PLOTS EXPORTED!")
    print(f"📁 Output file: {output_file}")
    print(f"📊 Contains: {len(plot_points_gdf):,} FIA plots")
    print(f"💾 File size: {file_size_mb:.1f} MB") 
    
    # Summary by presence/absence
    with_target = len(plot_points_gdf[plot_points_gdf['HAS_TARGET_SPECIES'] == 1])
    without_target = len(plot_points_gdf) - with_target
    
    print(f"\n📋 Plot Categories:")
    print(f"  🌲 Plots WITH target species: {with_target:,} ({with_target/len(plot_points_gdf)*100:.1f}%)")
    print(f"  🚫 Plots WITHOUT target species: {without_target:,} ({without_target/len(plot_points_gdf)*100:.1f}%)")
    
    print(f"\n🗺️  QGIS Usage Tips:")
    print(f"  1. Load this file alongside your species habitat polygons")
    print(f"  2. Style by 'HAS_TARGET_SPECIES' field:")
    print(f"     • Green dots = plots with target species (should align with habitat polygons)")
    print(f"     • Red dots = plots without target species (gaps in habitat or survey bias)")
    print(f"  3. Use 'TARGET_SPECIES_NAMES' field to see which species were found")
    print(f"  4. Check 'SURVEY_YEAR' to identify temporal coverage gaps")
    print(f"  5. Look for patterns:")
    print(f"     • Red dots inside habitat polygons = species not found in suitable habitat")
    print(f"     • Green dots outside habitat polygons = habitat model may be too restrictive")
    print(f"     • Areas with no dots = no FIA survey coverage")
    
    return output_file

def create_plot_buffers(plots_df):
    """Create buffer polygons around each plot location"""
    print(f"  🔵 Creating {BUFFER_RADIUS_KM}km buffers around {len(plots_df):,} plots...")
    
    # Create GeoDataFrame from plot coordinates
    geometry = [Point(lon, lat) for lon, lat in zip(plots_df['LON'], plots_df['LAT'])]
    gdf = gpd.GeoDataFrame(plots_df, geometry=geometry, crs='EPSG:4326')
    
    # Create buffers (in degrees)
    print(f"  📐 Buffer radius: {BUFFER_RADIUS_DEG:.4f} degrees (~{BUFFER_RADIUS_KM}km)")
    
    with tqdm(total=len(gdf), desc="  Creating buffers", unit="plots") as pbar:
        buffers = []
        for idx, row in gdf.iterrows():
            buffer_poly = row.geometry.buffer(BUFFER_RADIUS_DEG)
            buffers.append(buffer_poly)
            pbar.update(1)
    
    # Create new GeoDataFrame with buffers
    buffer_gdf = gpd.GeoDataFrame(
        plots_df[['TREE_COUNT', 'TOTAL_CARBON_AG', 'BASAL_AREA_SQFT']], 
        geometry=buffers, 
        crs='EPSG:4326'
    )
    
    print(f"  ✅ Created {len(buffer_gdf):,} buffer polygons")
    return buffer_gdf

def merge_overlapping_buffers(buffer_gdf, species_code, species_name):
    """Merge overlapping buffers into continuous areas"""
    print(f"  🔗 Merging overlapping buffers into continuous areas...")
    
    # For large datasets, use chunked processing to avoid memory issues
    num_buffers = len(buffer_gdf)
    print(f"  ⏳ Dissolving {num_buffers:,} buffer polygons...")
    
    if num_buffers > 5000:
        print(f"    🚀 Using chunked processing for large dataset...")
        # Process in chunks to avoid memory issues
        chunk_size = 2000
        chunks = []
        
        for i in range(0, num_buffers, chunk_size):
            chunk_end = min(i + chunk_size, num_buffers)
            chunk_gdf = buffer_gdf.iloc[i:chunk_end]
            print(f"    Processing chunk {i//chunk_size + 1}/{(num_buffers-1)//chunk_size + 1} ({len(chunk_gdf)} polygons)...")
            
            chunk_union = unary_union(chunk_gdf.geometry.tolist())
            chunks.append(chunk_union)
        
        print(f"    🔗 Merging {len(chunks)} chunks...")
        merged_geometry = unary_union(chunks)
    else:
        # Standard processing for smaller datasets
        merged_geometry = unary_union(buffer_gdf.geometry.tolist())
    
    # Handle both single polygon and multipolygon results
    if hasattr(merged_geometry, 'geoms'):
        # MultiPolygon
        polygons = list(merged_geometry.geoms)
        print(f"  🔢 Result: {len(polygons)} separate polygon areas")
    else:
        # Single Polygon
        polygons = [merged_geometry]
        print(f"  🔢 Result: 1 continuous polygon area")
    
    # Create summary statistics
    total_tree_count = buffer_gdf['TREE_COUNT'].sum()
    total_carbon = buffer_gdf['TOTAL_CARBON_AG'].sum()
    total_basal_area = buffer_gdf['BASAL_AREA_SQFT'].sum()
    
    # Create output GeoDataFrame
    merged_gdf = gpd.GeoDataFrame([{
        'species': species_name,
        'species_code': species_code,
        'total_plots': len(buffer_gdf),
        'total_trees': int(total_tree_count),
        'total_carbon_ag': round(total_carbon, 1),
        'total_basal_area_sqft': round(total_basal_area, 1),
        'buffer_radius_km': BUFFER_RADIUS_KM,
        'polygon_count': len(polygons),
        'geometry': merged_geometry
    }], crs='EPSG:4326')
    
    return merged_gdf

def clip_to_canada_border(gdf):
    """Clip distribution polygons to Canada-US border at 49°N"""
    print(f"  ✂️  Clipping to Canada-US border at {US_CANADA_BORDER}°N...")
    
    from shapely.geometry import box
    
    # Create clipping box at Pacific Northwest bounds
    clip_box = box(*PACIFIC_NW_BOUNDS)
    
    clipped_geometries = []
    for geom in gdf.geometry:
        # Clip geometry to the bounding box
        clipped_geom = geom.intersection(clip_box)
        clipped_geometries.append(clipped_geom)
    
    # Update geometry
    clipped_gdf = gdf.copy()
    clipped_gdf.geometry = clipped_geometries
    
    # Check final bounds
    final_bounds = clipped_gdf.total_bounds
    max_lat = final_bounds[3]
    
    if max_lat > US_CANADA_BORDER + 0.001:  # Small tolerance for floating point
        print(f"    ⚠️  Warning: Still extends {max_lat - US_CANADA_BORDER:.3f}° north of border")
    else:
        print(f"    ✅ Successfully clipped to {max_lat:.3f}°N (≤{US_CANADA_BORDER}°N)")
    
    return clipped_gdf

def simplify_geometry(gdf, tolerance=SIMPLIFICATION_TOLERANCE):
    """Simplify polygon geometry to reduce file size"""
    print(f"  🎨 Simplifying geometry with tolerance {tolerance} degrees...")
    
    from shapely.geometry import Polygon, MultiPolygon, LineString, Point
    
    original_vertices = 0
    simplified_vertices = 0
    
    def count_vertices(geom):
        """Count vertices in a geometry, handling different types"""
        if isinstance(geom, (Polygon,)):
            return len(geom.exterior.coords)
        elif isinstance(geom, (MultiPolygon,)):
            return sum(len(poly.exterior.coords) for poly in geom.geoms if isinstance(poly, Polygon))
        elif isinstance(geom, (LineString,)):
            return len(geom.coords)
        elif isinstance(geom, (Point,)):
            return 1
        elif hasattr(geom, 'geoms'):
            # Mixed geometry collection
            return sum(count_vertices(sub_geom) for sub_geom in geom.geoms)
        else:
            return 0
    
    simplified_geometries = []
    for geom in gdf.geometry:
        # Count original vertices
        original_vertices += count_vertices(geom)
        
        # Filter out non-polygon geometries after clipping
        if isinstance(geom, (Polygon, MultiPolygon)):
            # Standard polygon simplification
            simplified_geom = geom.simplify(tolerance, preserve_topology=True)
        elif hasattr(geom, 'geoms'):
            # Mixed geometry collection - keep only polygons
            polygon_parts = [g for g in geom.geoms if isinstance(g, (Polygon, MultiPolygon))]
            if polygon_parts:
                if len(polygon_parts) == 1:
                    simplified_geom = polygon_parts[0].simplify(tolerance, preserve_topology=True)
                else:
                    from shapely.ops import unary_union
                    combined = unary_union(polygon_parts)
                    simplified_geom = combined.simplify(tolerance, preserve_topology=True)
            else:
                # No polygons found, skip this geometry
                print(f"    ⚠️  Skipping non-polygon geometry after border clipping")
                continue
        else:
            # Skip non-polygon geometries (LineString, Point, etc.)
            print(f"    ⚠️  Skipping {type(geom).__name__} geometry after border clipping")
            continue
        
        simplified_geometries.append(simplified_geom)
        
        # Count simplified vertices
        simplified_vertices += count_vertices(simplified_geom)
    
    # Update geometry
    gdf = gdf.copy()
    gdf = gdf.iloc[:len(simplified_geometries)]  # Adjust for any skipped geometries
    gdf.geometry = simplified_geometries
    
    if original_vertices > 0:
        reduction_pct = 100 * (1 - simplified_vertices / original_vertices)
        print(f"  📉 Reduced vertices: {original_vertices:,} → {simplified_vertices:,} ({reduction_pct:.1f}% reduction)")
    else:
        print(f"  📉 Processed geometries: {len(simplified_geometries)} polygon features")
    
    return gdf

def remove_small_artifacts(gdf, min_area_sqkm=0.1):
    """Remove small polygon artifacts and holes to clean up the output"""
    print(f"  🧹 Removing small artifacts (< {min_area_sqkm} sq km)...")
    
    from shapely.geometry import Polygon, MultiPolygon
    
    cleaned_geometries = []
    for geom in gdf.geometry:
        if isinstance(geom, MultiPolygon):
            # Keep only large polygons from multipolygon
            large_polygons = []
            for poly in geom.geoms:
                if isinstance(poly, Polygon):
                    # Calculate approximate area in sq km
                    area_sqkm = poly.area * 12321  # Rough conversion at 45°N
                    if area_sqkm >= min_area_sqkm:
                        large_polygons.append(poly)
            
            if large_polygons:
                if len(large_polygons) == 1:
                    cleaned_geom = large_polygons[0]
                else:
                    from shapely.ops import unary_union
                    cleaned_geom = unary_union(large_polygons)
            else:
                # No large polygons, skip this geometry
                continue
                
        elif isinstance(geom, Polygon):
            # Check if polygon is large enough
            area_sqkm = geom.area * 12321  # Rough conversion at 45°N
            if area_sqkm >= min_area_sqkm:
                cleaned_geom = geom
            else:
                # Skip small polygons
                continue
        else:
            # Keep non-polygon geometries as-is
            cleaned_geom = geom
        
        # Apply buffer smoothing to remove small protrusions
        try:
            # Small positive buffer to fill gaps, then negative to restore size
            buffer_size = 0.0005  # ~55m buffer
            smoothed_geom = cleaned_geom.buffer(buffer_size, resolution=8)
            smoothed_geom = smoothed_geom.buffer(-buffer_size, resolution=8)
            
            # Ensure result is valid
            if smoothed_geom.is_valid and not smoothed_geom.is_empty:
                cleaned_geometries.append(smoothed_geom)
            else:
                # Fallback to original if smoothing failed
                cleaned_geometries.append(cleaned_geom)
        except:
            # Fallback to original if smoothing failed
            cleaned_geometries.append(cleaned_geom)
    
    if not cleaned_geometries:
        print(f"    ⚠️  No polygons remain after artifact removal")
        return gdf
    
    # Update geometry
    cleaned_gdf = gdf.copy()
    cleaned_gdf = cleaned_gdf.iloc[:len(cleaned_geometries)]
    cleaned_gdf.geometry = cleaned_geometries
    
    original_count = len(gdf)
    cleaned_count = len(cleaned_gdf)
    removed_count = original_count - cleaned_count
    
    print(f"  ✅ Removed {removed_count} small artifacts ({original_count} → {cleaned_count} polygons)")
    
    return cleaned_gdf

def create_elevation_mask(bounds, dem_file=DEM_FILE, min_elevation_ft=None, max_elevation_ft=MAX_ELEVATION_FT, use_cache=True):
    """Create a mask of areas outside the elevation range"""
    if min_elevation_ft is not None:
        print(f"  🏔️  Creating elevation mask for areas outside {min_elevation_ft:,} - {max_elevation_ft:,} feet...")
        cache_key = f"{min_elevation_ft}_{max_elevation_ft}"
    else:
        print(f"  🏔️  Creating elevation mask for areas above {max_elevation_ft:,} feet...")
        cache_key = str(max_elevation_ft)
    
    # Check cache first
    if use_cache:
        cached_mask = load_cached_mask("elevation", elevation_ft=cache_key)
        if cached_mask is not None:
            return cached_mask
    
    if not HAS_RASTERIO:
        print("  ⚠️  Skipping elevation masking - rasterio not available")
        return None
        
    if not Path(dem_file).exists():
        print(f"  ⚠️  Skipping elevation masking - DEM file not found: {dem_file}")
        return None
        
    print(f"    🔨 Building elevation mask from scratch...")
    
    try:
        with rasterio.open(dem_file) as dem_src:
            print(f"    📍 DEM CRS: {dem_src.crs}")
            print(f"    📏 DEM resolution: ~{dem_src.res[0]*111000:.0f}m")
            
            # Convert elevation to meters (DEM is typically in meters)
            max_elevation_m = max_elevation_ft * 0.3048
            min_elevation_m = min_elevation_ft * 0.3048 if min_elevation_ft is not None else None
            
            # Create a bounding box for Pacific Northwest (with US-Canada border)
            print(f"    📦 Clipping to Pacific Northwest bounds: {PACIFIC_NW_BOUNDS}")
            
            # Read DEM data for Pacific Northwest area
            from rasterio.windows import from_bounds
            window = from_bounds(*PACIFIC_NW_BOUNDS, dem_src.transform)
            
            # Clip window to DEM bounds
            window = window.intersection(
                rasterio.windows.Window(0, 0, dem_src.width, dem_src.height)
            )
            
            elevation_data = dem_src.read(1, window=window)
            window_transform = dem_src.window_transform(window)
            
            print(f"    🗺️  Read elevation data: {elevation_data.shape} ({elevation_data.nbytes/1024/1024:.1f} MB)")
            
            # Create mask for unsuitable elevation areas (True = outside range, to be removed)
            valid_data = ~np.isnan(elevation_data)
            
            if min_elevation_ft is not None:
                # Two-sided mask: too low OR too high
                # Handle negative elevations: when min_elevation_ft is 0, extend to -10ft to include coastal/river areas
                actual_min_elevation_m = min_elevation_m
                if min_elevation_ft == 0:
                    actual_min_elevation_m = -10 * 0.3048  # Extend to -10ft (-3.05m) for coastal areas
                    print(f"    🌊 Extending 0ft minimum to -10ft to include coastal/river areas")
                elif min_elevation_ft > 0 and min_elevation_ft <= 100:
                    # For low elevations (1-100ft), extend by 10ft to handle DEM negative values
                    actual_min_elevation_m = (min_elevation_ft - 10) * 0.3048
                    print(f"    🌊 Extending {min_elevation_ft}ft minimum to {min_elevation_ft-10}ft to handle negative DEM values")
                
                elevation_mask = ((elevation_data < actual_min_elevation_m) | (elevation_data > max_elevation_m)) & valid_data
                excluded_pixels = np.sum(elevation_mask)
                print(f"    📊 Excluded elevation pixels: {excluded_pixels:,} / {np.sum(valid_data):,} ({excluded_pixels/np.sum(valid_data)*100:.1f}%)")
                
                if excluded_pixels == 0:
                    print(f"    ✅ All areas within {min_elevation_ft:,} - {max_elevation_ft:,} feet - no masking needed")
                    return None
            else:
                # One-sided mask: too high only
                elevation_mask = (elevation_data > max_elevation_m) & valid_data
                excluded_pixels = np.sum(elevation_mask)
                print(f"    📊 High elevation pixels: {excluded_pixels:,} / {np.sum(valid_data):,} ({excluded_pixels/np.sum(valid_data)*100:.1f}%)")
                
                if excluded_pixels == 0:
                    print(f"    ✅ No areas above {max_elevation_ft:,} feet - no masking needed")
                    return None
            
            # Convert elevation mask to polygons
            print(f"    🔄 Converting elevation mask to polygons...")
            
            elevation_shapes = []
            for geom, value in rasterio.features.shapes(
                elevation_mask.astype(np.uint8), 
                mask=elevation_mask,
                transform=window_transform
            ):
                if value == 1:  # Excluded elevation areas
                    elevation_shapes.append(geom)
            
            if not elevation_shapes:
                print(f"    ✅ No unsuitable elevation polygons created")
                return None
            
            print(f"    📐 Created {len(elevation_shapes)} unsuitable elevation polygons")
            
            # Create GeoDataFrame of unsuitable elevation areas
            from shapely.geometry import shape
            if min_elevation_ft is not None:
                elevation_label = f'<{min_elevation_ft} or >{max_elevation_ft}'
            else:
                elevation_label = f'>{max_elevation_ft}'
                
            elevation_gdf = gpd.GeoDataFrame(
                [{'elevation_ft': elevation_label, 'mask_type': 'unsuitable_elevation'}] * len(elevation_shapes),
                geometry=[shape(geom) for geom in elevation_shapes],
                crs=dem_src.crs
            )
            
            # Reproject to WGS84 if needed
            if elevation_gdf.crs != "EPSG:4326":
                elevation_gdf = elevation_gdf.to_crs("EPSG:4326")
            
            # Merge overlapping elevation polygons
            print(f"    🔗 Merging overlapping unsuitable elevation areas...")
            merged_elevation = unary_union(elevation_gdf.geometry.tolist())
            
            # Save to cache
            if use_cache:
                save_cached_mask(merged_elevation, "elevation", elevation_ft=cache_key)
            
            print(f"    ✅ Elevation mask created successfully")
            return merged_elevation
            
    except Exception as e:
        print(f"    ❌ Error creating elevation mask: {e}")
        return None

def apply_elevation_constraints(distribution_gdf, elevation_mask, min_elevation_ft=None, max_elevation_ft=None):
    """Apply elevation constraints by subtracting unsuitable elevation areas"""
    if elevation_mask is None:
        print("  ⚠️  No elevation mask to apply")
        return distribution_gdf
    
    if min_elevation_ft is not None:
        print(f"  ✂️  Applying elevation constraints ({min_elevation_ft:,} - {max_elevation_ft:,} ft)...")
    else:
        print(f"  ✂️  Applying elevation constraints (>{max_elevation_ft:,} ft)...")
    
    try:
        # Get the distribution geometry
        distribution_geom = distribution_gdf.iloc[0].geometry
        
        # Log detailed area change
        constrained_geom = distribution_geom.difference(elevation_mask)
        modified_area_km2 = log_area_change("ELEVATION MASKING", distribution_geom, constrained_geom)
        
        # Calculate areas for backward compatibility
        if hasattr(distribution_geom, 'geoms'):
            original_area = sum(poly.area for poly in distribution_geom.geoms)
            original_count = len(distribution_geom.geoms)
        else:
            original_area = distribution_geom.area
            original_count = 1
        
        if hasattr(constrained_geom, 'geoms'):
            constrained_area = sum(poly.area for poly in constrained_geom.geoms)
            constrained_count = len(constrained_geom.geoms)
        else:
            constrained_area = constrained_geom.area
            constrained_count = 1
        
        area_reduction_pct = ((original_area - constrained_area) / original_area * 100) if original_area > 0 else 0
        
        # Update the GeoDataFrame
        constrained_gdf = distribution_gdf.copy()
        constrained_gdf.loc[0, 'geometry'] = constrained_geom
        constrained_gdf.loc[0, 'elevation_constrained'] = True
        constrained_gdf.loc[0, 'min_elevation_ft'] = min_elevation_ft
        constrained_gdf.loc[0, 'max_elevation_ft'] = max_elevation_ft
        constrained_gdf.loc[0, 'area_reduction_pct'] = round(area_reduction_pct, 1)
        constrained_gdf.loc[0, 'elevation_area_km2'] = round(modified_area_km2, 1)
        
        # Update polygon count
        constrained_gdf.loc[0, 'polygon_count'] = constrained_count
        
        print(f"    ✅ Elevation constraints applied successfully")
        return constrained_gdf
        
    except Exception as e:
        print(f"    ❌ Error applying elevation constraints: {e}")
        return distribution_gdf

def get_regional_tcc_mask_filename(threshold_pct, resolution_m=TCC_PROCESSING_RESOLUTION, min_area_sqkm=MIN_LOW_COVER_AREA_SQKM, version='processing'):
    """Generate filename for regional TCC mask with version preference"""
    base_filename = f"pnw_tcc_mask_{threshold_pct}pct_{resolution_m}m_{min_area_sqkm}sqkm"
    
    if version and version != 'original':
        filename = f"{base_filename}_{version}.geojson"
    else:
        filename = f"{base_filename}.geojson"
    
    return CACHE_DIR / filename

def find_best_tcc_version(threshold_pct, resolution_m=TCC_PROCESSING_RESOLUTION, min_area_sqkm=MIN_LOW_COVER_AREA_SQKM):
    """Find the best available version of a TCC mask configuration"""
    # Priority order: processing > analysis > display > original > qgis
    version_priority = ['processing', 'analysis', 'display', 'original', 'qgis']
    
    for version in version_priority:
        mask_file = get_regional_tcc_mask_filename(threshold_pct, resolution_m, min_area_sqkm, version)
        if mask_file.exists():
            return mask_file, version
    
    return None, None

def find_best_regional_tcc_mask(threshold_pct=None, resolution_m=None, min_area_sqkm=None):
    """Find the best available regional TCC mask that matches the target parameters"""
    if threshold_pct is None:
        threshold_pct = MIN_TREE_COVER_PCT
    if resolution_m is None:
        resolution_m = TCC_PROCESSING_RESOLUTION
    if min_area_sqkm is None:
        min_area_sqkm = MIN_LOW_COVER_AREA_SQKM
    
    print(f"🔍 Finding best regional TCC mask for: {threshold_pct}% threshold, {resolution_m}m resolution, {min_area_sqkm} sq km min area")
    
    # Get all possible TCC mask files
    mask_files = list(CACHE_DIR.glob("pnw_tcc_mask_*.geojson"))
    
    if not mask_files:
        print(f"    ❌ No regional TCC masks found")
        return None
    
    # Parse each mask file to extract parameters
    available_masks = []
    for mask_file in mask_files:
        try:
            # Parse filename: pnw_tcc_mask_5pct_30m_2.0sqkm.geojson
            name_parts = mask_file.stem.split('_')
            if len(name_parts) >= 6:
                file_threshold = int(name_parts[3].replace('pct', ''))
                file_resolution = int(name_parts[4].replace('m', ''))
                file_min_area = float(name_parts[5].replace('sqkm', ''))
                
                file_size_mb = mask_file.stat().st_size / (1024 * 1024)
                
                available_masks.append({
                    'file': mask_file,
                    'threshold': file_threshold,
                    'resolution': file_resolution,
                    'min_area': file_min_area,
                    'size_mb': file_size_mb
                })
        except:
            continue  # Skip files with unexpected naming
    
    if not available_masks:
        print(f"    ❌ No valid regional TCC mask files found")
        return None
    
    # Score each mask based on how well it matches the target
    scored_masks = []
    
    for mask in available_masks:
        score = 0
        
        # Threshold scoring (most important - 100 points max)
        threshold_diff = abs(mask['threshold'] - threshold_pct)
        if threshold_diff == 0:
            score += 100  # Perfect match
        elif threshold_diff <= 2:
            score += 80   # Very close
        elif threshold_diff <= 5:
            score += 60   # Close
        else:
            score += max(0, 40 - threshold_diff)  # Penalize large differences
        
        # Resolution scoring (30 points max)
        res_diff = abs(mask['resolution'] - resolution_m)
        if res_diff == 0:
            score += 30  # Perfect match
        elif res_diff <= 30:
            score += 20  # Close
        else:
            score += max(0, 10 - res_diff/10)  # Penalize large differences
        
        # Min area scoring (20 points max)
        area_diff = abs(mask['min_area'] - min_area_sqkm)
        if area_diff == 0:
            score += 20  # Perfect match
        elif area_diff <= 1:
            score += 15  # Close
        else:
            score += max(0, 10 - area_diff)  # Penalize large differences
        
        scored_masks.append({
            'mask': mask,
            'score': score,
            'threshold_diff': threshold_diff,
            'res_diff': res_diff,
            'area_diff': area_diff
        })
    
    # Sort by score (highest first)
    scored_masks.sort(key=lambda x: x['score'], reverse=True)
    
    best_match = scored_masks[0]
    best_mask = best_match['mask']
    
    print(f"    ✅ Best match: {best_mask['file'].name}")
    print(f"       Threshold: {best_mask['threshold']}% (diff: {best_match['threshold_diff']})")
    print(f"       Resolution: {best_mask['resolution']}m (diff: {best_match['res_diff']})")
    print(f"       Min area: {best_mask['min_area']} sq km (diff: {best_match['area_diff']:.1f})")
    print(f"       Match score: {best_match['score']}/150")
    
    # Warn if match is poor
    if best_match['score'] < 80:
        print(f"    ⚠️  WARNING: Poor match (score < 80)")
        print(f"       Consider creating a better mask:")
        print(f"       python scripts/create_regional_tcc_mask.py --threshold {threshold_pct} --resolution {resolution_m} --min-area {min_area_sqkm}")
    
    return best_mask['file']

def auto_create_regional_tcc_mask(threshold_pct, resolution_m=None, min_area_sqkm=None):
    """Automatically create a regional TCC mask if a suitable one doesn't exist"""
    if resolution_m is None:
        resolution_m = TCC_PROCESSING_RESOLUTION
    if min_area_sqkm is None:
        min_area_sqkm = MIN_LOW_COVER_AREA_SQKM
    
    # Check if exact match exists
    cache_file = get_regional_tcc_mask_filename(threshold_pct, resolution_m, min_area_sqkm)
    if cache_file.exists():
        return cache_file
    
    print(f"🔨 Auto-creating regional TCC mask: {cache_file.name}")
    print(f"    This may take 10-20 minutes...")
    
    try:
        import subprocess
        import sys
        
        # Call the regional TCC mask script
        cmd = [
            sys.executable, 
            "scripts/create_regional_tcc_mask.py",
            "--threshold", str(threshold_pct),
            "--resolution", str(resolution_m),
            "--min-area", str(min_area_sqkm),
            "--non-interactive"
        ]
        
        print(f"    Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode == 0 and cache_file.exists():
            print(f"    ✅ Auto-created regional TCC mask successfully")
            return cache_file
        else:
            print(f"    ❌ Failed to auto-create mask")
            if result.stderr:
                print(f"    Error: {result.stderr}")
            return None
            
    except Exception as e:
        print(f"    ❌ Auto-creation failed: {e}")
        return None

def load_regional_tcc_mask(threshold_pct=None, mask_file=None):
    """Load a pre-created regional TCC mask with intelligent fallback options"""
    if threshold_pct is None:
        threshold_pct = MIN_TREE_COVER_PCT
    
    if mask_file is not None:
        # Use specified mask file
        mask_path = Path(mask_file)
        if not mask_path.exists():
            print(f"❌ Specified regional TCC mask not found: {mask_path}")
            return None
    else:
        # Try to find best version for the configuration
        best_file, version_used = find_best_tcc_version(threshold_pct, TCC_PROCESSING_RESOLUTION, MIN_LOW_COVER_AREA_SQKM)
        
        if best_file is not None:
            mask_path = best_file
            print(f"✅ Found {version_used} version: {mask_path.name}")
        else:
            print(f"⚠️  No TCC mask found for {threshold_pct}% threshold")
            
            # Try to find the best available alternative configuration
            print(f"🔍 Looking for alternative regional TCC masks...")
            best_match_file = find_best_regional_tcc_mask(threshold_pct, TCC_PROCESSING_RESOLUTION, MIN_LOW_COVER_AREA_SQKM)
            
            if best_match_file is not None:
                mask_path = best_match_file
                print(f"✅ Using alternative mask: {mask_path.name}")
            else:
                # Try to auto-create the needed mask
                print(f"🔨 No suitable alternatives found - attempting auto-creation...")
                auto_created_file = auto_create_regional_tcc_mask(threshold_pct, TCC_PROCESSING_RESOLUTION, MIN_LOW_COVER_AREA_SQKM)
                
                if auto_created_file is not None:
                    mask_path = auto_created_file
                    print(f"✅ Auto-created mask: {mask_path.name}")
                else:
                    print(f"❌ Could not find or create suitable regional TCC mask")
                    print(f"💡 Manual creation:")
                    print(f"   python scripts/create_regional_tcc_mask.py --threshold {threshold_pct} --resolution {TCC_PROCESSING_RESOLUTION} --min-area {MIN_LOW_COVER_AREA_SQKM}")
                    return None
    
    print(f"📦 Loading regional TCC mask: {mask_path.name}")
    file_size_mb = mask_path.stat().st_size / (1024 * 1024)
    
    try:
        gdf = gpd.read_file(mask_path)
        if len(gdf) == 0:
            print(f"⚠️  Regional mask file is empty")
            return None
        
        mask_geometry = gdf.iloc[0].geometry
        if mask_geometry is None or mask_geometry.is_empty:
            print(f"⚠️  Regional mask has invalid geometry")
            return None
        
        # Display mask metadata if available
        row = gdf.iloc[0]
        if 'threshold_pct' in row:
            print(f"    🌳 Threshold: <{row['threshold_pct']}% tree cover")
        if 'resolution_m' in row:
            print(f"    📏 Resolution: {row['resolution_m']}m")
        if 'min_area_sqkm' in row:
            print(f"    📐 Min area: {row['min_area_sqkm']} sq km")
        if 'created_date' in row:
            created_date = row['created_date'][:10]  # Just the date part
            print(f"    📅 Created: {created_date}")
        
        print(f"✅ Loaded regional TCC mask ({file_size_mb:.1f} MB)")
        return mask_geometry
        
    except Exception as e:
        print(f"❌ Error loading regional TCC mask: {e}")
        return None

def create_tree_cover_mask(bounds, distribution_geom=None, tcc_file=TCC_FILE, min_tree_cover_pct=None, use_cache=True):
    if min_tree_cover_pct is None:
        min_tree_cover_pct = MIN_TREE_COVER_PCT
    """Create a high-resolution mask of areas with insufficient tree canopy cover"""
    print(f"  🌳 Creating HIGH-RESOLUTION tree cover mask for areas with <{min_tree_cover_pct}% canopy cover...")
    
    # Check cache first with resolution and area parameters
    cache_params = {
        'tcc_pct': min_tree_cover_pct,
        'resolution': TCC_PROCESSING_RESOLUTION,
        'min_area_sqkm': MIN_LOW_COVER_AREA_SQKM
    }
    
    if use_cache:
        cached_mask = load_cached_mask("tcc", **cache_params)
        if cached_mask is not None:
            return cached_mask
    
    if not HAS_RASTERIO:
        print("  ⚠️  Skipping tree cover masking - rasterio not available")
        return None
        
    if not Path(tcc_file).exists():
        print(f"  ⚠️  Skipping tree cover masking - TCC file not found: {tcc_file}")
        return None
    
    print(f"    🔨 Building HIGH-RESOLUTION tree cover mask from scratch...")
    print(f"    🎯 High-res settings: {TCC_PROCESSING_RESOLUTION}m resolution + ≥{MIN_LOW_COVER_AREA_SQKM} sq km areas + fine-grained detail")
    print(f"    📊 Target: Capture small patches of {min_tree_cover_pct}%+ tree cover with maximum spatial detail")
    
    try:
        with rasterio.open(tcc_file) as tcc_src:
            print(f"    📍 TCC CRS: {tcc_src.crs}")
            print(f"    📏 TCC original resolution: {tcc_src.res[0]:.0f}m")
            
            # OPTIMIZATION 1: Use species distribution bounds instead of entire Washington
            if distribution_geom is not None:
                # Get bounds of the species distribution (more focused area)
                species_bounds = distribution_geom.bounds  # (minx, miny, maxx, maxy) in WGS84
                print(f"    🎯 Using species distribution bounds instead of full Washington")
                print(f"    📦 Species bounds: {[round(x, 3) for x in species_bounds]}")
                
                # Add small buffer to ensure we don't clip too tightly
                buffer = 0.01  # ~1km buffer
                clipping_bounds = [
                    species_bounds[0] - buffer,  # minx
                    species_bounds[1] - buffer,  # miny  
                    species_bounds[2] + buffer,  # maxx
                    species_bounds[3] + buffer   # maxy
                ]
            else:
                # Fallback to Pacific Northwest bounds
                clipping_bounds = PACIFIC_NW_BOUNDS
                print(f"    📦 Using Pacific Northwest bounds (no species distribution provided)")
            
            # Reproject bounds to TCC CRS
            from rasterio.warp import transform_bounds
            clipping_bounds_tcc = transform_bounds("EPSG:4326", tcc_src.crs, *clipping_bounds)
            
            print(f"    📦 Clipping bounds (TCC CRS): {[round(x) for x in clipping_bounds_tcc]}")
            
            # HIGH-RESOLUTION PROCESSING: Minimize downsampling for maximum detail
            from rasterio.windows import from_bounds
            from rasterio.enums import Resampling
            
            # Calculate minimal downsampling factor for high-resolution processing
            original_res = tcc_src.res[0]  # meters
            if TCC_HIGH_RES_MODE and TCC_PROCESSING_RESOLUTION <= 30:
                # Use native resolution or minimal downsampling for high-res mode
                downsample_factor = max(1, int(TCC_PROCESSING_RESOLUTION / original_res))
                print(f"    🔍 HIGH-RES MODE: {original_res:.0f}m → {original_res * downsample_factor:.0f}m (factor: {downsample_factor})")
                print(f"    💡 Preserving maximum spatial detail for fine-grained tree cover analysis")
            else:
                # Standard processing
                downsample_factor = max(1, int(TCC_PROCESSING_RESOLUTION / original_res))
                print(f"    📉 Standard processing: {original_res:.0f}m → {original_res * downsample_factor:.0f}m (factor: {downsample_factor})")
            
            # Read downsampled data
            window = from_bounds(*clipping_bounds_tcc, tcc_src.transform)
            window = window.intersection(
                rasterio.windows.Window(0, 0, tcc_src.width, tcc_src.height)
            )
            
            if window.width <= 0 or window.height <= 0:
                print(f"    ⚠️  No overlap between bounds and TCC data")
                return None
            
            # Read with downsampling
            out_height = max(1, int(window.height // downsample_factor))
            out_width = max(1, int(window.width // downsample_factor))
            
            tcc_data = tcc_src.read(
                1, 
                window=window,
                out_shape=(out_height, out_width),
                resampling=Resampling.average
            )
            
            # Update transform for downsampled data
            window_transform = tcc_src.window_transform(window)
            # Scale transform for downsampling
            from rasterio.transform import Affine
            window_transform = Affine(
                window_transform.a * downsample_factor,  # x pixel size
                window_transform.b,                      # rotation
                window_transform.c,                      # x offset
                window_transform.d,                      # rotation  
                window_transform.e * downsample_factor,  # y pixel size
                window_transform.f                       # y offset
            )
            
            original_mb = (window.height * window.width * 1) / (1024 * 1024)
            actual_mb = tcc_data.nbytes / (1024 * 1024)
            reduction_pct = ((original_mb - actual_mb) / original_mb * 100) if original_mb > 0 else 0
            
            print(f"    🗺️  Read TCC data: {tcc_data.shape} ({actual_mb:.1f} MB, {reduction_pct:.1f}% reduction)")
            
            # Create mask for low tree cover areas
            valid_data_mask = (tcc_data >= 0) & (tcc_data <= 100)
            low_tree_cover_mask = (tcc_data < min_tree_cover_pct) & valid_data_mask
            
            low_cover_pixels = np.sum(low_tree_cover_mask)
            total_pixels = np.sum(valid_data_mask)
            low_cover_pct = (low_cover_pixels / total_pixels * 100) if total_pixels > 0 else 0
            
            print(f"    📊 Low tree cover pixels: {low_cover_pixels:,} / {total_pixels:,} ({low_cover_pct:.1f}%)")
            
            if low_cover_pixels == 0:
                print(f"    ✅ All areas have ≥{min_tree_cover_pct}% tree cover - no masking needed")
                return None
            
            # HIGH-RESOLUTION MORPHOLOGICAL PROCESSING: Fine-tuned for detail preservation
            if HAS_SCIPY and TCC_MORPHOLOGICAL_SMOOTHING:
                print(f"    🎯 Applying fine-tuned morphological processing for high-resolution data...")
                import time
                morph_start = time.time()
                
                if TCC_HIGH_RES_MODE:
                    # More aggressive processing to reduce small artifacts
                    structure = np.ones((5, 5))  # Larger structure to eliminate more artifacts
                    # More aggressive processing to reduce small artifacts
                    print(f"    🔄 Step 1/3: Binary closing (noise reduction)...")
                    low_tree_cover_mask_smoothed = ndimage.binary_closing(
                        low_tree_cover_mask, structure=structure, iterations=2
                    )
                    print(f"    🔄 Step 2/3: Binary opening (small artifact removal)...")
                    # More aggressive opening to remove small artifacts
                    low_tree_cover_mask_smoothed = ndimage.binary_opening(
                        low_tree_cover_mask_smoothed, structure=structure, iterations=2
                    )
                    print(f"    🔄 Step 3/3: Additional closing (final smoothing)...")
                    # Final closing to ensure connectivity
                    low_tree_cover_mask_smoothed = ndimage.binary_closing(
                        low_tree_cover_mask_smoothed, structure=structure, iterations=1
                    )
                    morph_elapsed = time.time() - morph_start
                    print(f"    ✅ Aggressive morphological processing applied - reduced artifacts ({morph_elapsed:.1f}s)")
                else:
                    # Standard processing for coarser resolutions
                    structure = np.ones((5, 5))
                    print(f"    🔄 Step 1/2: Binary closing (iterations=3)...")
                    low_tree_cover_mask_smoothed = ndimage.binary_closing(
                        low_tree_cover_mask, structure=structure, iterations=3
                    )
                    print(f"    🔄 Step 2/2: Binary opening (iterations=2)...")
                    low_tree_cover_mask_smoothed = ndimage.binary_opening(
                        low_tree_cover_mask_smoothed, structure=structure, iterations=2
                    )
                    morph_elapsed = time.time() - morph_start
                    print(f"    ✅ Standard morphological processing applied ({morph_elapsed:.1f}s)")
            else:
                low_tree_cover_mask_smoothed = low_tree_cover_mask
                if not HAS_SCIPY:
                    print(f"    ⚠️  Scipy not available - skipping morphological processing")
                else:
                    print(f"    ⚠️  Morphological smoothing disabled - using raw mask data")
            
            # Convert to polygons with size filtering
            print(f"    🔄 Converting to polygons with size filtering...")
            
            low_cover_shapes = []
            pixel_area_sqm = abs(window_transform.a * window_transform.e)  # pixel area in square meters
            min_pixels = int((MIN_LOW_COVER_AREA_SQKM * 1_000_000) / pixel_area_sqm)  # Convert sq km to pixels
            
            print(f"    📏 Pixel area: {pixel_area_sqm:.0f} sq m, Min area: {MIN_LOW_COVER_AREA_SQKM} sq km = {min_pixels} pixels")
            
            shape_count = 0
            processed_count = 0
            
            # Add progress tracking
            import time
            start_time = time.time()
            last_progress_time = start_time
            
            print(f"    🎯 Starting polygon extraction with progress tracking...")
            
            for geom, value in rasterio.features.shapes(
                low_tree_cover_mask_smoothed.astype(np.uint8), 
                mask=low_tree_cover_mask_smoothed,
                transform=window_transform
            ):
                processed_count += 1
                
                # Progress reporting every 5 seconds
                current_time = time.time()
                if current_time - last_progress_time >= 5.0:
                    elapsed = current_time - start_time
                    print(f"    ⏱️  Progress: {processed_count:,} polygons processed, {shape_count:,} kept in {elapsed:.1f}s")
                    last_progress_time = current_time
                
                if value == 1:  # Low tree cover areas
                    # HIGH-RESOLUTION SIZE FILTERING: Smaller minimum areas to capture fine detail
                    from shapely.geometry import shape as shapely_shape
                    polygon = shapely_shape(geom)
                    area_sqm = polygon.area  # This is approximate in projected coordinates
                    area_sqkm = area_sqm / 1_000_000
                    
                    if area_sqkm >= MIN_LOW_COVER_AREA_SQKM:
                        low_cover_shapes.append(geom)
                        shape_count += 1
                    
                    # Higher safety limit for high-resolution processing
                    safety_limit = 50000 if TCC_HIGH_RES_MODE else 10000
                    if shape_count >= safety_limit:
                        print(f"    ⚠️  Reached {safety_limit:,} polygons - stopping for performance")
                        print(f"    💡 Consider increasing MIN_LOW_COVER_AREA_SQKM if too many small polygons")
                        break
            
            total_elapsed = time.time() - start_time
            print(f"    ✅ Polygon extraction complete: {processed_count:,} total processed, {shape_count:,} kept in {total_elapsed:.1f}s")
            
            if not low_cover_shapes:
                print(f"    ✅ No large low tree cover areas found (all < {MIN_LOW_COVER_AREA_SQKM} sq km)")
                return None
            
            print(f"    📐 Created {len(low_cover_shapes):,} low tree cover polygons (≥{MIN_LOW_COVER_AREA_SQKM} sq km)")
            print(f"    🎯 High-resolution processing captured fine-grained spatial detail")
            
            # Create GeoDataFrame with high-resolution metadata
            from shapely.geometry import shape
            geometries = [shape(geom) for geom in low_cover_shapes]
            
            low_cover_gdf = gpd.GeoDataFrame(
                [{
                    'tree_cover_pct': f'<{min_tree_cover_pct}', 
                    'mask_type': 'low_tree_cover', 
                    'min_area_sqkm': MIN_LOW_COVER_AREA_SQKM,
                    'resolution_m': TCC_PROCESSING_RESOLUTION,
                    'high_res_mode': TCC_HIGH_RES_MODE
                }] * len(geometries),
                geometry=geometries,
                crs=tcc_src.crs
            )
            
            # Reproject to WGS84
            if low_cover_gdf.crs != "EPSG:4326":
                print(f"    🔄 Reprojecting to WGS84...")
                low_cover_gdf = low_cover_gdf.to_crs("EPSG:4326")
            
            # Simple merge since we have far fewer polygons now
            print(f"    🔗 Merging {len(low_cover_gdf):,} low tree cover areas...")
            try:
                merged_low_cover = unary_union(low_cover_gdf.geometry.tolist())
            except Exception as e:
                print(f"    ⚠️  Geometry error during merge: {e}")
                print(f"    🔧 Attempting repair with buffer(0) fix...")
                try:
                    # Fix invalid geometries
                    valid_geometries = []
                    for geom in low_cover_gdf.geometry:
                        if geom.is_valid:
                            valid_geometries.append(geom)
                        else:
                            try:
                                fixed_geom = geom.buffer(0)
                                if fixed_geom.is_valid:
                                    valid_geometries.append(fixed_geom)
                            except:
                                continue
                    
                    if valid_geometries:
                        merged_low_cover = unary_union(valid_geometries)
                        print(f"    ✅ Repair successful: {len(valid_geometries)}/{len(low_cover_gdf)} geometries kept")
                    else:
                        print(f"    ❌ Could not repair geometries - skipping tree cover masking")
                        return None
                except Exception as e2:
                    print(f"    ❌ Repair failed: {e2}")
                    return None
            
            # Save to cache with resolution parameters
            if use_cache:
                save_cached_mask(merged_low_cover, "tcc", **cache_params)
            
            print(f"    ✅ HIGH-RESOLUTION tree cover mask created successfully!")
            print(f"    📊 Result: {min_tree_cover_pct}% threshold at {TCC_PROCESSING_RESOLUTION}m resolution with {MIN_LOW_COVER_AREA_SQKM} sq km minimum area")
            return merged_low_cover
            
    except Exception as e:
        print(f"    ❌ Error creating tree cover mask: {e}")
        import traceback
        traceback.print_exc()
        return None

def apply_tree_cover_constraints(distribution_gdf, tree_cover_mask):
    """Apply tree cover constraints by subtracting low tree cover areas"""
    if tree_cover_mask is None:
        print("  ⚠️  No tree cover mask to apply")
        return distribution_gdf
    
    print(f"  ✂️  Applying tree cover constraints...")
    
    try:
        # Get the distribution geometry
        distribution_geom = distribution_gdf.iloc[0].geometry
        
        # Log detailed area change
        constrained_geom = distribution_geom.difference(tree_cover_mask)
        modified_area_km2 = log_area_change("TREE COVER MASKING", distribution_geom, constrained_geom)
        
        # Calculate areas for backward compatibility
        if hasattr(distribution_geom, 'geoms'):
            original_area = sum(poly.area for poly in distribution_geom.geoms)
            original_count = len(distribution_geom.geoms)
        else:
            original_area = distribution_geom.area
            original_count = 1
        
        if hasattr(constrained_geom, 'geoms'):
            constrained_area = sum(poly.area for poly in constrained_geom.geoms)
            constrained_count = len(constrained_geom.geoms)
        else:
            constrained_area = constrained_geom.area
            constrained_count = 1
        
        area_reduction_pct = ((original_area - constrained_area) / original_area * 100) if original_area > 0 else 0
        
        # Update the GeoDataFrame
        constrained_gdf = distribution_gdf.copy()
        constrained_gdf.loc[0, 'geometry'] = constrained_geom
        constrained_gdf.loc[0, 'tree_cover_constrained'] = True
        constrained_gdf.loc[0, 'min_tree_cover_pct'] = MIN_TREE_COVER_PCT
        
        # Add or update area reduction percentage
        elevation_reduction = constrained_gdf.loc[0].get('area_reduction_pct', 0)
        constrained_gdf.loc[0, 'tcc_area_reduction_pct'] = round(area_reduction_pct, 1)
        constrained_gdf.loc[0, 'tcc_area_km2'] = round(modified_area_km2, 1)
        
        # Fix total reduction calculation
        original_area_km2 = constrained_gdf.loc[0].get('original_area_km2', calculate_area_km2(distribution_geom))
        total_reduction_pct = ((original_area_km2 - modified_area_km2) / original_area_km2 * 100) if original_area_km2 > 0 else 0
        constrained_gdf.loc[0, 'total_area_reduction_pct'] = round(total_reduction_pct, 1)
        
        # Update polygon count
        constrained_gdf.loc[0, 'polygon_count'] = constrained_count
        
        print(f"    ✅ Tree cover constraints applied successfully")
        return constrained_gdf
        
    except Exception as e:
        print(f"    ❌ Error applying tree cover constraints: {e}")
        return distribution_gdf

def create_hydrography_mask(bounds, use_cache=True):
    """Create a mask of water bodies from available hydrography data"""
    if not APPLY_HYDROGRAPHY_MASK:
        return None
        
    print(f"  🌊 Creating hydrography mask from available water data...")
    
    # Check cache first
    if use_cache:
        cached_mask = load_cached_mask("hydrography")
        if cached_mask is not None:
            return cached_mask
    
    print(f"    🔨 Building hydrography mask from scratch...")
    
    # Initialize collections for multi-state processing
    all_water_features = []
    successful_sources = []
    
    # Try each hydrography source in order
    for source in HYDROGRAPHY_SOURCES:
        print(f"    Trying: {source['description']}")
        
        if source['type'] == 'built-in':
            # Use simple built-in water polygons
            return create_builtin_water_polygons()
        
        elif source['type'] == 'raster':
            file_path = source['file']
            
            if not Path(file_path).exists():
                print(f"      ⚠️  File not found: {file_path}")
                continue
            
            if not HAS_RASTERIO:
                print(f"      ⚠️  Rasterio not available for raster processing")
                continue
            
            try:
                # Load raster water data (e.g., NHDPlus stream network)
                print(f"      📁 Loading raster: {file_path}")
                
                with rasterio.open(file_path) as src:
                    print(f"      📍 Raster CRS: {src.crs}")
                    print(f"      📏 Raster resolution: {src.res[0]:.1f}m")
                    
                    # Create a bounding box for Pacific Northwest (with US-Canada border)
                    # Reproject PNW bounds to raster CRS
                    from rasterio.warp import transform_bounds
                    pnw_bounds_raster = transform_bounds("EPSG:4326", src.crs, *PACIFIC_NW_BOUNDS)
                    
                    print(f"      📦 Clipping to Pacific Northwest bounds...")
                    
                    # Read raster data for Pacific Northwest area
                    from rasterio.windows import from_bounds
                    window = from_bounds(*pnw_bounds_raster, src.transform)
                    
                    # Clip window to raster bounds
                    window = window.intersection(
                        rasterio.windows.Window(0, 0, src.width, src.height)
                    )
                    
                    if window.width <= 0 or window.height <= 0:
                        print(f"      ⚠️  No overlap between Pacific Northwest bounds and raster data")
                        continue
                    
                    water_data = src.read(1, window=window)
                    window_transform = src.window_transform(window)
                    
                    print(f"      🗺️  Read raster data: {water_data.shape} ({water_data.nbytes/1024/1024:.1f} MB)")
                    
                    # For NHDPlus stream network: value 1 = water, others = no water
                    # Create mask for water areas (True = water, to be removed)
                    water_mask = (water_data == 1)  # Adjust this based on your raster values
                    
                    water_pixels = np.sum(water_mask)
                    total_pixels = np.sum(water_data != src.nodata)
                    water_pct = (water_pixels / total_pixels * 100) if total_pixels > 0 else 0
                    
                    print(f"      💧 Water pixels: {water_pixels:,} / {total_pixels:,} ({water_pct:.1f}%)")
                    
                    if water_pixels == 0:
                        print(f"      ⚠️  No water pixels found in data")
                        continue
                    
                    # Convert water mask to polygons
                    print(f"      🔄 Converting water mask to polygons...")
                    
                    water_shapes = []
                    for geom, value in rasterio.features.shapes(
                        water_mask.astype(np.uint8), 
                        mask=water_mask,
                        transform=window_transform
                    ):
                        if value == 1:  # Water areas
                            water_shapes.append(geom)
                    
                    if not water_shapes:
                        print(f"      ⚠️  No water polygons created")
                        continue
                    
                    print(f"      📐 Created {len(water_shapes)} water polygons")
                    
                    # Create GeoDataFrame of water areas
                    from shapely.geometry import shape
                    water_gdf = gpd.GeoDataFrame(
                        [{'water_source': 'stream_network', 'mask_type': 'hydrography'}] * len(water_shapes),
                        geometry=[shape(geom) for geom in water_shapes],
                        crs=src.crs
                    )
                    
                    # Reproject to WGS84
                    if water_gdf.crs != "EPSG:4326":
                        water_gdf = water_gdf.to_crs("EPSG:4326")
                    
                    # Merge overlapping water polygons
                    print(f"      🔗 Merging water bodies...")
                    merged_water = unary_union(water_gdf.geometry.tolist())
                    
                    print(f"      ✅ Hydrography mask created from {source['description']}")
                    return merged_water
                    
            except Exception as e:
                print(f"      ❌ Error processing raster {file_path}: {e}")
                continue
        
        elif source['type'] == 'vector':
            file_path = source['file']
            
            if not Path(file_path).exists():
                print(f"      ⚠️  File not found: {file_path}")
                continue
            
            try:
                # Load vector water data
                print(f"      📁 Loading: {file_path}")
                water_gdf = gpd.read_file(file_path)
                
                print(f"      📊 Loaded {len(water_gdf):,} water features")
                print(f"      📍 CRS: {water_gdf.crs}")
                
                # Filter for major water features only
                print(f"      🔍 Filtering for major water features...")
                
                original_count = len(water_gdf)
                
                # Filter 1: Size-based filtering
                if 'areasqkm' in water_gdf.columns:
                    # Keep features larger than 0.1 sq km (reasonable lakes/ponds)
                    size_mask = water_gdf['areasqkm'] > 0.1
                    print(f"        💧 Size filter: {size_mask.sum():,} features > 0.1 sq km")
                else:
                    size_mask = pd.Series([True] * len(water_gdf))
                    print(f"        ⚠️  No area column - skipping size filter")
                
                # Filter 2: Feature type filtering
                if 'ftype' in water_gdf.columns:
                    # Keep major water body types:
                    # 390 = Lake/Pond, 466 = Reservoir, 493 = Estuary 
                    # Exclude: 336 = Canal/Ditch, 460 = StreamRiver, etc.
                    major_types = [390, 466, 493]  # Lake/Pond, Reservoir, Estuary
                    type_mask = water_gdf['ftype'].isin(major_types)
                    print(f"        🏞️  Type filter: {type_mask.sum():,} major water body types (Lake/Pond/Reservoir/Estuary)")
                else:
                    type_mask = pd.Series([True] * len(water_gdf))
                    print(f"        ⚠️  No ftype column - skipping type filter")
                
                # Filter 3: Named important features (keep major sounds/bays even if smaller)
                if 'gnis_name' in water_gdf.columns:
                    # Important named water bodies to always include
                    important_names = [
                        'Puget Sound', 'Hood Canal', 'Strait of Juan de Fuca',
                        'Columbia River', 'Snake River', 'Lake Washington', 
                        'Lake Chelan', 'Lake Crescent', 'Franklin D. Roosevelt Lake',
                        'Lake Sammamish', 'Green Lake', 'Lake Union', 'Ross Lake',
                        'Diablo Lake', 'Gorge Lake', 'Lake Cushman', 'Alder Lake'
                    ]
                    
                    # Create case-insensitive matching
                    name_mask = pd.Series([False] * len(water_gdf))
                    if not water_gdf['gnis_name'].isna().all():
                        for name in important_names:
                            name_mask |= water_gdf['gnis_name'].str.contains(
                                name, case=False, na=False
                            )
                    
                    print(f"        🏷️  Named features: {name_mask.sum():,} important named water bodies")
                else:
                    name_mask = pd.Series([False] * len(water_gdf))
                    print(f"        ⚠️  No gnis_name column - skipping named features")
                
                # Combine filters: (size AND type) OR important names
                combined_mask = (size_mask & type_mask) | name_mask
                water_gdf = water_gdf[combined_mask]
                
                filtered_count = len(water_gdf)
                reduction_pct = ((original_count - filtered_count) / original_count) * 100
                
                print(f"      ✨ Filtered: {original_count:,} → {filtered_count:,} features ({reduction_pct:.1f}% reduction)")
                print(f"      🎯 Focus: Major lakes, reservoirs, estuaries, and important named water bodies")
                
                # Reproject to WGS84 if needed
                if water_gdf.crs != "EPSG:4326":
                    water_gdf = water_gdf.to_crs("EPSG:4326")
                
                # Clip to Pacific Northwest bounds
                from shapely.geometry import box
                pnw_bbox = box(*PACIFIC_NW_BOUNDS)
                
                # Filter to features that intersect Pacific Northwest
                water_in_pnw = water_gdf[water_gdf.intersects(pnw_bbox)]
                print(f"      🗺️  {len(water_in_pnw):,} major water features intersect Pacific Northwest")
                
                if len(water_in_pnw) == 0:
                    print(f"      ⚠️  No major water features in Pacific Northwest bounds")
                    continue
                
                # Add to collection instead of returning immediately
                all_water_features.append(water_in_pnw)
                successful_sources.append(source['description'])
                
                print(f"      ✅ Collected water features from {source['description']} (major features only)")
                continue  # Continue to next source instead of returning
                
            except Exception as e:
                print(f"      ❌ Error loading {file_path}: {e}")
                continue
    
    # Process collected water features from all sources
    if not all_water_features:
        print(f"    ⚠️  No hydrography data sources available")
        return None
    
    print(f"  📊 Successfully collected water features from {len(successful_sources)} sources:")
    for source_desc in successful_sources:
        print(f"    ✅ {source_desc}")
    
    # Combine all water features
    print(f"  🔗 Combining water features from all Pacific Northwest states...")
    combined_water_gdf = pd.concat(all_water_features, ignore_index=True)
    
    print(f"  📐 Total water features: {len(combined_water_gdf):,}")
    
    # Merge all water bodies into a single geometry
    print(f"    🔗 Merging all Pacific Northwest water bodies...")
    merged_water = unary_union(combined_water_gdf.geometry.tolist())
    
    # Save to cache
    if use_cache:
        save_cached_mask(merged_water, "hydrography")
    
    print(f"    ✅ Pacific Northwest hydrography mask created successfully!")
    return merged_water

def create_builtin_water_polygons():
    """Create simple built-in water body polygons for Washington state"""
    print(f"      🔧 Using built-in water polygons...")
    
    from shapely.geometry import Polygon
    
    # Major water bodies in Washington (approximate)
    water_bodies = [
        # Puget Sound (simplified)
        [(-123.2, 47.0), (-122.2, 47.0), (-122.2, 48.5), (-123.2, 48.5)],
        
        # Lake Washington area
        [(-122.3, 47.5), (-122.2, 47.5), (-122.2, 47.7), (-122.3, 47.7)],
        
        # Columbia River - southern border
        [(-124.0, 45.5), (-121.0, 45.5), (-121.0, 46.3), (-124.0, 46.3)],
        [(-121.0, 45.5), (-117.0, 45.5), (-117.0, 46.1), (-121.0, 46.1)],
        
        # Lake Chelan
        [(-120.4, 47.8), (-120.1, 47.8), (-120.1, 48.1), (-120.4, 48.1)],
        
        # Franklin D. Roosevelt Lake (Columbia River reservoir)
        [(-118.5, 47.5), (-118.0, 47.5), (-118.0, 49.0), (-118.5, 49.0)],
        
        # Strait of Juan de Fuca
        [(-124.8, 48.3), (-123.0, 48.3), (-123.0, 48.6), (-124.8, 48.6)],
    ]
    
    water_polygons = []
    for coords in water_bodies:
        if len(coords) >= 3:
            try:
                water_polygons.append(Polygon(coords))
            except Exception:
                continue
    
    if not water_polygons:
        print(f"      ⚠️  No valid built-in water polygons created")
        return None
    
    # Merge all water polygons
    merged_water = unary_union(water_polygons)
    
    print(f"      📐 Created {len(water_polygons)} built-in water body polygons")
    return merged_water

def apply_hydrography_constraints(distribution_gdf, hydrography_mask):
    """Apply hydrography constraints by subtracting water bodies"""
    if hydrography_mask is None:
        print("  ⚠️  No hydrography mask to apply")
        return distribution_gdf
    
    print(f"  ✂️  Applying hydrography constraints...")
    
    try:
        # Get the distribution geometry
        distribution_geom = distribution_gdf.iloc[0].geometry
        
        # Calculate areas before masking
        if hasattr(distribution_geom, 'geoms'):
            original_area = sum(poly.area for poly in distribution_geom.geoms)
            original_count = len(distribution_geom.geoms)
        else:
            original_area = distribution_geom.area
            original_count = 1
        
        print(f"    📏 Before hydrography masking: {original_count} polygons, {original_area:.4f} sq degrees")
        
        # Subtract water bodies from distribution
        constrained_geom = distribution_geom.difference(hydrography_mask)
        
        # Calculate areas after masking
        if hasattr(constrained_geom, 'geoms'):
            constrained_area = sum(poly.area for poly in constrained_geom.geoms)
            constrained_count = len(constrained_geom.geoms)
        else:
            constrained_area = constrained_geom.area
            constrained_count = 1
        
        area_reduction_pct = ((original_area - constrained_area) / original_area * 100) if original_area > 0 else 0
        
        print(f"    📏 After hydrography masking: {constrained_count} polygons, {constrained_area:.4f} sq degrees")
        print(f"    📉 Hydrography area reduction: {area_reduction_pct:.1f}%")
        
        # Update the GeoDataFrame
        constrained_gdf = distribution_gdf.copy()
        constrained_gdf.loc[0, 'geometry'] = constrained_geom
        constrained_gdf.loc[0, 'hydrography_constrained'] = True
        
        # Add or update area reduction percentage
        total_reduction = constrained_gdf.loc[0].get('total_area_reduction_pct', 0) + area_reduction_pct
        constrained_gdf.loc[0, 'hydro_area_reduction_pct'] = round(area_reduction_pct, 1)
        constrained_gdf.loc[0, 'total_area_reduction_pct'] = round(total_reduction, 1)
        
        # Update polygon count
        constrained_gdf.loc[0, 'polygon_count'] = constrained_count
        
        print(f"    ✅ Hydrography constraints applied successfully")
        return constrained_gdf
        
    except Exception as e:
        print(f"    ❌ Error applying hydrography constraints: {e}")
        return distribution_gdf

def export_intermediate_files(plots_df, buffer_gdf, merged_gdf, species_name, buffer_size_desc):
    """Export intermediate files for QGIS editing"""
    print(f"  💾 Exporting intermediate files for QGIS editing...")
    
    species_filename = species_name.lower().replace(" ", "_")
    buffer_desc = buffer_size_desc.replace(" ", "_").replace("(", "").replace(")", "").lower()
    
    # Export 1: Plot points with attributes
    print(f"    📍 Exporting plot points...")
    
    # Calculate percentile rank for total carbon
    plots_df_enhanced = plots_df.copy()
    plots_df_enhanced['CARBON_PERCENTILE'] = plots_df_enhanced['TOTAL_CARBON_AG'].rank(pct=True) * 100
    plots_df_enhanced['CARBON_PERCENTILE'] = plots_df_enhanced['CARBON_PERCENTILE'].round(1)
    
    plot_points_gdf = gpd.GeoDataFrame(
        plots_df_enhanced, 
        geometry=[Point(lon, lat) for lon, lat in zip(plots_df_enhanced['LON'], plots_df_enhanced['LAT'])],
        crs='EPSG:4326'
    )
    
    plot_points_file = OUTPUT_DIR / f"{species_filename}_plots_{buffer_desc}_points.geojson"
    plot_points_gdf.to_file(plot_points_file)
    
    # Export 2: Individual buffer circles (for size adjustment)
    print(f"    🔵 Exporting individual buffer circles...")
    individual_buffers_file = OUTPUT_DIR / f"{species_filename}_individual_buffers_{buffer_desc}.geojson"
    
    # Add buffer metadata to each circle
    buffer_export_gdf = buffer_gdf.copy()
    buffer_export_gdf['buffer_km'] = BUFFER_RADIUS_KM
    buffer_export_gdf['species'] = species_name
    buffer_export_gdf['edit_instructions'] = 'Resize circles as needed, then dissolve overlaps'
    
    buffer_export_gdf.to_file(individual_buffers_file)
    
    # Export 3: Merged habitat areas (for gap filling/anomaly removal) 
    print(f"    🌲 Exporting merged habitat areas...")
    merged_habitat_file = OUTPUT_DIR / f"{species_filename}_merged_habitat_{buffer_desc}_EDIT_ME.geojson"
    
    # Add editing instructions
    merged_export_gdf = merged_gdf.copy()
    merged_export_gdf['edit_instructions'] = 'Fill gaps, remove anomalies, refine boundaries'
    merged_export_gdf['next_step'] = 'Save as _edited.geojson when done'
    
    merged_export_gdf.to_file(merged_habitat_file)
    
    # Export 4: QGIS project setup instructions
    instructions_file = OUTPUT_DIR / f"{species_filename}_{buffer_desc}_QGIS_EDITING_INSTRUCTIONS.md"
    
    instructions = f"""# {species_name} Habitat Map Editing Instructions

## Files Created:
1. **{plot_points_file.name}** - Original FIA plot locations (reference)
2. **{individual_buffers_file.name}** - Individual buffer circles (adjust sizes)
3. **{merged_habitat_file.name}** - Merged habitat areas (edit boundaries)

## QGIS Editing Workflow:

### Step 1: Load Files
1. Open QGIS
2. Load all 3 files above
3. Set CRS to EPSG:4326 (WGS84)

### Step 2: Review Current Buffers
1. Style plot points as small red dots
2. Style individual buffers as transparent circles with colored borders
3. Style merged habitat as semi-transparent green fill
4. Add basemap (OpenStreetMap or satellite imagery)

### Step 3: Edit Options

**Option A - Adjust Individual Buffers:**
1. Edit individual_buffers layer
2. Select circles and resize using scale tool
3. Move circles if needed to better align with habitat
4. When satisfied, use Vector > Geoprocessing > Dissolve to merge overlaps
5. Save result as: `{species_filename}_custom_buffers_{buffer_desc}_edited.geojson`

**Option B - Edit Merged Habitat:**
1. Edit merged_habitat layer  
2. Use Add Feature tool to fill gaps between suitable areas
3. Use Delete Feature or Edit tool to remove anomalous areas
4. Smooth boundaries using Simplify tool if needed
5. Save as: `{species_filename}_merged_habitat_{buffer_desc}_edited.geojson`

### Step 4: Final Output
- Save your edited version with `_edited` suffix
- Run the script again with `--use-edited-file` to apply environmental constraints

## Buffer Size Reference:
- Current buffer: {BUFFER_RADIUS_KM}km radius ({BUFFER_RADIUS_KM*2}km diameter)
- Buffer area per plot: ~{3.14159 * BUFFER_RADIUS_KM**2:.1f} sq km
- Use this as reference when manually adjusting sizes

## Tips:
- Zoom in close when editing boundaries
- Use satellite imagery to identify actual forest vs non-forest
- Consider topography - avoid placing habitat across ridges/valleys inappropriately
- Fill small gaps (<1km) between similar habitat areas
- Remove isolated small patches that seem unrealistic
"""
    
    with open(instructions_file, 'w') as f:
        f.write(instructions)
    
    print(f"  ✅ Exported intermediate files:")
    print(f"    📍 Plot points: {plot_points_file.name}")
    print(f"    🔵 Individual buffers: {individual_buffers_file.name}")  
    print(f"    🌲 Merged habitat: {merged_habitat_file.name}")
    print(f"    📝 Instructions: {instructions_file.name}")
    
    return {
        'plot_points': plot_points_file,
        'individual_buffers': individual_buffers_file,
        'merged_habitat': merged_habitat_file,
        'instructions': instructions_file
    }

def load_edited_habitat_file(species_name, buffer_desc):
    """Load a manually edited habitat file from QGIS"""
    species_filename = species_name.lower().replace(" ", "_")
    
    # Look for edited files
    possible_edited_files = [
        OUTPUT_DIR / f"{species_filename}_merged_habitat_{buffer_desc}_edited.geojson",
        OUTPUT_DIR / f"{species_filename}_custom_buffers_{buffer_desc}_edited.geojson",
        OUTPUT_DIR / f"{species_filename}_habitat_{buffer_desc}_edited.geojson"
    ]
    
    for edited_file in possible_edited_files:
        if edited_file.exists():
            print(f"  📂 Loading edited habitat file: {edited_file.name}")
            try:
                edited_gdf = gpd.read_file(edited_file)
                print(f"    ✅ Loaded {len(edited_gdf)} edited habitat features")
                return edited_gdf
            except Exception as e:
                print(f"    ❌ Error loading {edited_file}: {e}")
                continue
    
    return None

def main():
    """Main execution function"""
    # Parse arguments
    parser = argparse.ArgumentParser(description='Cached Species Distribution Mapping')
    parser.add_argument('--no-cache', action='store_true', help='Skip using cached masks')
    parser.add_argument('--clear-cache', action='store_true', help='Clear all cached masks')
    parser.add_argument('--species', choices=['douglas-fir', 'ponderosa-pine', 'black-cottonwood'], default='ponderosa-pine', 
                       help='Species to process (default: ponderosa-pine)')
    parser.add_argument('--buffer-size', choices=list(BUFFER_OPTIONS.keys()), default=DEFAULT_BUFFER_SIZE,
                       help=f'Buffer size preset (default: {DEFAULT_BUFFER_SIZE})')
    parser.add_argument('--buffer-km', type=float, help='Custom buffer radius in km (overrides --buffer-size)')
    parser.add_argument('--workflow', choices=['buffers_only', 'staged', 'full'], default='full',
                       help='Workflow mode: buffers_only=export for editing, staged=check for edited files, full=complete process')
    parser.add_argument('--use-edited-file', action='store_true', help='Look for and use manually edited habitat file')
    parser.add_argument('--export-all-plots', action='store_true', help='Export ALL plots (with/without target species) for comparison analysis')
    parser.add_argument('--pause-on-high-loss', action='store_true', help='Pause processing when area loss exceeds 70% to review constraints')
    parser.add_argument('--tcc-threshold', type=int, help='Custom tree canopy cover threshold percentage (overrides default)', metavar='PCT')
    parser.add_argument('--min-artifact-size', type=float, default=0.1, help='Minimum polygon size in sq km (smaller polygons removed as artifacts)', metavar='SQKM')
    parser.add_argument('--use-species-tcc-mask', action='store_true', help='Use species-specific TCC mask creation instead of regional masks')
    parser.add_argument('--regional-tcc-mask', type=str, help='Path to specific regional TCC mask file to use')
    parser.add_argument('--auto-create-tcc-mask', action='store_true', help='Auto-create missing regional TCC masks when needed')
    parser.add_argument('--find-best-tcc-mask', action='store_true', help='Find and use best available regional TCC mask')
    parser.add_argument('--list-tcc-masks', action='store_true', help='List available regional TCC masks and exit')
    
    args = parser.parse_args()
    
    use_cache = not args.no_cache
    
    # List TCC masks and exit if requested
    if args.list_tcc_masks:
        print("📋 AVAILABLE REGIONAL TCC MASKS")
        print("="*50)
        
        mask_files = list(CACHE_DIR.glob("pnw_tcc_mask_*.geojson"))
        if mask_files:
            for mask_file in sorted(mask_files):
                file_size_mb = mask_file.stat().st_size / (1024 * 1024)
                print(f"  📁 {mask_file.name} ({file_size_mb:.1f} MB)")
                
                # Try to parse parameters from filename
                try:
                    name_parts = mask_file.stem.split('_')
                    if len(name_parts) >= 6:
                        threshold = name_parts[3].replace('pct', '')
                        resolution = name_parts[4].replace('m', '')
                        min_area = name_parts[5].replace('sqkm', '')
                        print(f"     Threshold: {threshold}%, Resolution: {resolution}m, Min area: {min_area} sq km")
                except:
                    pass
        else:
            print("  (No regional TCC mask files found)")
            print(f"\n💡 Create some with:")
            print(f"  python scripts/create_regional_tcc_mask.py --create-batch")
        
        print(f"\n🔍 Find best match for specific parameters:")
        print(f"  python scripts/create_regional_tcc_mask.py --find-best-match 5 30 2.0")
        return
    
    # Set custom tree canopy threshold if provided
    global MIN_TREE_COVER_PCT, USE_REGIONAL_TCC_MASK, REGIONAL_TCC_MASK_FILE
    if args.tcc_threshold:
        MIN_TREE_COVER_PCT = args.tcc_threshold
        print(f"🌳 Using custom tree canopy threshold: {MIN_TREE_COVER_PCT}%")
    
    # Set TCC mask mode
    if args.use_species_tcc_mask:
        USE_REGIONAL_TCC_MASK = False
        print(f"🔧 Using species-specific TCC mask creation")
    
    if args.regional_tcc_mask:
        REGIONAL_TCC_MASK_FILE = args.regional_tcc_mask
        print(f"🔧 Using specified regional TCC mask: {args.regional_tcc_mask}")
    
    if args.auto_create_tcc_mask:
        print(f"🔧 Auto-creation of missing TCC masks enabled")
    
    if args.find_best_tcc_mask:
        print(f"🔧 Best-match TCC mask selection enabled")
    
    # Set buffer size
    global BUFFER_RADIUS_KM, BUFFER_RADIUS_DEG
    if args.buffer_km:
        BUFFER_RADIUS_KM = args.buffer_km
        buffer_desc = f"custom {args.buffer_km}km radius"
    else:
        BUFFER_RADIUS_KM = BUFFER_OPTIONS[args.buffer_size]['km']
        buffer_desc = BUFFER_OPTIONS[args.buffer_size]['desc']
    
    BUFFER_RADIUS_DEG = BUFFER_RADIUS_KM / 111.0
    
    # Set species configuration
    # Check for custom species via environment variables first
    import os
    if 'CUSTOM_SPECIES_CODE' in os.environ and 'CUSTOM_SPECIES_NAME' in os.environ:
        SPECIES_CODE = int(os.environ['CUSTOM_SPECIES_CODE'])
        SPECIES_NAME = os.environ['CUSTOM_SPECIES_NAME']
        MAX_ELEVATION_FT = 6000  # Default elevation limit for custom species
        print(f"🔧 Using custom species from environment: {SPECIES_NAME} (Code: {SPECIES_CODE})")
    elif args.species == 'douglas-fir':
        SPECIES_CODE = DOUGLAS_FIR_CODE
        SPECIES_NAME = DOUGLAS_FIR_NAME
        MAX_ELEVATION_FT = DOUGLAS_FIR_MAX_ELEVATION_FT
    elif args.species == 'black-cottonwood':
        SPECIES_CODE = BLACK_COTTONWOOD_CODE
        SPECIES_NAME = BLACK_COTTONWOOD_NAME
        MAX_ELEVATION_FT = 4000  # Black Cottonwood typically found at lower elevations
    else:  # ponderosa-pine
        SPECIES_CODE = PONDEROSA_PINE_CODE
        SPECIES_NAME = PONDEROSA_PINE_NAME
        MAX_ELEVATION_FT = PONDEROSA_PINE_MAX_ELEVATION_FT
    
    # Clear cache if requested
    if args.clear_cache:
        print("🗑️  Clearing cache...")
        cache_files = list(CACHE_DIR.glob("*.geojson"))
        for cache_file in cache_files:
            cache_file.unlink()
            print(f"    Deleted: {cache_file.name}")
        print(f"✅ Cleared {len(cache_files)} cache files")
        return
    
    # Export all plots if requested
    if args.export_all_plots:
        target_codes = [SPECIES_CODE] if args.species in ['douglas-fir', 'ponderosa-pine', 'black-cottonwood'] else [DOUGLAS_FIR_CODE, PONDEROSA_PINE_CODE]
        export_all_plots_comparison(target_codes)
        return
    
    print(f"🌲 CREATING PACIFIC NORTHWEST {SPECIES_NAME.upper()} HABITAT MAP")
    print("="*77)
    print(f"Species: {SPECIES_NAME}")
    print(f"Buffer: {buffer_desc}")
    print(f"Workflow: {args.workflow}")
    print(f"Coverage: Washington, Oregon, Idaho (Pacific Northwest)")
    print(f"Bounds: {PACIFIC_NW_BOUNDS} (clipped at US-Canada border)")
    print(f"Caching: {'Enabled' if use_cache else 'Disabled'}")
    
    if args.workflow == 'buffers_only':
        print(f"🎯 Mode: Export buffers for QGIS editing (no environmental constraints)")
    elif args.workflow == 'staged':
        print(f"🎯 Mode: Check for edited files, then apply environmental constraints")
    elif args.use_edited_file:
        print(f"🎯 Mode: Use manually edited habitat file + apply environmental constraints") 
    else:
        print(f"🎯 Mode: Full automated processing")
    
    # Show constraints that will be applied (unless buffers_only mode)
    if args.workflow != 'buffers_only':
        constraints = []
        if APPLY_ELEVATION_MASK:
            constraints.append(f"Elevation: exclude >{MAX_ELEVATION_FT:,} feet")
        if APPLY_TCC_MASK:
            constraints.append(f"Tree cover: exclude <{MIN_TREE_COVER_PCT}% canopy")
        if APPLY_HYDROGRAPHY_MASK:
            constraints.append(f"Hydrography: exclude water bodies")
        
        if constraints:
            print(f"Environmental Constraints: {' | '.join(constraints)}")
        else:
            print(f"Environmental Constraints: None")
    
    print()
    
    try:
        # Step 1: Calculate species-specific elevation limits from data
        min_elevation_ft, max_elevation_ft = calculate_species_elevation_limits(SPECIES_CODE, SPECIES_NAME)
        
        # Override the default MAX_ELEVATION_FT if we have data-driven limits
        if max_elevation_ft is not None:
            MAX_ELEVATION_FT = max_elevation_ft
            print(f"  ✅ Using data-driven elevation limits: {min_elevation_ft or 0:,} - {max_elevation_ft:,} ft")
        else:
            min_elevation_ft = None  # Use default one-sided constraint
            print(f"  ⚠️  Using default elevation limit: >{MAX_ELEVATION_FT:,} ft")
        
        # Step 2: Get species plot data
        plots_df = get_species_plots(SPECIES_CODE, SPECIES_NAME)
        
        if len(plots_df) == 0:
            print(f"❌ No {SPECIES_NAME} plots found in database")
            return
        
        # Step 2: Create buffers around plots
        buffer_gdf = create_plot_buffers(plots_df)
        
        # Step 3: Check workflow mode
        if args.workflow == 'buffers_only':
            # Export files for editing and exit
            print(f"📤 EXPORTING FILES FOR QGIS EDITING")
            print("="*50)
            
            # Create a quick merge for the export
            merged_gdf = merge_overlapping_buffers(buffer_gdf, SPECIES_CODE, SPECIES_NAME)
            
            # Export intermediate files
            export_files = export_intermediate_files(plots_df, buffer_gdf, merged_gdf, SPECIES_NAME, buffer_desc)
            
            print(f"\n✅ EXPORT COMPLETE!")
            print(f"📁 Files exported to: {OUTPUT_DIR}")
            print(f"📝 Next steps:")
            print(f"   1. Open QGIS and load the exported files")
            print(f"   2. Edit habitat boundaries as needed")
            print(f"   3. Save edited version with '_edited' suffix")
            print(f"   4. Run script with --workflow staged to apply environmental constraints")
            print(f"\n📖 See detailed instructions in: {export_files['instructions'].name}")
            return
        
        # Step 4: Check for edited habitat file (staged or use-edited-file modes)
        edited_habitat_gdf = None
        if args.workflow == 'staged' or args.use_edited_file:
            buffer_desc_clean = buffer_desc.replace(" ", "_").replace("(", "").replace(")", "").lower()
            edited_habitat_gdf = load_edited_habitat_file(SPECIES_NAME, buffer_desc_clean)
            
            if edited_habitat_gdf is not None:
                print(f"✅ Using manually edited habitat file")
                distribution_gdf = edited_habitat_gdf
                
                # Ensure required columns exist with proper values
                if 'species' not in distribution_gdf.columns or distribution_gdf['species'].iloc[0] is None:
                    distribution_gdf['species'] = SPECIES_NAME
                if 'species_code' not in distribution_gdf.columns or distribution_gdf['species_code'].iloc[0] is None:
                    distribution_gdf['species_code'] = SPECIES_CODE
                if 'total_plots' not in distribution_gdf.columns or distribution_gdf['total_plots'].iloc[0] is None:
                    distribution_gdf['total_plots'] = len(plots_df)
                if 'buffer_radius_km' not in distribution_gdf.columns or distribution_gdf['buffer_radius_km'].iloc[0] is None:
                    distribution_gdf['buffer_radius_km'] = BUFFER_RADIUS_KM
                    
            elif args.workflow == 'staged':
                print(f"⚠️  No edited habitat file found - creating fresh buffers for editing")
                
                # Create merged buffers and export for editing
                merged_gdf = merge_overlapping_buffers(buffer_gdf, SPECIES_CODE, SPECIES_NAME)
                export_files = export_intermediate_files(plots_df, buffer_gdf, merged_gdf, SPECIES_NAME, buffer_desc)
                
                print(f"\n📤 FILES EXPORTED FOR EDITING")
                print(f"📝 Edit the habitat boundaries in QGIS, then run again with --workflow staged")
                print(f"📖 See instructions in: {export_files['instructions'].name}")
                return
        
        # Step 5: Create merged habitat areas (if not using edited file)
        if edited_habitat_gdf is None:
            distribution_gdf = merge_overlapping_buffers(buffer_gdf, SPECIES_CODE, SPECIES_NAME)
            
            # Export intermediate files if requested
            if EXPORT_INTERMEDIATE_FILES and args.workflow == 'full':
                export_intermediate_files(plots_df, buffer_gdf, distribution_gdf, SPECIES_NAME, buffer_desc)
        
        # Step 6: Apply environmental constraints (if not buffers_only mode)
        if args.workflow != 'buffers_only':
            if APPLY_ELEVATION_MASK or APPLY_TCC_MASK or APPLY_HYDROGRAPHY_MASK:
                bounds = distribution_gdf.total_bounds
                
                # Track cumulative area loss
                original_geometry = distribution_gdf.iloc[0].geometry
                original_area_km2 = calculate_area_km2(original_geometry)
                distribution_gdf.loc[0, 'original_area_km2'] = round(original_area_km2, 1)
                
                print(f"\n  📊 STARTING ENVIRONMENTAL CONSTRAINTS")
                print(f"    🏞️  Original edited habitat: {original_area_km2:,.0f} km²")
                
                steps_completed = []
                
                # Apply elevation constraints
                if APPLY_ELEVATION_MASK:
                    elevation_mask = create_elevation_mask(bounds, min_elevation_ft=min_elevation_ft, max_elevation_ft=MAX_ELEVATION_FT, use_cache=use_cache)
                    distribution_gdf = apply_elevation_constraints(distribution_gdf, elevation_mask, min_elevation_ft, MAX_ELEVATION_FT)
                    steps_completed.append("elevation")
                    
                    # Log cumulative loss
                    current_area_km2 = distribution_gdf.loc[0, 'elevation_area_km2']
                    log_cumulative_loss(original_area_km2, current_area_km2, steps_completed, args.pause_on_high_loss)
                
                # Apply tree cover constraints
                if APPLY_TCC_MASK:
                    if USE_REGIONAL_TCC_MASK:
                        # Use pre-created regional TCC mask (recommended for multiple species)
                        print(f"  🌳 Using regional TCC mask approach...")
                        
                        if args.find_best_tcc_mask:
                            # Find best available mask instead of exact match
                            print(f"    🔍 Finding best available regional TCC mask...")
                            best_mask_file = find_best_regional_tcc_mask(MIN_TREE_COVER_PCT, TCC_PROCESSING_RESOLUTION, MIN_LOW_COVER_AREA_SQKM)
                            if best_mask_file:
                                tree_cover_mask = load_regional_tcc_mask(threshold_pct=None, mask_file=str(best_mask_file))
                            else:
                                tree_cover_mask = None
                        elif args.auto_create_tcc_mask:
                            # Auto-create if missing
                            print(f"    🔨 Auto-creating TCC mask if needed...")
                            auto_created_file = auto_create_regional_tcc_mask(MIN_TREE_COVER_PCT, TCC_PROCESSING_RESOLUTION, MIN_LOW_COVER_AREA_SQKM)
                            if auto_created_file:
                                tree_cover_mask = load_regional_tcc_mask(threshold_pct=None, mask_file=str(auto_created_file))
                            else:
                                tree_cover_mask = None
                        else:
                            # Standard loading with intelligent fallback
                            tree_cover_mask = load_regional_tcc_mask(
                                threshold_pct=MIN_TREE_COVER_PCT, 
                                mask_file=REGIONAL_TCC_MASK_FILE
                            )
                        
                        if tree_cover_mask is None:
                            # Critical error - do not fall back to species-specific mask creation
                            error_msg = f"❌ CRITICAL ERROR: Regional TCC mask could not be loaded!"
                            print(f"  {error_msg}")
                            print(f"  📁 Attempted to load: {REGIONAL_TCC_MASK_FILE}")
                            print(f"  💡 Available TCC masks in cache/species_masks/:")
                            
                            # List available TCC masks
                            tcc_masks = list(Path("cache/species_masks").glob("pnw_tcc_mask_*.geojson"))
                            if tcc_masks:
                                for mask_file in sorted(tcc_masks):
                                    file_size_mb = mask_file.stat().st_size / (1024 * 1024)
                                    print(f"    📄 {mask_file.name} ({file_size_mb:.1f} MB)")
                            else:
                                print(f"    (No TCC mask files found)")
                            
                            print(f"  🔧 Fix by using one of these options:")
                            print(f"    --regional-tcc-mask cache/species_masks/[filename]")
                            print(f"    --find-best-tcc-mask")
                            print(f"    --auto-create-tcc-mask")
                            
                            raise ValueError(error_msg)
                    else:
                        # Create species-specific mask (legacy mode)
                        print(f"  🌳 Creating species-specific TCC mask...")
                        current_distribution_geom = distribution_gdf.iloc[0].geometry
                        tree_cover_mask = create_tree_cover_mask(bounds, distribution_geom=current_distribution_geom, min_tree_cover_pct=MIN_TREE_COVER_PCT, use_cache=use_cache)
                    
                    distribution_gdf = apply_tree_cover_constraints(distribution_gdf, tree_cover_mask)
                    steps_completed.append("tree_cover")
                    
                    # Log cumulative loss
                    if 'tcc_area_km2' in distribution_gdf.columns:
                        current_area_km2 = distribution_gdf.loc[0, 'tcc_area_km2']
                        log_cumulative_loss(original_area_km2, current_area_km2, steps_completed, args.pause_on_high_loss)
                    else:
                        print(f"    ⚠️  Tree cover masking failed - skipping cumulative loss calculation")
                
                # Apply hydrography constraints
                if APPLY_HYDROGRAPHY_MASK:
                    hydrography_mask = create_hydrography_mask(bounds, use_cache=use_cache)
                    distribution_gdf = apply_hydrography_constraints(distribution_gdf, hydrography_mask)
                    steps_completed.append("hydrography")
                    
                    # Log cumulative loss
                    current_area_km2 = calculate_area_km2(distribution_gdf.iloc[0].geometry)
                    log_cumulative_loss(original_area_km2, current_area_km2, steps_completed, args.pause_on_high_loss)
        
        # Step 7: Apply final processing
        distribution_gdf = clip_to_canada_border(distribution_gdf)
        distribution_gdf = simplify_geometry(distribution_gdf)
        distribution_gdf = remove_small_artifacts(distribution_gdf, min_area_sqkm=args.min_artifact_size)
        
        # Step 8: Save output
        suffix_parts = []
        if APPLY_ELEVATION_MASK and args.workflow != 'buffers_only':
            suffix_parts.append("elevation")
        if APPLY_TCC_MASK and args.workflow != 'buffers_only':
            suffix_parts.append("tcc")
        if APPLY_HYDROGRAPHY_MASK and args.workflow != 'buffers_only':
            suffix_parts.append("hydro")
        
        output_suffix = "_" + "_".join(suffix_parts) + "_constrained" if suffix_parts else ""
        cache_suffix = "_cached" if use_cache else ""
        buffer_suffix = f"_{args.buffer_size}buf" if args.buffer_size != DEFAULT_BUFFER_SIZE else ""
        edited_suffix = "_from_edited" if edited_habitat_gdf is not None else ""
        
        species_filename = SPECIES_NAME.lower().replace(" ", "_")
        output_file = OUTPUT_DIR / f"{species_filename}_continuous_distribution_pnw{buffer_suffix}{output_suffix}{edited_suffix}{cache_suffix}.geojson"
        
        # Save optimized version for MBTiles (smaller file size, fewer properties)
        original_size_mb = 0
        try:
            save_optimized_geojson(distribution_gdf, output_file, precision=6)
            
            # Calculate size reduction if we can
            optimized_size_mb = output_file.stat().st_size / (1024 * 1024)
            print(f"  📁 File size: {optimized_size_mb:.1f} MB (optimized for MBTiles)")
            
        except Exception as e:
            print(f"  ⚠️  Optimization failed, using standard GeoJSON: {e}")
            distribution_gdf.to_file(output_file, driver='GeoJSON')
            optimized_size_mb = output_file.stat().st_size / (1024 * 1024)
            print(f"  📁 File size: {optimized_size_mb:.1f} MB")
        
        # Success output
        print(f"\n✅ SUCCESS!")
        print(f"📁 Output file: {output_file}")
        print(f"📊 Summary:")
        row = distribution_gdf.iloc[0]
        print(f"   • Species: {row.get('species', SPECIES_NAME) or SPECIES_NAME}")
        
        # Handle None values in formatting
        total_plots = row.get('total_plots', len(plots_df) if 'plots_df' in locals() else 'N/A')
        if isinstance(total_plots, (int, float)):
            print(f"   • Source plots: {total_plots:,}")
        else:
            print(f"   • Source plots: {total_plots}")
            
        buffer_radius = row.get('buffer_radius_km', BUFFER_RADIUS_KM)
        print(f"   • Buffer radius: {buffer_radius}km")
        
        polygon_count = row.get('polygon_count', 'N/A')
        if isinstance(polygon_count, (int, float)):
            print(f"   • Polygon areas: {polygon_count:,}")
        else:
            print(f"   • Polygon areas: {polygon_count}")
        
        if edited_habitat_gdf is not None:
            print(f"   • Source: Manually edited habitat boundaries")
        
        # Display constraint information
        if args.workflow != 'buffers_only':
            if APPLY_ELEVATION_MASK and 'elevation_constrained' in row:
                min_elev = row.get('min_elevation_ft', None)
                max_elev = row.get('max_elevation_ft', 'N/A')
                if min_elev is not None:
                    print(f"   • Elevation range: {min_elev:g} - {max_elev:g} feet")
                else:
                    print(f"   • Elevation limit: >{max_elev:g} feet")
                
            if APPLY_TCC_MASK and 'tree_cover_constrained' in row:
                print(f"   • Tree cover limit: ≥{row.get('min_tree_cover_pct', 'N/A')}% canopy")
                
            if APPLY_HYDROGRAPHY_MASK and 'hydrography_constrained' in row:
                print(f"   • Water body exclusions applied")
                
            # Show detailed area breakdown
            original_area_km2 = row.get('original_area_km2', 'N/A')
            final_area_km2 = calculate_area_km2(row.geometry) if row.geometry is not None else 'N/A'
            
            if isinstance(original_area_km2, (int, float)) and isinstance(final_area_km2, (int, float)):
                total_reduction_pct = ((original_area_km2 - final_area_km2) / original_area_km2 * 100)
                print(f"   • Area change: {original_area_km2:,.0f} km² → {final_area_km2:,.0f} km² ({total_reduction_pct:.1f}% reduction)")
            else:
                print(f"   • Total area reduction: {row.get('total_area_reduction_pct', 'N/A')}%")
        
        # Calculate file size
        file_size_mb = output_file.stat().st_size / (1024 * 1024)
        print(f"   • File size: {file_size_mb:.1f} MB")
        
        print(f"\n🗺️  Next steps:")
        if args.workflow == 'buffers_only':
            print(f"   1. Edit the exported files in QGIS")
            print(f"   2. Run with --workflow staged to apply environmental constraints")
        else:
            print(f"   1. Load {output_file.name} in QGIS or MapBox")
            print(f"   2. Style as green/forest colored polygon")
            print(f"   3. Use as foundation layer for additional analysis")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    import sys
    
    # Check for command line arguments
    if len(sys.argv) > 1:
        arg = sys.argv[1].lower()
        if arg in ['--no-constraints', '--original', '--basic']:
            print("🔧 All constraints disabled via command line argument")
            APPLY_ELEVATION_MASK = False
            APPLY_TCC_MASK = False
            APPLY_HYDROGRAPHY_MASK = False
        elif arg in ['--elevation-only']:
            print("🔧 Only elevation masking enabled via command line argument")
            APPLY_ELEVATION_MASK = True
            APPLY_TCC_MASK = False
            APPLY_HYDROGRAPHY_MASK = False
        elif arg in ['--tcc-only', '--tree-cover-only']:
            print("🔧 Only tree cover masking enabled via command line argument")
            APPLY_ELEVATION_MASK = False
            APPLY_TCC_MASK = True
            APPLY_HYDROGRAPHY_MASK = False
        elif arg in ['--hydro-only', '--hydrography-only']:
            print("🔧 Only hydrography masking enabled via command line argument")
            APPLY_ELEVATION_MASK = False
            APPLY_TCC_MASK = False
            APPLY_HYDROGRAPHY_MASK = True
        elif arg in ['--elevation-hydro', '--elevation-hydrography']:
            print("🔧 Elevation and hydrography masking enabled via command line argument")
            APPLY_ELEVATION_MASK = True
            APPLY_TCC_MASK = False
            APPLY_HYDROGRAPHY_MASK = True
        elif arg in ['--all-constraints', '--full']:
            print("🔧 All constraints enabled via command line argument")
            APPLY_ELEVATION_MASK = True
            APPLY_TCC_MASK = True
            APPLY_HYDROGRAPHY_MASK = True
        elif arg in ['--help', '-h']:
            print("SPECIES DISTRIBUTION MAPPING - Usage Examples")
            print("="*60)
            print()
            print("🎯 BASIC USAGE:")
            print("  python create_species_distribution_maps_cached.py")
            print("  # Creates full habitat map with default settings")
            print()
            print("🔧 BUFFER SIZE OPTIONS:")
            print("  --buffer-size small    # 1km radius (local habitat)")
            print("  --buffer-size medium   # 2km radius (regional habitat) [DEFAULT]")
            print("  --buffer-size large    # 4km radius (landscape connectivity)")
            print("  --buffer-size xlarge   # 6km radius (broad corridors)")
            print("  --buffer-km 3.5        # Custom radius in km")
            print()
            print("🔄 WORKFLOW OPTIONS:")
            print("  --workflow buffers_only  # Export files for QGIS editing")
            print("  --workflow staged        # Check for edited files, then process")
            print("  --workflow full          # Complete automated processing [DEFAULT]")
            print()
            print("📊 SPECIES OPTIONS:")
            print("  --species douglas-fir     # Douglas Fir habitat")
            print("  --species ponderosa-pine  # Ponderosa Pine habitat [DEFAULT]")
            print()
            print("🗺️  RECOMMENDED WORKFLOW FOR CUSTOM EDITING:")
            print("  1. python create_species_distribution_maps_cached.py --workflow buffers_only --buffer-size large")
            print("  2. Open QGIS, load exported files, edit habitat boundaries")
            print("  3. Save edited version with '_edited' suffix")
            print("  4. python create_species_distribution_maps_cached.py --workflow staged --buffer-size large")
            print()
            print("🔧 CONSTRAINT OPTIONS:")
            print("  --no-constraints     # No environmental constraints")
            print("  --elevation-only     # Only elevation masking")
            print("  --tcc-only          # Only tree cover masking")
            print("  --hydro-only        # Only hydrography masking")
            print("  --all-constraints   # All environmental constraints [DEFAULT]")
            print()
            print("💾 CACHE OPTIONS:")
            print("  --no-cache          # Skip using cached environmental masks")
            print("  --clear-cache       # Clear all cached masks")
            print()
            print("EXAMPLE COMMANDS:")
            print("-" * 40)
            print("# Export large buffers for manual editing:")
            print("python create_species_distribution_maps_cached.py --workflow buffers_only --buffer-size large")
            print()
            print("# Process edited habitat with all constraints:")
            print("python create_species_distribution_maps_cached.py --workflow staged --buffer-size large")
            print()
            print("# Quick Douglas Fir map with small buffers:")
            print("python create_species_distribution_maps_cached.py --species douglas-fir --buffer-size small")
            print()
            print("# Custom buffer size with no constraints:")
            print("python create_species_distribution_maps_cached.py --buffer-km 5.0 --no-constraints")
            print()
            print("# Export ALL plots for survey coverage analysis:")
            print("python create_species_distribution_maps_cached.py --export-all-plots --species ponderosa-pine")
            sys.exit(0)
    
    main() 