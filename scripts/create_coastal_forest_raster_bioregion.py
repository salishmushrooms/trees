#!/usr/bin/env python3
"""
Create Coastal Forest Bioregion using Raster-Based Approach

This script creates a coastal forest bioregion polygon using a raster-based approach:
1. Applies longitude constraint (west of -123.0°) for coastal definition
2. Applies elevation constraint (<500ft for shore pine habitat)
3. Applies tree cover threshold (>10% for forest definition)
4. Creates natural boundaries following topography and forest coverage
5. Validates against shore pine (Pinus contorta) plot locations

Data Assets Used:
- Elevation: pnw_elevation_120m_mapbox.tif (120m resolution)
- Tree Cover: pnw_tree_cover_30m_full.tif (30m resolution)
- Validation: Shore pine plots from FIA data
"""

import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.geometry import Point, Polygon, MultiPolygon
from shapely.ops import unary_union
import rasterio
from rasterio.mask import mask
from rasterio.features import shapes
from rasterio.warp import calculate_default_transform, reproject, Resampling
import warnings
from pathlib import Path
from tqdm import tqdm
import json

warnings.filterwarnings('ignore')

# Configuration - Coastal Forest Bioregion Parameters (Updated)
LONGITUDE_CUTOFF = -123.5     # More restrictive - closer to coast
MAX_ELEVATION_FT = 300        # Lower elevation for true coastal forests
MIN_TREE_COVER_PCT = 15       # Higher threshold for forest definition
OUTPUT_RESOLUTION = 120       # Match elevation data resolution
DISTANCE_FROM_COAST_KM = 30   # Maximum distance from coastline

# Geographic boundaries
NORTH_BOUNDARY = 49.0         # Northern extent (southern BC border)
SOUTH_BOUNDARY = 41.5         # Southern extent (northern CA)
WEST_BOUNDARY = -125.0        # Pacific Ocean boundary
EAST_BOUNDARY = -121.0        # Inland boundary (generous buffer)

# Exclusion zones (to remove Puget Sound and inland areas)
PUGET_SOUND_EXCLUSION = {
    'north': 48.5,
    'south': 47.0, 
    'west': -123.2,
    'east': -122.0
}

# Species validation
SHORE_PINE_SPCD = 108         # Lodgepole pine (includes shore pine)

# File paths
ELEVATION_RASTER = "outputs/mapbox_masks/pnw_elevation_120m_mapbox.tif"
TREE_COVER_RASTER = "outputs/mapbox_masks/pnw_tree_cover_30m_full.tif"
PLOT_DATA = "outputs/plot_carbon_percentiles_latest_surveys.geojson"
OUTPUT_DIR = Path("outputs/bioregions")

def convert_elevation_to_feet(elevation_meters):
    """Convert elevation from meters to feet"""
    return elevation_meters * 3.28084

def create_coastal_mask_raster():
    """Create coastal forest mask by combining elevation, tree cover, and longitude constraints"""
    print("=== Creating Coastal Forest Raster Mask ===")
    
    # Step 1: Load and process elevation data
    print(f"Loading elevation data from {ELEVATION_RASTER}")
    with rasterio.open(ELEVATION_RASTER) as elev_src:
        elev_data = elev_src.read(1)
        elev_transform = elev_src.transform
        elev_crs = elev_src.crs
        elev_profile = elev_src.profile
        
        # Convert to feet and apply elevation constraint
        elev_data_ft = convert_elevation_to_feet(elev_data)
        elevation_mask = (elev_data_ft <= MAX_ELEVATION_FT) & (elev_data_ft >= 0)
        
        print(f"Elevation constraint: {np.sum(elevation_mask)} pixels ≤ {MAX_ELEVATION_FT}ft")
    
    # Step 2: Use 120m tree cover data to match elevation resolution
    print(f"Loading tree cover data from 120m resolution file")
    tcc_120m_path = "outputs/mapbox_masks/pnw_tree_cover_120m_quarter.tif"
    
    with rasterio.open(tcc_120m_path) as tcc_src:
        # Check if resolutions match
        if tcc_src.shape == elev_data.shape:
            print("Using matching resolution data")
            tcc_data = tcc_src.read(1)
        else:
            print("Resampling tree cover to match elevation...")
            tcc_data = np.zeros_like(elev_data, dtype=np.float32)
            reproject(
                source=rasterio.band(tcc_src, 1),
                destination=tcc_data,
                src_transform=tcc_src.transform,
                src_crs=tcc_src.crs,
                dst_transform=elev_transform,
                dst_crs=elev_crs,
                resampling=Resampling.average
            )
        
        # Apply tree cover constraint
        tree_cover_mask = tcc_data >= MIN_TREE_COVER_PCT
        print(f"Tree cover constraint: {np.sum(tree_cover_mask)} pixels ≥ {MIN_TREE_COVER_PCT}%")
    
    # Step 3: Create geographic constraints
    print(f"Creating geographic constraints...")
    rows, cols = elev_data.shape
    
    # Create coordinate arrays
    row_indices, col_indices = np.meshgrid(np.arange(rows), np.arange(cols), indexing='ij')
    
    # Convert to coordinates
    lons, lats = rasterio.transform.xy(elev_transform, row_indices.flatten(), col_indices.flatten())
    lons = np.array(lons).reshape(rows, cols)
    lats = np.array(lats).reshape(rows, cols)
    
    # Apply longitude constraint (more restrictive)
    longitude_mask = lons <= LONGITUDE_CUTOFF
    print(f"Longitude constraint: {np.sum(longitude_mask)} pixels west of {LONGITUDE_CUTOFF}°")
    
    # Exclude Puget Sound area
    puget_sound_mask = ~(
        (lats >= PUGET_SOUND_EXCLUSION['south']) & 
        (lats <= PUGET_SOUND_EXCLUSION['north']) &
        (lons >= PUGET_SOUND_EXCLUSION['west']) & 
        (lons <= PUGET_SOUND_EXCLUSION['east'])
    )
    print(f"Puget Sound exclusion: {np.sum(~puget_sound_mask)} pixels excluded")
    
    # Combine geographic constraints
    geographic_mask = longitude_mask & puget_sound_mask
    
    # Step 4: Combine all constraints
    print("Combining all constraints...")
    coastal_forest_mask = elevation_mask & tree_cover_mask & geographic_mask
    
    print(f"Final coastal forest mask: {np.sum(coastal_forest_mask)} pixels")
    print(f"Approximate area: {np.sum(coastal_forest_mask) * (OUTPUT_RESOLUTION**2) / 1000000:.0f} km²")
    
    return coastal_forest_mask, elev_transform, elev_crs, elev_profile

def raster_to_polygons(mask_array, transform, crs):
    """Convert raster mask to vector polygons"""
    print("Converting raster mask to vector polygons...")
    
    # Convert boolean mask to uint8
    mask_uint8 = mask_array.astype(np.uint8)
    
    # Extract polygons
    from shapely.geometry import shape
    polygons = []
    for geom, value in shapes(mask_uint8, mask=mask_uint8, transform=transform):
        if value == 1:  # Only keep forest pixels
            # Convert dict to shapely geometry
            poly = shape(geom)
            polygons.append(poly)
    
    print(f"Extracted {len(polygons)} polygon features")
    
    if not polygons:
        print("Warning: No polygons extracted from mask")
        return gpd.GeoDataFrame()
    
    # Create GeoDataFrame
    gdf = gpd.GeoDataFrame(
        [{'coastal_forest': 1} for _ in polygons],
        geometry=polygons,
        crs=crs
    )
    
    return gdf

def simplify_and_merge_polygons(gdf, tolerance=0.005, min_area_sqkm=2.0):
    """Simplify and merge small polygons with more aggressive filtering"""
    print(f"Simplifying polygons (tolerance={tolerance}, min_area={min_area_sqkm}km²)...")
    
    if gdf.empty:
        return gdf
    
    # Convert to geographic CRS for area calculation
    gdf_geo = gdf.to_crs('EPSG:4326')
    
    # Calculate areas in km²
    gdf_geo['area_sqkm'] = gdf_geo.geometry.area * 111 * 111  # Approximate conversion
    
    # Filter by minimum area (more aggressive to remove offshore bands)
    large_polygons = gdf_geo[gdf_geo['area_sqkm'] >= min_area_sqkm]
    
    if large_polygons.empty:
        print("Warning: No polygons meet minimum area requirement")
        return gdf_geo
    
    print(f"Kept {len(large_polygons)} polygons ≥ {min_area_sqkm}km²")
    
    # Apply more aggressive simplification
    large_polygons['geometry'] = large_polygons.geometry.simplify(tolerance)
    
    # Additional filtering: remove very narrow polygons (likely offshore bands)
    filtered_polygons = []
    for idx, poly in large_polygons.iterrows():
        bounds = poly.geometry.bounds
        width = bounds[2] - bounds[0]  # longitude width  
        height = bounds[3] - bounds[1]  # latitude height
        
        # Skip very narrow coastal bands (aspect ratio > 20:1)
        aspect_ratio = max(width, height) / min(width, height) if min(width, height) > 0 else 0
        if aspect_ratio < 20:  # Keep polygons with reasonable aspect ratios
            filtered_polygons.append(poly)
    
    if not filtered_polygons:
        print("Warning: No polygons remain after aspect ratio filtering")
        large_polygons = gdf_geo[gdf_geo['area_sqkm'] >= min_area_sqkm]  # Fall back
    else:
        large_polygons = gpd.GeoDataFrame(filtered_polygons)
        print(f"After aspect ratio filtering: {len(large_polygons)} polygons")
    
    # Merge all polygons into single multipolygon
    try:
        merged_geometry = unary_union(large_polygons.geometry)
    except Exception as e:
        print(f"Warning: Error merging polygons: {e}")
        print("Using individual polygons instead of merged geometry")
        # Create MultiPolygon from individual polygons
        from shapely.geometry import MultiPolygon
        valid_polygons = []
        for geom in large_polygons.geometry:
            if geom.is_valid:
                valid_polygons.append(geom)
            else:
                # Try to fix invalid geometry
                try:
                    fixed_geom = geom.buffer(0)
                    if fixed_geom.is_valid:
                        valid_polygons.append(fixed_geom)
                except:
                    continue
        
        if valid_polygons:
            merged_geometry = MultiPolygon(valid_polygons) if len(valid_polygons) > 1 else valid_polygons[0]
        else:
            print("Error: No valid polygons to merge")
            return gpd.GeoDataFrame()
    
    # Create final bioregion
    coastal_bioregion = gpd.GeoDataFrame(
        [{
            'region_name': 'Coastal Forests',
            'region_code': 'coastal_forest',
            'description': 'Pacific coastal forests <300ft elevation with >15% tree cover, west of -123.5°, excluding Puget Sound',
            'elevation_max_ft': MAX_ELEVATION_FT,
            'longitude_cutoff': LONGITUDE_CUTOFF,
            'min_tree_cover_pct': MIN_TREE_COVER_PCT,
            'area_sqkm': merged_geometry.area * 111 * 111,
            'created_date': pd.Timestamp.now().isoformat(),
            'method': 'raster-based_constraints',
            'validation_species': 'Shore pine (Pinus contorta var. contorta)',
            'data_sources': 'FIA plots, NLCD tree cover, USGS elevation'
        }],
        geometry=[merged_geometry],
        crs='EPSG:4326'
    )
    
    return coastal_bioregion

def validate_with_shore_pine_plots(bioregion_gdf):
    """Validate bioregion against shore pine plot locations"""
    print("Validating bioregion against shore pine plots...")
    
    try:
        # Load plot data
        plots_gdf = gpd.read_file(PLOT_DATA)
        
        # Filter for shore pine plots (Pinus contorta)
        # Look for carbon columns indicating presence
        shore_pine_cols = [col for col in plots_gdf.columns if 'LODGEPOLE_PINE' in col and 'CARBON' in col]
        
        if not shore_pine_cols:
            print("Warning: No shore pine carbon columns found in plot data")
            return {}
        
        # Get plots with shore pine presence
        shore_pine_col = shore_pine_cols[0]  # Use first match
        shore_pine_plots = plots_gdf[plots_gdf[shore_pine_col] > 0]
        
        # Filter for coastal elevation range
        elev_col = 'ELEVATION_FT' if 'ELEVATION_FT' in shore_pine_plots.columns else 'ELEV'
        if elev_col in shore_pine_plots.columns:
            coastal_pine_plots = shore_pine_plots[shore_pine_plots[elev_col] <= MAX_ELEVATION_FT]
        else:
            coastal_pine_plots = shore_pine_plots
        
        # Filter for longitude constraint
        coastal_pine_plots = coastal_pine_plots[coastal_pine_plots.geometry.x <= LONGITUDE_CUTOFF]
        
        print(f"Found {len(coastal_pine_plots)} shore pine plots in coastal range")
        
        if len(coastal_pine_plots) == 0:
            return {'validation': 'no_plots_found'}
        
        # Check how many plots fall within bioregion
        plots_in_bioregion = gpd.sjoin(coastal_pine_plots, bioregion_gdf, how='inner', predicate='within')
        
        validation_stats = {
            'total_shore_pine_plots': len(shore_pine_plots),
            'coastal_shore_pine_plots': len(coastal_pine_plots),
            'plots_captured_by_bioregion': len(plots_in_bioregion),
            'capture_rate': len(plots_in_bioregion) / len(coastal_pine_plots) if len(coastal_pine_plots) > 0 else 0,
            'validation_success': len(plots_in_bioregion) >= len(coastal_pine_plots) * 0.8  # 80% capture rate
        }
        
        print(f"Validation results:")
        print(f"  Total shore pine plots: {validation_stats['total_shore_pine_plots']}")
        print(f"  Coastal shore pine plots: {validation_stats['coastal_shore_pine_plots']}")
        print(f"  Plots captured by bioregion: {validation_stats['plots_captured_by_bioregion']}")
        print(f"  Capture rate: {validation_stats['capture_rate']:.1%}")
        
        # Save validation plots for inspection
        if len(coastal_pine_plots) > 0:
            plot_output = OUTPUT_DIR / "shore_pine_validation_plots.geojson"
            coastal_pine_plots.to_file(plot_output, driver='GeoJSON')
            print(f"Saved validation plots to {plot_output}")
        
        return validation_stats
        
    except Exception as e:
        print(f"Warning: Could not validate with plot data: {e}")
        return {'validation': 'error', 'error': str(e)}

def main():
    """Main processing function"""
    print("=== Creating Coastal Forest Bioregion (Raster-Based) ===")
    
    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Create raster mask combining all constraints
    mask_array, transform, crs, profile = create_coastal_mask_raster()
    
    # Step 2: Convert raster to vector polygons
    coastal_polygons = raster_to_polygons(mask_array, transform, crs)
    
    # Step 3: Simplify and merge polygons
    coastal_bioregion = simplify_and_merge_polygons(coastal_polygons)
    
    if coastal_bioregion.empty:
        print("Error: No coastal forest bioregion created")
        return
    
    # Step 4: Validate against shore pine plots
    validation_stats = validate_with_shore_pine_plots(coastal_bioregion)
    
    # Add validation results to bioregion attributes
    if validation_stats:
        for key, value in validation_stats.items():
            coastal_bioregion[key] = value
    
    # Step 5: Save outputs
    output_file = OUTPUT_DIR / "coastal_forest_raster_bioregion_v2.geojson"
    print(f"\nSaving coastal forest bioregion to {output_file}")
    coastal_bioregion.to_file(output_file, driver='GeoJSON')
    
    # Save processing summary
    summary = {
        'method': 'raster_based_constraints_v2',
        'parameters': {
            'longitude_cutoff': LONGITUDE_CUTOFF,
            'max_elevation_ft': MAX_ELEVATION_FT,
            'min_tree_cover_pct': MIN_TREE_COVER_PCT,
            'output_resolution_m': OUTPUT_RESOLUTION,
            'distance_from_coast_km': DISTANCE_FROM_COAST_KM,
            'puget_sound_excluded': True,
            'min_area_sqkm': 2.0,
            'simplification_tolerance': 0.005
        },
        'data_sources': {
            'elevation': ELEVATION_RASTER,
            'tree_cover': TREE_COVER_RASTER,
            'validation': PLOT_DATA
        },
        'validation': validation_stats,
        'created': pd.Timestamp.now().isoformat()
    }
    
    summary_file = OUTPUT_DIR / "coastal_forest_raster_summary_v2.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    
    # Print final summary
    print("\n=== Final Summary ===")
    if not coastal_bioregion.empty:
        area_sqkm = coastal_bioregion.iloc[0]['area_sqkm']
        print(f"Coastal forest bioregion area: {area_sqkm:.0f} km²")
        
        if validation_stats and 'capture_rate' in validation_stats:
            print(f"Shore pine plot validation: {validation_stats['capture_rate']:.1%} capture rate")
            if validation_stats.get('validation_success', False):
                print("✅ Validation successful (≥80% plot capture)")
            else:
                print("⚠️  Validation needs improvement (<80% plot capture)")
    
    print(f"\n✅ Coastal forest bioregion created successfully!")
    print(f"   Output: {output_file}")
    print(f"   Summary: {summary_file}")

if __name__ == "__main__":
    main()