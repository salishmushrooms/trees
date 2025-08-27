#!/usr/bin/env python3
"""
Create Eastern Cascades Forest Bioregion using Shapefile Boundary + Raster Constraints

This script creates an Eastern Cascades forest bioregion by:
1. Loading user-defined Eastern Cascades boundary shapefile as base area
2. Applying elevation constraints (≤3500ft) within the boundary
3. Applying tree cover constraints (≥20%) within the boundary
4. Creating natural boundaries following topography and forest coverage
5. Generating output suitable for Mapbox integration

Adapted from Olympic Mountains bioregion implementation for Eastern Cascades forests.
"""

import geopandas as gpd
import pandas as pd
import numpy as np
from scipy import ndimage
from shapely.geometry import Point, Polygon, MultiPolygon
from shapely.ops import unary_union
from bioregion_geometry_utils import apply_geometry_optimization
import rasterio
from rasterio.mask import mask
from rasterio.features import shapes
from rasterio.warp import calculate_default_transform, reproject, Resampling
import warnings
from pathlib import Path
from tqdm import tqdm
import json

warnings.filterwarnings('ignore')

# Configuration - Eastern Cascades Forest Bioregion Parameters
MAX_ELEVATION_FT = 5000       # Increased upper limit (was 3500ft)
MIN_TREE_COVER_PCT = 20       # Drier forest requirement (appropriate for Eastern Cascades)
OUTPUT_RESOLUTION = 180       # Optimized resolution for balance of detail and performance

# Processing parameters - balanced filtering to preserve natural gaps
min_area_sqkm = 0.1           # Allow smaller forest patches
morphological_closing_pixels = 2  # Fill holes up to ~2 pixels (360m radius) to reduce small gaps

# File paths
BASE_BOUNDARY_SHP = "/Users/JJC/trees/outputs/bioregions/eastern_cascades_broad.shp"
# Try to use 180m elevation data if available, otherwise fall back to 120m
ELEVATION_RASTER_180M = "outputs/mapbox_masks/pnw_elevation_180m_mapbox.tif"  
ELEVATION_RASTER_120M = "outputs/mapbox_masks/pnw_elevation_120m_mapbox.tif"
ELEVATION_RASTER = ELEVATION_RASTER_180M if Path(ELEVATION_RASTER_180M).exists() else ELEVATION_RASTER_120M
TREE_COVER_RASTER = "data/raw/nlcd_tcc_conus_2021_v2021-4.tif"
# High Cascades area to subtract from Eastern Cascades
HIGH_CASCADES_GEOJSON = "/Users/JJC/trees/outputs/bioregions/high_cascades_constrained.geojson"
OUTPUT_DIR = Path("outputs/bioregions")

def load_base_boundary():
    """Load and validate the base boundary shapefile"""
    print(f"Loading base boundary from {BASE_BOUNDARY_SHP}")
    
    if not Path(BASE_BOUNDARY_SHP).exists():
        print(f"Error: Base boundary shapefile not found: {BASE_BOUNDARY_SHP}")
        return None
    
    boundary = gpd.read_file(BASE_BOUNDARY_SHP)
    
    print(f"Loaded boundary with {len(boundary)} features")
    print(f"Columns: {list(boundary.columns)}")
    print(f"CRS: {boundary.crs}")
    
    # Ensure proper CRS
    if boundary.crs != 'EPSG:4326':
        print(f"Converting boundary from {boundary.crs} to EPSG:4326")
        boundary = boundary.to_crs('EPSG:4326')
    
    # Validate geometry
    boundary_gdf = boundary.copy()
    
    # Ensure we have valid, non-empty geometries
    valid_geoms = boundary_gdf[boundary_gdf.geometry.is_valid & ~boundary_gdf.geometry.is_empty]
    
    if len(valid_geoms) < len(boundary_gdf):
        print(f"Warning: Removed {len(boundary_gdf) - len(valid_geoms)} invalid geometries")
        boundary_gdf = valid_geoms
    
    # Merge all features into single boundary
    if len(boundary_gdf) > 1:
        try:
            merged_boundary = unary_union(boundary_gdf.geometry)
            boundary_gdf = gpd.GeoDataFrame([{'boundary': 1}], geometry=[merged_boundary], crs='EPSG:4326')
            print("Merged multiple features into single boundary")
        except Exception as e:
            print(f"Warning: Could not merge features: {e}")
    
    # Calculate area
    area_sqkm = boundary_gdf.geometry.area.sum() * 111 * 111  # Approximate conversion
    print(f"Base boundary area: {area_sqkm:.0f} km²")
    
    # Print bounds for reference
    bounds = boundary_gdf.total_bounds
    print(f"Boundary extent: {bounds[0]:.3f}, {bounds[1]:.3f} to {bounds[2]:.3f}, {bounds[3]:.3f}")
    
    return boundary_gdf

def subtract_high_cascades_area(eastern_boundary_gdf):
    """Subtract High Cascades area from Eastern Cascades boundary"""
    print(f"Loading High Cascades area to subtract from {HIGH_CASCADES_GEOJSON}")
    
    if not Path(HIGH_CASCADES_GEOJSON).exists():
        print(f"Warning: High Cascades file not found: {HIGH_CASCADES_GEOJSON}")
        print("Proceeding without subtraction...")
        return eastern_boundary_gdf
    
    try:
        # Load High Cascades area
        high_cascades_gdf = gpd.read_file(HIGH_CASCADES_GEOJSON)
        high_cascades_gdf = high_cascades_gdf.to_crs('EPSG:4326')  # Ensure same CRS
        
        print(f"Loaded High Cascades: {len(high_cascades_gdf)} features")
        
        # Get the geometries
        eastern_geom = eastern_boundary_gdf.geometry.iloc[0]
        high_cascades_geom = unary_union(high_cascades_gdf.geometry)
        
        # Subtract High Cascades from Eastern Cascades
        subtracted_geom = eastern_geom.difference(high_cascades_geom)
        
        # Calculate areas for reporting
        original_area = eastern_geom.area * 111 * 111
        high_cascades_area = high_cascades_geom.area * 111 * 111
        remaining_area = subtracted_geom.area * 111 * 111
        
        print(f"Original Eastern Cascades area: {original_area:.0f} km²")
        print(f"High Cascades area to subtract: {high_cascades_area:.0f} km²")
        print(f"Remaining Eastern Cascades area: {remaining_area:.0f} km²")
        print(f"Subtracted: {original_area - remaining_area:.0f} km² ({((original_area - remaining_area)/original_area)*100:.1f}%)")
        
        # Create new GeoDataFrame with subtracted geometry
        subtracted_gdf = gpd.GeoDataFrame([{'id': 1}], geometry=[subtracted_geom], crs='EPSG:4326')
        
        return subtracted_gdf
        
    except Exception as e:
        print(f"Error subtracting High Cascades: {e}")
        print("Proceeding with original Eastern Cascades boundary...")
        return eastern_boundary_gdf

def apply_raster_constraints_within_boundary(boundary_gdf):
    """Apply elevation and tree cover constraints within the boundary"""
    print("=== Applying Raster Constraints Within Boundary ===")
    
    # Get boundary geometry for masking
    boundary_geom = boundary_gdf.geometry.iloc[0]
    
    # Apply buffer (Latest 2025: no buffer, use two-step clipping instead)
    boundary_geom = boundary_geom.buffer(0)  
    print("Applied 0m buffer (Latest 2025: no buffer, rely on two-step clipping)")
    
    # Step 1: Load and process elevation data
    print(f"Loading elevation data from {ELEVATION_RASTER}")
    with rasterio.open(ELEVATION_RASTER) as elev_src:
        # Mask elevation data to boundary
        try:
            elev_masked, elev_transform = mask(elev_src, [boundary_geom], crop=True, nodata=elev_src.nodata)
            elev_data = elev_masked[0]  # Get first band
            elev_crs = elev_src.crs
            elev_profile = elev_src.profile.copy()
            elev_profile.update({
                'height': elev_data.shape[0],
                'width': elev_data.shape[1],
                'transform': elev_transform
            })
        except Exception as e:
            print(f"Error masking elevation data: {e}")
            return None, None, None, None
        
        # Apply elevation constraint - Eastern Cascades (data already in feet)
        elev_data_ft = elev_data  # Data is already in feet
        # No minimum elevation for Eastern Cascades
        elevation_mask = (elev_data_ft <= MAX_ELEVATION_FT) & (~np.isnan(elev_data_ft))
        
        print(f"Elevation constraint: {np.sum(elevation_mask)} pixels ≤ {MAX_ELEVATION_FT}ft")
        
        # Debug: Check elevation values where mask is True
        if np.sum(elevation_mask) > 0:
            valid_elevations = elev_data_ft[elevation_mask]
            print(f"Elevation range in kept pixels: {np.min(valid_elevations):.1f}ft to {np.max(valid_elevations):.1f}ft")
    
    # Step 2: Load and resample tree cover data to match elevation grid
    print(f"Loading tree cover data from {TREE_COVER_RASTER}")
    tcc_path = Path(TREE_COVER_RASTER)
    
    with rasterio.open(tcc_path) as tcc_src:
        try:
            # Mask tree cover data to the same extent as elevation
            if tcc_src.shape == elev_data.shape and tcc_src.transform == elev_transform:
                print("Using matching resolution tree cover data")
                tcc_masked, _ = mask(tcc_src, [boundary_geom], crop=True, nodata=tcc_src.nodata)
                tcc_data = tcc_masked[0]
            else:
                print("Resampling tree cover to match elevation grid...")
                # Initialize with NaN to ensure NoData areas are properly handled
                tcc_data = np.full_like(elev_data, np.nan, dtype=np.float32)
                
                # First mask the tree cover to boundary
                tcc_masked, tcc_transform = mask(tcc_src, [boundary_geom], crop=True, nodata=tcc_src.nodata)
                
                # Then resample to elevation grid
                reproject(
                    source=tcc_masked,
                    destination=tcc_data,
                    src_transform=tcc_transform,
                    src_crs=tcc_src.crs,
                    dst_transform=elev_transform, 
                    dst_crs=elev_crs,
                    resampling=Resampling.average,
                    dst_nodata=np.nan
                )
        except Exception as e:
            print(f"Error processing tree cover data: {e}")
            return None, None, None, None
        
        # Apply tree cover constraint with pre-smoothing to reduce small gaps
        valid_tcc_mask = ~np.isnan(tcc_data)
        
        if morphological_closing_pixels > 0:
            print(f"Pre-smoothing tree cover data to reduce small gaps...")
            from scipy.ndimage import binary_closing, generate_binary_structure
            
            # Create initial mask
            initial_tcc_mask = (tcc_data >= MIN_TREE_COVER_PCT) & valid_tcc_mask
            
            # Apply morphological closing to fill small gaps
            structure = generate_binary_structure(2, 2)  # 8-connectivity
            smoothed_mask = binary_closing(initial_tcc_mask, structure=structure, iterations=morphological_closing_pixels)
            
            # Use smoothed mask instead of original
            tree_cover_mask = smoothed_mask
            gaps_filled = np.sum(smoothed_mask) - np.sum(initial_tcc_mask)
            print(f"Pre-smoothing filled {gaps_filled} pixels ({gaps_filled * (OUTPUT_RESOLUTION**2) / 10000:.1f} hectares) of small gaps")
        else:
            tree_cover_mask = (tcc_data >= MIN_TREE_COVER_PCT) & valid_tcc_mask
        
        print(f"Tree cover constraint: {np.sum(tree_cover_mask)} pixels ≥ {MIN_TREE_COVER_PCT}% within boundary")
        print(f"Excluding {np.sum(~valid_tcc_mask)} pixels with NoData (likely water/non-land areas)")
    
    # Step 3: Combine constraints
    print("Combining elevation and tree cover constraints...")
    eastern_cascades_mask = elevation_mask & tree_cover_mask
    
    print(f"Final Eastern Cascades forest mask: {np.sum(eastern_cascades_mask)} pixels")
    print(f"Approximate area: {np.sum(eastern_cascades_mask) * (OUTPUT_RESOLUTION**2) / 1000000:.0f} km²")
    
    # Apply morphological closing to fill small holes in the mask (if enabled)
    if morphological_closing_pixels > 0:
        print(f"Applying morphological closing to fill small holes ({morphological_closing_pixels} pixel radius)...")
        structure = ndimage.generate_binary_structure(2, 2)  # 8-connectivity
        eastern_cascades_mask_cleaned = ndimage.binary_closing(
            eastern_cascades_mask, 
            structure=structure, 
            iterations=morphological_closing_pixels
        )
        
        holes_filled = np.sum(eastern_cascades_mask_cleaned) - np.sum(eastern_cascades_mask)
        print(f"Filled {holes_filled} pixels ({holes_filled * (OUTPUT_RESOLUTION**2) / 10000:.1f} hectares) of small holes")
        
        return eastern_cascades_mask_cleaned, elev_transform, elev_crs, elev_profile
    else:
        print("Morphological closing disabled - preserving all natural gaps from low tree density areas")
        return eastern_cascades_mask, elev_transform, elev_crs, elev_profile

def raster_to_polygons_optimized(mask_array, transform, crs):
    """Convert raster mask to vector polygons with optimized processing"""
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
            # Basic validity check
            if poly.is_valid and poly.area > 0:
                polygons.append(poly)
    
    print(f"Extracted {len(polygons)} valid polygon features")
    
    if not polygons:
        print("Warning: No polygons extracted from mask")
        return gpd.GeoDataFrame()
    
    # Create GeoDataFrame
    gdf = gpd.GeoDataFrame(
        [{'eastern_cascades_forest': 1} for _ in polygons],
        geometry=polygons,
        crs=crs
    )
    
    return gdf

def remove_artifacts_and_simplify(gdf, min_area_sqkm=min_area_sqkm, max_perimeter_area_ratio=None):
    """Remove artifacts and simplify geometry"""
    print(f"Removing artifacts and simplifying (min_area={min_area_sqkm}km², max_P/A_ratio={max_perimeter_area_ratio})...")
    
    if gdf.empty:
        return gdf
    
    # Convert to geographic CRS for area calculation
    gdf_geo = gdf.to_crs('EPSG:4326')
    
    # Skip water exclusion
    print("Skipping water exclusion")
    
    # Calculate areas and perimeters in km first
    gdf_geo['area_sqkm'] = gdf_geo.geometry.area * 111 * 111  # Approximate conversion
    gdf_geo['perimeter_km'] = gdf_geo.geometry.length * 111   # Approximate conversion
    
    # Remove isolated single pixels (likely edge artifacts)
    print("Removing isolated single-pixel artifacts...")
    pixel_area_km2 = (120 * 120) / 1000000  # 120m pixels in km²
    single_pixel_threshold = pixel_area_km2 * 1.5  # Remove anything smaller than 1.5 pixels
    
    isolated_mask = gdf_geo['area_sqkm'] < single_pixel_threshold
    isolated_count = isolated_mask.sum()
    if isolated_count > 0:
        gdf_geo = gdf_geo[~isolated_mask]
        print(f"Removed {isolated_count} isolated single-pixel polygons")
    
    # Remove rectangular edge artifacts (thin strips along boundaries)
    print("Removing rectangular edge artifacts...")
    # Calculate bounding box dimensions for each polygon
    bounds = gdf_geo.bounds
    gdf_geo['bbox_width'] = (bounds['maxx'] - bounds['minx']) * 111000  # Convert to meters
    gdf_geo['bbox_height'] = (bounds['maxy'] - bounds['miny']) * 111000  # Convert to meters
    gdf_geo['aspect_ratio'] = gdf_geo[['bbox_width', 'bbox_height']].max(axis=1) / gdf_geo[['bbox_width', 'bbox_height']].min(axis=1)
    
    # Remove thin rectangular strips (high aspect ratio + one dimension ~pixel size)
    min_dimension = gdf_geo[['bbox_width', 'bbox_height']].min(axis=1)
    is_thin_strip = (gdf_geo['aspect_ratio'] > 10) & (min_dimension < 150)  # <150m width and >10:1 ratio
    
    strip_count = is_thin_strip.sum()
    if strip_count > 0:
        gdf_geo = gdf_geo[~is_thin_strip]
        print(f"Removed {strip_count} thin rectangular edge artifacts")
    gdf_geo['perimeter_area_ratio'] = gdf_geo['perimeter_km'] / gdf_geo['area_sqkm']
    
    # Filter by minimum area
    large_polygons = gdf_geo[gdf_geo['area_sqkm'] >= min_area_sqkm]
    print(f"After area filter: {len(large_polygons)} polygons ≥ {min_area_sqkm}km²")
    
    if large_polygons.empty:
        print("Warning: No polygons meet minimum area requirement")
        return gdf_geo
    
    # Skip perimeter-to-area ratio filter for now
    compact_polygons = large_polygons
    print("Skipping perimeter/area ratio filter - keeping all polygons above minimum area")
    
    # Apply geometry optimization (simplify and remove holes)
    compact_polygons = apply_geometry_optimization(compact_polygons)
    
    # Try to merge polygons
    try:
        merged_geometry = unary_union(compact_polygons.geometry)
        print("Successfully merged polygons")
    except Exception as e:
        print(f"Warning: Error merging polygons: {e}")
        print("Using MultiPolygon of individual features")
        from shapely.geometry import MultiPolygon
        valid_polygons = [geom for geom in compact_polygons.geometry if geom.is_valid]
        merged_geometry = MultiPolygon(valid_polygons) if len(valid_polygons) > 1 else valid_polygons[0]
    
    # Calculate final area
    total_area = compact_polygons.geometry.area.sum() * 111 * 111  # km²
    
    # Create final GeoDataFrame
    eastern_cascades_bioregion = gpd.GeoDataFrame(
        [{
            'region_name': 'Eastern Cascades Forest',
            'region_code': 'eastern_cascades',
            'description': f'Eastern Cascades forest ≤{MAX_ELEVATION_FT}ft elevation with ≥{MIN_TREE_COVER_PCT}% tree cover',
            'elevation_max_ft': MAX_ELEVATION_FT,
            'min_tree_cover_pct': MIN_TREE_COVER_PCT,
            'area_sqkm': total_area,
            'created_date': pd.Timestamp.now().isoformat(),
            'source_boundary': BASE_BOUNDARY_SHP
        }],
        geometry=[merged_geometry],
        crs='EPSG:4326'
    )
    
    return eastern_cascades_bioregion

def main():
    """Main processing function"""
    print("=== Creating Eastern Cascades Forest Bioregion (Shapefile + Raster Constraints) ===")
    
    # Load base boundary
    boundary = load_base_boundary()
    if boundary is None:
        print("Failed to load base boundary")
        return
    
    # Subtract High Cascades area from Eastern Cascades boundary
    boundary = subtract_high_cascades_area(boundary)
    
    # Apply raster constraints
    forest_mask, transform, crs, _ = apply_raster_constraints_within_boundary(boundary)
    if forest_mask is None:
        print("Failed to apply raster constraints")
        return
    
    # Convert to polygons
    eastern_cascades_polygons = raster_to_polygons_optimized(forest_mask, transform, crs)
    
    if eastern_cascades_polygons.empty:
        print("No forest polygons extracted")
        return
    
    # Clean up and simplify
    eastern_cascades_bioregion = remove_artifacts_and_simplify(eastern_cascades_polygons)
    
    if eastern_cascades_bioregion.empty:
        print("No polygons remain after filtering")
        return
    
    
    # Save outputs
    output_file = OUTPUT_DIR / "eastern_cascades_constrained.geojson"
    summary_file = OUTPUT_DIR / "eastern_cascades_summary.json"
    
    print(f"\nSaving Eastern Cascades forest bioregion to {output_file}")
    eastern_cascades_bioregion.to_file(output_file, driver='GeoJSON')
    
    # Create summary
    summary = {
        'bioregion_name': 'Eastern Cascades Forest',
        'creation_date': pd.Timestamp.now().isoformat(),
        'parameters': {
            'base_boundary': BASE_BOUNDARY_SHP,
            'elevation_range_ft': f"0-{MAX_ELEVATION_FT}",
            'min_tree_cover_pct': MIN_TREE_COVER_PCT,
            'output_resolution_m': OUTPUT_RESOLUTION,
            'min_area_sqkm': 0.08,
            'simplification_tolerance': 0.0005
        },
        'data_sources': {
            'base_boundary': BASE_BOUNDARY_SHP,
            'elevation': ELEVATION_RASTER,
            'tree_cover': TREE_COVER_RASTER
        },
        'results': {
            'total_area_sqkm': float(eastern_cascades_bioregion.area_sqkm.iloc[0]),
            'polygon_count': len(eastern_cascades_bioregion),
        }
    }
    
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"\n=== Final Summary ===")
    print(f"Eastern Cascades forest bioregion area: {eastern_cascades_bioregion.area_sqkm.iloc[0]:.0f} km²")
    print(f"Elevation range: 0-{MAX_ELEVATION_FT}ft")
    print(f"Minimum tree cover: {MIN_TREE_COVER_PCT}%")
    
    print(f"\n✅ Eastern Cascades forest bioregion created successfully!")
    print(f"   Output: {output_file}")
    print(f"   Summary: {summary_file}")

if __name__ == "__main__":
    main()