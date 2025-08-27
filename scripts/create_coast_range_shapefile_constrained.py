#!/usr/bin/env python3
"""
Create Coast Range Bioregion using Shapefile Boundary + Raster Constraints

This script creates a coast range bioregion by:
1. Loading user-defined coast range boundary shapefile as base area
2. Applying elevation constraints within the boundary (≥500ft)
3. Applying tree cover constraints (≥30%) within the boundary
4. Creating natural boundaries following topography and forest coverage

This approach provides much better control than pure raster-based boundary creation.
"""

import geopandas as gpd
import pandas as pd
import numpy as np
from scipy import ndimage
from shapely.geometry import Point, Polygon, MultiPolygon, GeometryCollection
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

# Configuration - Coast Range Bioregion Parameters
MIN_ELEVATION_FT = 500        # Coast Range minimum elevation
MIN_TREE_COVER_PCT = 30       # Forest definition threshold
OUTPUT_RESOLUTION = 240       # Reduced resolution to prevent small holes (was 120)

# Processing parameters - more aggressive filtering to prevent small holes
min_area_sqkm = 0.5           # Increased minimum area to remove fragments (was 0.08)
inward_buffer = 0             # No buffer - use precise boundary clipping instead
morphological_closing_pixels = 2  # Fill small holes in raster before vectorization

# File paths
BASE_BOUNDARY_SHP = "/Users/JJC/trees/outputs/bioregions/coastal_range_broad.shp"
# Try to use 240m elevation data if available, otherwise fall back to 120m
ELEVATION_RASTER_240M = "outputs/mapbox_masks/pnw_elevation_240m_mapbox.tif"
ELEVATION_RASTER_120M = "outputs/mapbox_masks/pnw_elevation_120m_mapbox.tif"
ELEVATION_RASTER = ELEVATION_RASTER_240M if Path(ELEVATION_RASTER_240M).exists() else ELEVATION_RASTER_120M
TREE_COVER_RASTER = "outputs/mapbox_masks/pnw_tree_cover_30m_full.tif"
OUTPUT_DIR = Path("outputs/bioregions")

def load_base_boundary():
    """Load the user-defined coast range boundary shapefile"""
    print(f"Loading base boundary from {BASE_BOUNDARY_SHP}")
    
    try:
        boundary_gdf = gpd.read_file(BASE_BOUNDARY_SHP)
        boundary_gdf = boundary_gdf.to_crs('EPSG:4326')  # Ensure WGS84
        
        print(f"Loaded boundary with {len(boundary_gdf)} features")
        print(f"Columns: {list(boundary_gdf.columns)}")
        print(f"CRS: {boundary_gdf.crs}")
        
        if len(boundary_gdf) == 0:
            print("\n❌ ERROR: The shapefile is empty (0 features)")
            print("\nTo fix this, you need to:")
            print("1. Open QGIS or another GIS software")
            print("2. Create a new polygon layer")
            print("3. Draw the coast range boundary you want")
            print("4. Save it as a shapefile to the path above")
            return None
        
        # Check for valid geometries
        valid_geoms = boundary_gdf[boundary_gdf.geometry.is_valid & ~boundary_gdf.geometry.is_empty]
        if len(valid_geoms) == 0:
            print("❌ ERROR: No valid geometries found in shapefile")
            return None
        
        if len(valid_geoms) < len(boundary_gdf):
            print(f"Warning: Using {len(valid_geoms)} valid geometries out of {len(boundary_gdf)} total")
            boundary_gdf = valid_geoms
        
        # Merge all features into single boundary if multiple
        if len(boundary_gdf) > 1:
            merged_boundary = unary_union(boundary_gdf.geometry)
            boundary_gdf = gpd.GeoDataFrame([{'boundary': 1}], geometry=[merged_boundary], crs='EPSG:4326')
            print("Merged multiple features into single boundary")
        
        # Calculate area
        area_sqkm = boundary_gdf.geometry.area.sum() * 111 * 111  # Approximate conversion
        print(f"Base boundary area: {area_sqkm:.0f} km²")
        
        # Print bounds for reference
        bounds = boundary_gdf.total_bounds
        print(f"Boundary extent: {bounds[0]:.3f}, {bounds[1]:.3f} to {bounds[2]:.3f}, {bounds[3]:.3f}")
        
        return boundary_gdf
        
    except Exception as e:
        print(f"Error loading boundary shapefile: {e}")
        print("Please ensure the shapefile exists and contains valid polygon features")
        return None

def apply_raster_constraints_within_boundary(boundary_gdf):
    """Apply elevation and tree cover constraints within the boundary"""
    print("=== Applying Raster Constraints Within Boundary ===")
    
    # Get boundary geometry for masking
    boundary_geom = boundary_gdf.geometry.iloc[0]
    
    # Buffer inward to avoid edge artifacts
    boundary_geom = boundary_geom.buffer(inward_buffer)
    print(f"Applied {abs(inward_buffer)*111:.0f}m buffer (Latest 2025: no buffer, rely on two-step clipping)")
    
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
        
        # Apply elevation constraint - data already in feet, keep areas ≥500ft
        elev_data_ft = elev_data  # Data is already in feet  
        elevation_mask = (elev_data_ft >= MIN_ELEVATION_FT) & (~np.isnan(elev_data_ft))
        
        print(f"Elevation constraint: {np.sum(elevation_mask)} pixels ≥ {MIN_ELEVATION_FT}ft within boundary")
        
        # Debug: Check elevation values where mask is True
        valid_elevations = elev_data_ft[elevation_mask]
        if len(valid_elevations) > 0:
            print(f"Elevation range in kept pixels: {np.min(valid_elevations):.1f}ft to {np.max(valid_elevations):.1f}ft")
        else:
            print("Warning: No pixels meet elevation criteria")
    
    # Step 2: Load and resample tree cover data to match elevation grid
    print(f"Loading tree cover data from {TREE_COVER_RASTER}")
    
    # Use lowest resolution tree cover data available to prevent small holes
    # Prefer 240m > 120m > 30m resolution
    tcc_240m_path = "outputs/mapbox_masks/pnw_tree_cover_240m.tif"
    tcc_120m_path = "outputs/mapbox_masks/pnw_tree_cover_120m_quarter.tif"
    
    if Path(tcc_240m_path).exists():
        tcc_path = tcc_240m_path
        print("Using 240m tree cover data to minimize holes")
    elif Path(tcc_120m_path).exists():
        tcc_path = tcc_120m_path
        print("Using 120m tree cover data")
    else:
        tcc_path = TREE_COVER_RASTER
        print("Using 30m tree cover data (will be resampled to 240m)")
    
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
        
        # Apply tree cover constraint
        valid_tcc_mask = ~np.isnan(tcc_data)
        tree_cover_mask = (tcc_data >= MIN_TREE_COVER_PCT) & valid_tcc_mask
        print(f"Tree cover constraint: {np.sum(tree_cover_mask)} pixels ≥ {MIN_TREE_COVER_PCT}% within boundary")
        print(f"Excluding {np.sum(~valid_tcc_mask)} pixels with NoData (likely water/non-land areas)")
    
    # Step 3: Combine constraints
    print("Combining elevation and tree cover constraints...")
    coast_range_mask = elevation_mask & tree_cover_mask
    
    print(f"Final coast range mask: {np.sum(coast_range_mask)} pixels")
    print(f"Approximate area: {np.sum(coast_range_mask) * (OUTPUT_RESOLUTION**2) / 1000000:.0f} km²")
    
    # Apply morphological closing to fill small holes in the mask
    print(f"Applying morphological closing to fill small holes ({morphological_closing_pixels} pixel radius)...")
    structure = ndimage.generate_binary_structure(2, 2)  # 8-connectivity
    coast_range_mask_cleaned = ndimage.binary_closing(
        coast_range_mask, 
        structure=structure, 
        iterations=morphological_closing_pixels
    )
    
    holes_filled = np.sum(coast_range_mask_cleaned) - np.sum(coast_range_mask)
    print(f"Filled {holes_filled} pixels ({holes_filled * (OUTPUT_RESOLUTION**2) / 10000:.1f} hectares) of small holes")
    
    return coast_range_mask_cleaned, elev_transform, elev_crs, elev_profile

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
        [{'coast_range': 1} for _ in polygons],
        geometry=polygons,
        crs=crs
    )
    
    return gdf

def remove_artifacts_and_simplify(gdf):
    """Remove artifacts and simplify geometry"""
    print(f"Removing artifacts and simplifying (min_area={min_area_sqkm}km²)...")
    
    if gdf.empty:
        return gdf
    
    # Convert to geographic CRS for area calculation
    gdf_geo = gdf.to_crs('EPSG:4326')
    
    # Calculate areas in km²
    gdf_geo['area_sqkm'] = gdf_geo.geometry.area * 111 * 111  # Approximate conversion
    
    # Remove isolated single pixels (likely edge artifacts)
    print("Removing isolated single-pixel artifacts...")
    pixel_area_km2 = (OUTPUT_RESOLUTION * OUTPUT_RESOLUTION) / 1000000  # 120m pixels in km²
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
    
    # Filter by minimum area
    large_polygons = gdf_geo[gdf_geo['area_sqkm'] >= min_area_sqkm]
    print(f"After area filter: {len(large_polygons)} polygons ≥ {min_area_sqkm}km²")
    
    if large_polygons.empty:
        print("Warning: No polygons meet minimum area requirement")
        return gdf_geo
    
    # Apply geometry optimization (simplify and remove holes)
    large_polygons = apply_geometry_optimization(large_polygons)
    
    # Skip merging for performance - keep individual polygons
    # Combine all polygons into a single MultiPolygon feature
    print(f"Combining {len(large_polygons)} polygons into single MultiPolygon feature...")
    
    # Fix any invalid geometries before combining
    print("Fixing invalid geometries before combining...")
    from shapely.validation import make_valid
    
    def fix_invalid_geometry(geom):
        if geom is None or geom.is_empty:
            return None
        if geom.is_valid:
            return geom
        try:
            return make_valid(geom)
        except:
            return None
    
    large_polygons['geometry'] = large_polygons['geometry'].apply(fix_invalid_geometry)
    large_polygons = large_polygons[large_polygons.geometry.notna()]
    
    # Merge all geometries into a single MultiPolygon
    combined_geometry = unary_union(large_polygons.geometry.tolist())
    
    # Ensure result is a MultiPolygon (not GeometryCollection)
    def extract_polygons(geom):
        """Extract all Polygon geometries from any geometry type"""
        if isinstance(geom, Polygon):
            return [geom]
        elif isinstance(geom, MultiPolygon):
            return list(geom.geoms)
        elif isinstance(geom, GeometryCollection):
            polygons = []
            for g in geom.geoms:
                if isinstance(g, Polygon):
                    polygons.append(g)
                elif isinstance(g, MultiPolygon):
                    polygons.extend(list(g.geoms))
            return polygons
        else:
            return []
    
    polygon_list = extract_polygons(combined_geometry)
    if not polygon_list:
        print("Error: No valid polygons found in combined geometry")
        return gpd.GeoDataFrame()
    
    # Create proper MultiPolygon
    if len(polygon_list) == 1:
        final_geometry = polygon_list[0]
    else:
        final_geometry = MultiPolygon(polygon_list)
    
    # Create single feature with only region_name for Mapbox styling
    coast_range_bioregion = gpd.GeoDataFrame([
        {'region_name': 'Coast Range'}
    ], geometry=[final_geometry], crs=large_polygons.crs)
    
    print(f"Created single MultiPolygon feature for Coast Range")
    
    return coast_range_bioregion

def main():
    """Main processing function"""
    print("=== Creating Coast Range Bioregion (Shapefile + Raster Constraints) ===")
    
    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Load base boundary shapefile
    boundary_gdf = load_base_boundary()
    if boundary_gdf is None:
        print("Error: Could not load base boundary shapefile")
        return
    
    # Step 2: Apply raster constraints within boundary
    mask_array, transform, crs, profile = apply_raster_constraints_within_boundary(boundary_gdf)
    if mask_array is None:
        print("Error: Could not apply raster constraints")
        return
    
    # Step 3: Convert raster to vector polygons
    coast_range_polygons = raster_to_polygons_optimized(mask_array, transform, crs)
    
    # Step 4: Remove artifacts and simplify
    coast_range_bioregion = remove_artifacts_and_simplify(coast_range_polygons)
    
    if coast_range_bioregion.empty:
        print("Error: No coast range bioregion created")
        return
    
    # Step 5: Save outputs
    output_file = OUTPUT_DIR / "coast_range_constrained.geojson"
    print(f"\nSaving coast range bioregion to {output_file}")
    coast_range_bioregion.to_file(output_file, driver='GeoJSON')
    
    # Save processing summary
    summary = {
        'method': 'shapefile_boundary_with_raster_constraints',
        'parameters': {
            'min_elevation_ft': MIN_ELEVATION_FT,
            'min_tree_cover_pct': MIN_TREE_COVER_PCT,
            'output_resolution_m': OUTPUT_RESOLUTION,
            'min_area_sqkm': min_area_sqkm,
            'inward_buffer_degrees': inward_buffer,
            'simplification_tolerance': 0.002,
            'coordinate_precision_decimal_places': 3,
            'morphological_closing_pixels': morphological_closing_pixels
        },
        'data_sources': {
            'base_boundary': BASE_BOUNDARY_SHP,
            'elevation': ELEVATION_RASTER,
            'tree_cover': TREE_COVER_RASTER
        },
        'created': pd.Timestamp.now().isoformat()
    }
    
    summary_file = OUTPUT_DIR / "coast_range_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    
    # Print final summary
    print("\n=== Final Summary ===")
    if not coast_range_bioregion.empty:
        # Calculate total area for summary
        total_area = coast_range_bioregion.geometry.area.sum() * 111 * 111  # Approximate conversion
        print(f"Coast Range bioregion: 1 MultiPolygon feature, {total_area:.0f} km² total area")
    
    print(f"\n✅ Coast range bioregion created successfully!")
    print(f"   Output: {output_file}")
    print(f"   Summary: {summary_file}")

if __name__ == "__main__":
    main()