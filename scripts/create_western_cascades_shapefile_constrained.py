#!/usr/bin/env python3
"""
Create Western Cascades Bioregion using Shapefile Boundary + Raster Constraints

This script creates a Western Cascades bioregion by:
1. Loading user-defined Western Cascades boundary shapefile as base area
2. Applying elevation constraints within the boundary (≥500ft)
3. Applying tree cover constraints (20-100%) within the boundary
4. Creating natural boundaries following topography and forest coverage

This approach provides much better control than pure raster-based boundary creation.
"""

import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.geometry import Point, Polygon, MultiPolygon, GeometryCollection
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

# Configuration - Western Cascades Bioregion Parameters
MIN_ELEVATION_FT = 500      # Western Cascades minimum elevation
MAX_ELEVATION_FT = 6000     # Western Cascades maximum elevation
MIN_TREE_COVER_PCT = 20     # Minimum tree cover for forest areas
MAX_TREE_COVER_PCT = 100    # Maximum tree cover (full forest)
OUTPUT_RESOLUTION = 120     # Match elevation data resolution

# Processing parameters
min_area_sqkm = 0.5          # Minimum polygon area to remove fragments
inward_buffer = 0            # No buffer - use precise boundary clipping instead

# File paths
BASE_BOUNDARY_SHP = "/Users/JJC/trees/outputs/bioregions/western_cascades_broad.shp"
ELEVATION_RASTER = "outputs/mapbox_masks/pnw_elevation_120m_mapbox.tif"
TREE_COVER_RASTER = "data/raw/nlcd_tcc_conus_2021_v2021-4.tif"
OUTPUT_DIR = Path("outputs/bioregions")

def load_base_boundary():
    """Load the user-defined Western Cascades boundary shapefile"""
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
            print("3. Draw the Western Cascades boundary you want")
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

def apply_raster_constraints_within_boundary(boundary_gdf, clipped_nlcd_path):
    """Apply elevation and tree cover constraints within the boundary"""
    print("=== Applying Raster Constraints Within Boundary ===")
    
    # Get boundary geometry for masking
    boundary_geom = boundary_gdf.geometry.iloc[0]
    
    # Buffer inward to avoid edge artifacts
    boundary_geom = boundary_geom.buffer(inward_buffer)
    print(f"Applied {abs(inward_buffer)*111:.0f}m inward buffer to avoid edge artifacts")
    
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
        
        # Apply elevation constraint - data already in feet, keep areas within range
        elev_data_ft = elev_data  # Data is already in feet  
        elevation_mask = (elev_data_ft >= MIN_ELEVATION_FT) & (elev_data_ft <= MAX_ELEVATION_FT) & (~np.isnan(elev_data_ft))
        
        print(f"Elevation constraint: {np.sum(elevation_mask)} pixels {MIN_ELEVATION_FT}-{MAX_ELEVATION_FT}ft within boundary")
        
        # Debug: Check elevation values where mask is True
        valid_elevations = elev_data_ft[elevation_mask]
        if len(valid_elevations) > 0:
            print(f"Elevation range in kept pixels: {np.min(valid_elevations):.1f}ft to {np.max(valid_elevations):.1f}ft")
        else:
            print("Warning: No pixels meet elevation criteria")
    
    # Step 2: Load and resample tree cover data to match elevation grid
    print(f"Loading tree cover data from {TREE_COVER_RASTER}")
    
    # Use the clipped NLCD data
    print(f"Using clipped NLCD data: {clipped_nlcd_path}")
    
    with rasterio.open(clipped_nlcd_path) as tcc_src:
        try:
            # Handle CRS mismatch - reproject boundary to match tree cover data if needed
            if tcc_src.crs != elev_crs:
                print(f"Reprojecting boundary from {elev_crs} to {tcc_src.crs} for tree cover processing")
                boundary_gdf_tcc = gpd.GeoDataFrame([1], geometry=[boundary_geom], crs=elev_crs)
                boundary_gdf_tcc = boundary_gdf_tcc.to_crs(tcc_src.crs)
                boundary_geom_tcc = boundary_gdf_tcc.geometry.iloc[0]
            else:
                boundary_geom_tcc = boundary_geom
                
            # Mask tree cover data to the same extent as elevation
            if tcc_src.shape == elev_data.shape and tcc_src.transform == elev_transform:
                print("Using matching resolution tree cover data")
                tcc_masked, _ = mask(tcc_src, [boundary_geom_tcc], crop=True, nodata=tcc_src.nodata)
                tcc_data = tcc_masked[0]
            else:
                print("Resampling tree cover to match elevation grid...")
                # Initialize with NaN to ensure NoData areas are properly handled
                tcc_data = np.full_like(elev_data, np.nan, dtype=np.float32)
                
                # First mask the tree cover to boundary
                tcc_masked, tcc_transform = mask(tcc_src, [boundary_geom_tcc], crop=True, nodata=tcc_src.nodata)
                
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
        tree_cover_mask = (tcc_data >= MIN_TREE_COVER_PCT) & (tcc_data <= MAX_TREE_COVER_PCT) & valid_tcc_mask
        print(f"Tree cover constraint: {np.sum(tree_cover_mask)} pixels {MIN_TREE_COVER_PCT}-{MAX_TREE_COVER_PCT}% within boundary")
        print(f"Excluding {np.sum(~valid_tcc_mask)} pixels with NoData (likely water/non-land areas)")
    
    # Step 3: Combine constraints
    print("Combining elevation and tree cover constraints...")
    western_cascades_mask = elevation_mask & tree_cover_mask
    
    print(f"Final Western Cascades mask: {np.sum(western_cascades_mask)} pixels")
    print(f"Approximate area: {np.sum(western_cascades_mask) * (OUTPUT_RESOLUTION**2) / 1000000:.0f} km²")
    
    return western_cascades_mask, elev_transform, elev_crs, elev_profile

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
        [{'western_cascades': 1} for _ in polygons],
        geometry=polygons,
        crs=crs
    )
    
    return gdf

def clip_to_original_boundary(gdf, boundary_gdf):
    """Clip polygons using two-step process to remove edge artifacts"""
    print("Clipping results to original boundary shapefile...")
    
    if gdf.empty:
        return gdf
    
    # Ensure both datasets are in same CRS
    gdf_clipped = gdf.to_crs(boundary_gdf.crs)
    original_boundary = boundary_gdf.geometry.iloc[0]  # Get the boundary geometry
    
    print("Step 1: Initial clipping to original boundary")
    
    # Find polygons that intersect with original boundary
    intersecting_mask = gdf_clipped.geometry.intersects(original_boundary)
    intersecting_polygons = gdf_clipped[intersecting_mask]
    
    print(f"Found {len(intersecting_polygons)} polygons intersecting original boundary (of {len(gdf)} total)")
    
    if intersecting_polygons.empty:
        print("Warning: No polygons intersect with original boundary")
        return gdf_clipped
    
    # First clip: Clip each intersecting polygon to the boundary
    first_clip_geometries = []
    for idx, row in intersecting_polygons.iterrows():
        try:
            clipped_geom = row.geometry.intersection(original_boundary)
            # Only keep non-empty geometries
            if not clipped_geom.is_empty and clipped_geom.area > 0:
                first_clip_geometries.append((row.to_dict(), clipped_geom))
        except Exception as e:
            print(f"Warning: Error in first clip for polygon {idx}: {e}")
            continue
    
    if not first_clip_geometries:
        print("Warning: No valid geometries after first clip")
        return gpd.GeoDataFrame()
    
    print(f"First clip completed: {len(first_clip_geometries)} polygons")
    
    print("Step 2: Inward buffer and re-clip to remove edge artifacts")
    
    # Create inward buffered boundary (50m buffer in degrees ≈ 0.00045)
    buffer_degrees = -0.00045  # Approximately -50m at this latitude
    buffered_boundary = original_boundary.buffer(buffer_degrees)
    
    print(f"Applied {abs(buffer_degrees * 111000):.0f}m inward buffer to remove edge artifacts")
    
    # Second clip: Re-clip to buffered boundary and apply area filter
    final_geometries = []
    edge_artifacts_removed = 0
    
    for data, geom in first_clip_geometries:
        try:
            # Check if geometry intersects with buffered boundary
            if geom.intersects(buffered_boundary):
                final_clipped_geom = geom.intersection(buffered_boundary)
                
                # Only keep non-empty geometries with sufficient area
                if not final_clipped_geom.is_empty and final_clipped_geom.area > 0:
                    # Convert area to square meters and apply size filter
                    area_m2 = final_clipped_geom.area * (111000**2)
                    if area_m2 >= 2000:  # Only keep parts >= 2000 sq meters
                        final_geometries.append((data, final_clipped_geom))
                    else:
                        edge_artifacts_removed += 1
                else:
                    edge_artifacts_removed += 1
            else:
                edge_artifacts_removed += 1  # Polygon was entirely in the edge zone
        except Exception as e:
            print(f"Warning: Error in second clip: {e}")
            edge_artifacts_removed += 1
            continue
    
    if edge_artifacts_removed > 0:
        print(f"Removed {edge_artifacts_removed} edge artifacts and small fragments")
    
    # Create final clipped geometries list for the common code below
    clipped_geometries = final_geometries
    
    if not clipped_geometries:
        print("Warning: No valid clipped geometries")
        return gpd.GeoDataFrame()
    
    # Create new GeoDataFrame with clipped geometries
    data_list = []
    geom_list = []
    for data, geom in clipped_geometries:
        data_copy = data.copy()
        data_copy.pop('geometry', None)  # Remove geometry from data dict
        data_list.append(data_copy)
        geom_list.append(geom)
    
    clipped_gdf = gpd.GeoDataFrame(data_list, geometry=geom_list, crs=boundary_gdf.crs)
    
    print(f"Successfully clipped to {len(clipped_gdf)} polygons within original boundary")
    return clipped_gdf

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
    
    # Simplify geometry for web use (consistent with bioregion combination)
    tolerance = 0.001  # ~100m simplification
    large_polygons['geometry'] = large_polygons.geometry.simplify(tolerance)
    
    # Round coordinates to 3 decimal places (~100m precision)
    from shapely.ops import transform
    def round_coordinates(geom):
        """Round coordinates to 3 decimal places (~100m precision)"""
        def round_coords(x, y, z=None):
            rounded_x = round(x, 3)
            rounded_y = round(y, 3)
            return (rounded_x, rounded_y) if z is None else (rounded_x, rounded_y, round(z, 3))
        return transform(round_coords, geom)
    
    large_polygons['geometry'] = large_polygons['geometry'].apply(round_coordinates)
    
    # Skip merging - keep individual polygons for faster processing
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
    western_cascades_bioregion = gpd.GeoDataFrame([
        {'region_name': 'Western Cascades'}
    ], geometry=[final_geometry], crs=large_polygons.crs)
    
    print(f"Created single MultiPolygon feature for Western Cascades")
    
    return western_cascades_bioregion

def create_clipped_nlcd_if_needed(boundary_gdf):
    """Create a clipped version of NLCD data for the boundary extent if it doesn't exist"""
    clipped_nlcd_path = "outputs/mapbox_masks/nlcd_tcc_western_cascades_extent.tif"
    
    # Check if clipped version already exists
    if Path(clipped_nlcd_path).exists():
        print(f"Using existing clipped NLCD data: {clipped_nlcd_path}")
        return clipped_nlcd_path
    
    print("Creating clipped NLCD tree cover data for boundary extent...")
    
    # Ensure output directory exists
    Path(clipped_nlcd_path).parent.mkdir(parents=True, exist_ok=True)
    
    # Get boundary geometry with small buffer for edge coverage
    boundary_geom = boundary_gdf.geometry.iloc[0]
    boundary_buffered = boundary_geom.buffer(0.01)  # ~1km buffer for edge coverage
    
    try:
        with rasterio.open(TREE_COVER_RASTER) as src:
            print(f"Clipping {TREE_COVER_RASTER} to boundary extent...")
            print(f"NLCD CRS: {src.crs}")
            print(f"Boundary CRS: {boundary_gdf.crs}")
            
            # Reproject boundary to match NLCD CRS
            boundary_reproj = boundary_gdf.to_crs(src.crs)
            boundary_geom_reproj = boundary_reproj.geometry.iloc[0]
            boundary_buffered_reproj = boundary_geom_reproj  # No buffer - precise clipping
            
            # Mask and crop to boundary
            masked_data, masked_transform = mask(
                src, 
                [boundary_buffered_reproj], 
                crop=True, 
                nodata=src.nodata
            )
            
            # Update profile for output
            output_profile = src.profile.copy()
            output_profile.update({
                'height': masked_data.shape[1],
                'width': masked_data.shape[2],
                'transform': masked_transform
            })
            
            # Write clipped data
            with rasterio.open(clipped_nlcd_path, 'w', **output_profile) as dst:
                dst.write(masked_data)
            
            print(f"Created clipped NLCD data: {clipped_nlcd_path}")
            print(f"Clipped size: {masked_data.shape[1]} x {masked_data.shape[2]} pixels")
            
    except Exception as e:
        print(f"Error creating clipped NLCD data: {e}")
        print(f"Falling back to original file: {TREE_COVER_RASTER}")
        return TREE_COVER_RASTER
    
    return clipped_nlcd_path

def main():
    """Main processing function"""
    print("=== Creating Western Cascades Bioregion (Shapefile + Raster Constraints) ===")
    
    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Load base boundary shapefile
    boundary_gdf = load_base_boundary()
    if boundary_gdf is None:
        print("Error: Could not load base boundary shapefile")
        return
    
    # Step 1.5: Create clipped NLCD data for efficient processing
    clipped_nlcd_path = create_clipped_nlcd_if_needed(boundary_gdf)
    
    # Step 2: Apply raster constraints within boundary
    mask_array, transform, crs, profile = apply_raster_constraints_within_boundary(boundary_gdf, clipped_nlcd_path)
    if mask_array is None:
        print("Error: Could not apply raster constraints")
        return
    
    # Step 3: Convert raster to vector polygons
    western_cascades_polygons = raster_to_polygons_optimized(mask_array, transform, crs)
    
    # Step 4: Clip to original boundary to remove areas outside shapefile
    clipped_polygons = clip_to_original_boundary(western_cascades_polygons, boundary_gdf)
    
    # Step 5: Remove artifacts and simplify
    western_cascades_bioregion = remove_artifacts_and_simplify(clipped_polygons)
    
    if western_cascades_bioregion.empty:
        print("Error: No Western Cascades bioregion created")
        return
    
    # Step 6: Save outputs
    output_file = OUTPUT_DIR / "western_cascades_constrained.geojson"
    print(f"\nSaving Western Cascades bioregion to {output_file}")
    western_cascades_bioregion.to_file(output_file, driver='GeoJSON')
    
    # Save processing summary
    summary = {
        'method': 'shapefile_boundary_with_raster_constraints',
        'parameters': {
            'min_elevation_ft': MIN_ELEVATION_FT,
            'max_elevation_ft': MAX_ELEVATION_FT,
            'min_tree_cover_pct': MIN_TREE_COVER_PCT,
            'max_tree_cover_pct': MAX_TREE_COVER_PCT,
            'output_resolution_m': OUTPUT_RESOLUTION,
            'min_area_sqkm': min_area_sqkm,
            'inward_buffer_degrees': inward_buffer,
            'simplification_tolerance': 0.001,
            'coordinate_precision_decimal_places': 3
        },
        'data_sources': {
            'base_boundary': BASE_BOUNDARY_SHP,
            'elevation': ELEVATION_RASTER,
            'tree_cover': TREE_COVER_RASTER
        },
        'created': pd.Timestamp.now().isoformat()
    }
    
    summary_file = OUTPUT_DIR / "western_cascades_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    
    # Print final summary
    print("\n=== Final Summary ===")
    if not western_cascades_bioregion.empty:
        polygon_count = len(western_cascades_bioregion)
        print(f"Western Cascades bioregion: {polygon_count} feature(s)")
    
    print(f"\n✅ Western Cascades bioregion created successfully!")
    print(f"   Output: {output_file}")
    print(f"   Summary: {summary_file}")

if __name__ == "__main__":
    main()