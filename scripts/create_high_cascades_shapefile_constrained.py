#!/usr/bin/env python3
"""
Create High Cascades Forest Bioregion using Shapefile Boundary + Raster Constraints

This script creates a High Cascades forest bioregion by:
1. Loading user-defined High Cascades boundary shapefile as base area
2. Applying elevation constraints (≥3500ft) within the boundary
3. Applying tree cover constraints (≥20%) within the boundary
4. Creating natural boundaries following topography and forest coverage
5. Generating output suitable for Mapbox integration

Adapted from Olympic Mountains bioregion implementation for High Cascades subalpine forests.
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

# Configuration - High Cascades Forest Bioregion Parameters
MIN_ELEVATION_FT = 3500       # Subalpine forest lower limit
MIN_TREE_COVER_PCT = 20       # Alpine/subalpine forest requirement
OUTPUT_RESOLUTION = 120       # Match elevation data resolution

# Performance optimization settings (Latest 2025)
min_area_sqkm = 0.5           # Larger minimum for fewer polygons
skip_merging = True           # Keep individual polygons for performance
inward_buffer = 0             # No buffer - use precise boundary clipping instead

# Geometry simplification settings (Latest 2025)
simplification_tolerance = 0.001      # ~100m simplification for web performance
coordinate_precision_decimal_places = 3  # 3 decimal places = ~100m precision

# File paths
BASE_BOUNDARY_SHP = "/Users/JJC/trees/outputs/bioregions/high_cascades_broad.shp"
ELEVATION_RASTER = "outputs/mapbox_masks/pnw_elevation_120m_mapbox.tif"
# Use full CONUS dataset for dynamic clipping
TREE_COVER_RASTER = "data/raw/nlcd_tcc_conus_2021_v2021-4.tif"
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

def create_clipped_nlcd_if_needed(boundary_gdf):
    """Create clipped NLCD file for boundary if it doesn't exist"""
    clipped_nlcd_path = "outputs/mapbox_masks/nlcd_tcc_high_cascades_extent.tif"
    
    if Path(clipped_nlcd_path).exists():
        print(f"Using existing clipped NLCD: {clipped_nlcd_path}")
        return clipped_nlcd_path
    
    print(f"Creating clipped NLCD from full CONUS dataset...")
    
    # Ensure output directory exists
    Path(clipped_nlcd_path).parent.mkdir(parents=True, exist_ok=True)
    
    # Get boundary geometry and buffer slightly for clipping
    boundary_geom = boundary_gdf.geometry.iloc[0]
    
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
            clipped_data, clipped_transform = mask(
                src, 
                [boundary_buffered_reproj], 
                crop=True, 
                nodata=src.nodata
            )
            
            # Update profile for output
            profile = src.profile.copy()
            profile.update({
                'height': clipped_data.shape[1],
                'width': clipped_data.shape[2],
                'transform': clipped_transform
            })
            
            # Write clipped file
            with rasterio.open(clipped_nlcd_path, 'w', **profile) as dst:
                dst.write(clipped_data)
            
            print(f"Created clipped NLCD: {clipped_nlcd_path}")
            print(f"Clipped size: {clipped_data.shape[1]} x {clipped_data.shape[2]} pixels")
            return clipped_nlcd_path
    
    except Exception as e:
        print(f"Error creating clipped NLCD: {e}")
        print("Falling back to original NLCD file")
        return TREE_COVER_RASTER

def apply_raster_constraints_within_boundary(boundary_gdf, clipped_nlcd_path=None):
    """Apply elevation and tree cover constraints within the boundary"""
    print("=== Applying Raster Constraints Within Boundary ===")
    
    # Get boundary geometry for masking
    boundary_geom = boundary_gdf.geometry.iloc[0]
    
    # Apply buffer (Latest 2025: no buffer, use two-step clipping instead)
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
        
        # Apply elevation constraint - High Cascades (data already in feet)
        elev_data_ft = elev_data  # Data is already in feet
        # No maximum elevation for High Cascades
        elevation_mask = (elev_data_ft >= MIN_ELEVATION_FT) & (~np.isnan(elev_data_ft))
        
        print(f"Elevation constraint: {np.sum(elevation_mask)} pixels ≥ {MIN_ELEVATION_FT}ft")
        
        # Debug: Check elevation values where mask is True
        if np.sum(elevation_mask) > 0:
            valid_elevations = elev_data_ft[elevation_mask]
            print(f"Elevation range in kept pixels: {np.min(valid_elevations):.1f}ft to {np.max(valid_elevations):.1f}ft")
    
    # Step 2: Load and resample tree cover data to match elevation grid
    # Use clipped NLCD for better performance
    if clipped_nlcd_path is None:
        tcc_path = create_clipped_nlcd_if_needed(boundary_gdf)
    else:
        tcc_path = clipped_nlcd_path
    print(f"Loading tree cover data from {tcc_path}")
    tcc_path = Path(tcc_path)
    
    with rasterio.open(tcc_path) as tcc_src:
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
        tree_cover_mask = (tcc_data >= MIN_TREE_COVER_PCT) & valid_tcc_mask
        print(f"Tree cover constraint: {np.sum(tree_cover_mask)} pixels ≥ {MIN_TREE_COVER_PCT}% within boundary")
        print(f"Excluding {np.sum(~valid_tcc_mask)} pixels with NoData (likely water/non-land areas)")
    
    # Step 3: Combine constraints
    print("Combining elevation and tree cover constraints...")
    high_cascades_mask = elevation_mask & tree_cover_mask
    
    print(f"Final High Cascades forest mask: {np.sum(high_cascades_mask)} pixels")
    print(f"Approximate area: {np.sum(high_cascades_mask) * (OUTPUT_RESOLUTION**2) / 1000000:.0f} km²")
    
    return high_cascades_mask, elev_transform, elev_crs, elev_profile

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
        [{'high_cascades_forest': 1} for _ in polygons],
        geometry=polygons,
        crs=crs
    )
    
    return gdf

def clip_to_original_boundary(gdf, boundary_gdf):
    """Apply two-step clipping to remove edge artifacts (Latest 2025 approach)"""
    print("Applying two-step clipping to remove edge artifacts...")
    
    if gdf.empty:
        return gdf
    
    # Get original boundary
    original_boundary = boundary_gdf.geometry.iloc[0]
    
    clipped_polygons = []
    
    for idx, row in gdf.iterrows():
        geom = row.geometry
        
        try:
            # Step 1: Clip to original boundary (creates partial pixel artifacts)
            clipped_geom = geom.intersection(original_boundary)
            
            if clipped_geom.is_empty:
                continue
                
            # Step 2: Apply 50m inward buffer and re-clip to remove edge artifacts
            buffer_degrees = -0.00045  # ~50m inward buffer
            buffered_boundary = original_boundary.buffer(buffer_degrees)
            
            # Re-clip to buffered boundary removes all edge artifacts
            final_clipped_geom = clipped_geom.intersection(buffered_boundary)
            
            if not final_clipped_geom.is_empty and final_clipped_geom.area > 0:
                # Filter out very small pieces during clipping (2000 sq meters minimum)
                min_area_degrees2 = 2000 / (111000 * 111000)  # Convert 2000m² to degrees²
                
                if hasattr(final_clipped_geom, 'geoms'):  # MultiPolygon
                    valid_parts = [g for g in final_clipped_geom.geoms if g.area >= min_area_degrees2]
                    if valid_parts:
                        if len(valid_parts) == 1:
                            final_clipped_geom = valid_parts[0]
                        else:
                            from shapely.geometry import MultiPolygon
                            final_clipped_geom = MultiPolygon(valid_parts)
                    else:
                        continue
                elif final_clipped_geom.area < min_area_degrees2:
                    continue
                
                # Create new row with clipped geometry
                new_row = row.copy()
                new_row.geometry = final_clipped_geom
                clipped_polygons.append(new_row)
                
        except Exception as e:
            print(f"Warning: Error clipping polygon {idx}: {e}")
            continue
    
    if clipped_polygons:
        clipped_gdf = gpd.GeoDataFrame(clipped_polygons, crs=gdf.crs)
        print(f"Two-step clipping: {len(gdf)} → {len(clipped_gdf)} polygons")
        return clipped_gdf
    else:
        print("Warning: No polygons remained after two-step clipping")
        return gpd.GeoDataFrame(columns=gdf.columns, crs=gdf.crs)

def remove_artifacts_and_simplify(gdf, boundary_gdf, min_area_sqkm=0.5, max_perimeter_area_ratio=None):
    """Remove artifacts and simplify geometry with latest optimizations"""
    if min_area_sqkm is None:
        min_area_sqkm = min_area_sqkm  # Use global setting
    
    print(f"Removing artifacts and simplifying (min_area={min_area_sqkm}km², max_P/A_ratio={max_perimeter_area_ratio})...")
    """Remove artifacts and simplify geometry"""
    print(f"Removing artifacts and simplifying (min_area={min_area_sqkm}km², max_P/A_ratio={max_perimeter_area_ratio})...")
    
    if gdf.empty:
        return gdf
    
    # Apply two-step clipping first
    gdf = clip_to_original_boundary(gdf, boundary_gdf)
    
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
    
    # Simplify geometry for web use (consistent with bioregion combination)
    tolerance = simplification_tolerance  # Use global setting (~100m)
    compact_polygons['geometry'] = compact_polygons.geometry.simplify(tolerance)
    
    # Round coordinates to 3 decimal places (~100m precision)
    from shapely.ops import transform
    def round_coordinates(geom):
        """Round coordinates to 3 decimal places (~100m precision)"""
        def round_coords(x, y, z=None):
            rounded_x = round(x, coordinate_precision_decimal_places)
            rounded_y = round(y, coordinate_precision_decimal_places)
            return (rounded_x, rounded_y) if z is None else (rounded_x, rounded_y, round(z, coordinate_precision_decimal_places))
        return transform(round_coords, geom)
    
    compact_polygons['geometry'] = compact_polygons['geometry'].apply(round_coordinates)
    
    # Combine all polygons into a single MultiPolygon feature  
    print(f"Combining {len(compact_polygons)} polygons into single MultiPolygon feature...")
    
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
    
    compact_polygons['geometry'] = compact_polygons['geometry'].apply(fix_invalid_geometry)
    compact_polygons = compact_polygons[compact_polygons.geometry.notna()]
    
    # Merge all geometries into a single MultiPolygon
    combined_geometry = unary_union(compact_polygons.geometry.tolist())
    
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
    high_cascades_bioregion = gpd.GeoDataFrame([
        {'region_name': 'High Cascades'}
    ], geometry=[final_geometry], crs=compact_polygons.crs)
    
    print(f"Created single MultiPolygon feature for High Cascades")
    
    return high_cascades_bioregion

def main():
    """Main processing function"""
    print("=== Creating High Cascades Forest Bioregion (Shapefile + Raster Constraints) ===")
    
    # Load base boundary
    boundary = load_base_boundary()
    if boundary is None:
        print("Failed to load base boundary")
        return
    
    # Create clipped NLCD first
    clipped_nlcd_path = create_clipped_nlcd_if_needed(boundary)
    
    # Apply raster constraints
    forest_mask, transform, crs, _ = apply_raster_constraints_within_boundary(boundary, clipped_nlcd_path)
    if forest_mask is None:
        print("Failed to apply raster constraints")
        return
    
    # Convert to polygons
    high_cascades_polygons = raster_to_polygons_optimized(forest_mask, transform, crs)
    
    if high_cascades_polygons.empty:
        print("No forest polygons extracted")
        return
    
    # Clean up and simplify
    high_cascades_bioregion = remove_artifacts_and_simplify(high_cascades_polygons, boundary)
    
    if high_cascades_bioregion.empty:
        print("No polygons remain after filtering")
        return
    
    
    # Save outputs
    output_file = OUTPUT_DIR / "high_cascades_constrained.geojson"
    summary_file = OUTPUT_DIR / "high_cascades_summary.json"
    
    print(f"\nSaving High Cascades forest bioregion to {output_file}")
    high_cascades_bioregion.to_file(output_file, driver='GeoJSON')
    
    # Create summary
    summary = {
        'bioregion_name': 'High Cascades Forest',
        'creation_date': pd.Timestamp.now().isoformat(),
        'parameters': {
            'base_boundary': BASE_BOUNDARY_SHP,
            'elevation_range_ft': f"{MIN_ELEVATION_FT}+",
            'min_tree_cover_pct': MIN_TREE_COVER_PCT,
            'output_resolution_m': OUTPUT_RESOLUTION,
            'min_area_sqkm': min_area_sqkm,
            'simplification_tolerance': simplification_tolerance,
            'coordinate_precision_decimal_places': coordinate_precision_decimal_places,
            'skip_merging': skip_merging,
            'two_step_clipping': True
        },
        'data_sources': {
            'base_boundary': BASE_BOUNDARY_SHP,
            'elevation': ELEVATION_RASTER,
            'tree_cover': TREE_COVER_RASTER
        },
        'results': {
            'total_area_sqkm': float(high_cascades_bioregion.iloc[0].get('total_area_sqkm', high_cascades_bioregion.geometry.area.sum() * 111 * 111)),
            'polygon_count': len(high_cascades_bioregion),
        }
    }
    
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"\n=== Final Summary ===")
    total_area = high_cascades_bioregion.iloc[0].get('total_area_sqkm', high_cascades_bioregion.geometry.area.sum() * 111 * 111)
    print(f"High Cascades forest bioregion area: {total_area:.0f} km²")
    print(f"Elevation range: {MIN_ELEVATION_FT}ft+")
    print(f"Minimum tree cover: {MIN_TREE_COVER_PCT}%")
    
    print(f"\n✅ High Cascades forest bioregion created successfully!")
    print(f"   Output: {output_file}")
    print(f"   Summary: {summary_file}")

if __name__ == "__main__":
    main()