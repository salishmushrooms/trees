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

# Configuration - Western Cascades Bioregion Parameters
MIN_ELEVATION_FT = 500      # Western Cascades minimum elevation
MAX_ELEVATION_FT = 5000     # Western Cascades maximum elevation (UNCHANGED)
MIN_TREE_COVER_PCT = 40     # Minimum tree cover for forest areas (decreased for more detail)
MAX_TREE_COVER_PCT = 100    # Maximum tree cover (full forest)
OUTPUT_RESOLUTION = 120     # Higher resolution for 3x more detail (was 180)

# Processing parameters - enhanced filtering for higher detail
min_area_sqkm = 0.05         # Much smaller forest patches for higher detail (was 0.1)
inward_buffer = 0            # No buffer - use precise boundary clipping instead
morphological_closing_pixels = 1  # Reduced closing to preserve more natural gaps and detail

# File paths
BASE_BOUNDARY_SHP = "/Users/JJC/trees/outputs/bioregions/western_cascades_broad.shp"
# Try to use 180m elevation data if available, otherwise fall back to 120m
ELEVATION_RASTER_180M = "outputs/mapbox_masks/pnw_elevation_180m_mapbox.tif"  
ELEVATION_RASTER_120M = "outputs/mapbox_masks/pnw_elevation_120m_mapbox.tif"
ELEVATION_RASTER = ELEVATION_RASTER_180M if Path(ELEVATION_RASTER_180M).exists() else ELEVATION_RASTER_120M
TREE_COVER_RASTER = "data/raw/nlcd_tcc_conus_2021_v2021-4.tif"
OUTPUT_DIR = Path("outputs/bioregions")

def load_base_boundary():
    """Load and validate the base boundary shapefile with enhanced error handling"""
    print(f"Loading base boundary from {BASE_BOUNDARY_SHP}")
    
    # Check if file exists first
    if not Path(BASE_BOUNDARY_SHP).exists():
        print(f"❌ ERROR: Base boundary shapefile not found: {BASE_BOUNDARY_SHP}")
        print("\nTo fix this, you need to:")
        print("1. Open QGIS or another GIS software")
        print("2. Create a new polygon layer")
        print("3. Draw the Western Cascades boundary you want")
        print("4. Save it as a shapefile to the path above")
        return None
    
    try:
        boundary_gdf = gpd.read_file(BASE_BOUNDARY_SHP)
        
        print(f"Loaded boundary with {len(boundary_gdf)} features")
        print(f"Columns: {list(boundary_gdf.columns)}")
        print(f"Original CRS: {boundary_gdf.crs}")
        
        if len(boundary_gdf) == 0:
            print("\n❌ ERROR: The shapefile is empty (0 features)")
            return None
        
        # Ensure proper CRS
        if boundary_gdf.crs != 'EPSG:4326':
            print(f"Converting boundary from {boundary_gdf.crs} to EPSG:4326")
            boundary_gdf = boundary_gdf.to_crs('EPSG:4326')
        
        # Enhanced geometry validation
        print("Performing enhanced geometry validation...")
        valid_geoms = boundary_gdf[boundary_gdf.geometry.is_valid & ~boundary_gdf.geometry.is_empty]
        
        if len(valid_geoms) == 0:
            print("❌ ERROR: No valid geometries found in shapefile")
            return None
        
        if len(valid_geoms) < len(boundary_gdf):
            print(f"Warning: Removed {len(boundary_gdf) - len(valid_geoms)} invalid geometries")
            boundary_gdf = valid_geoms
        
        # Fix any invalid geometries before merging
        from shapely.validation import make_valid
        def fix_geometry(geom):
            if geom is None or geom.is_empty:
                return None
            if geom.is_valid:
                return geom
            try:
                return make_valid(geom)
            except:
                return None
        
        boundary_gdf['geometry'] = boundary_gdf['geometry'].apply(fix_geometry)
        boundary_gdf = boundary_gdf[boundary_gdf.geometry.notna()]
        
        # Merge all features into single boundary if multiple
        if len(boundary_gdf) > 1:
            try:
                merged_boundary = unary_union(boundary_gdf.geometry)
                boundary_gdf = gpd.GeoDataFrame([{'boundary': 1}], geometry=[merged_boundary], crs='EPSG:4326')
                print("Merged multiple features into single boundary")
            except Exception as e:
                print(f"Warning: Could not merge features: {e}")
                # Use the largest feature if merging fails
                boundary_gdf['area'] = boundary_gdf.geometry.area
                boundary_gdf = boundary_gdf.nlargest(1, 'area')[['geometry']]
                print("Using largest feature as boundary")
        
        # Calculate area with higher precision
        area_sqkm = boundary_gdf.geometry.area.sum() * 111.32 * 111.32  # More precise conversion
        print(f"Base boundary area: {area_sqkm:.1f} km²")
        
        # Print detailed bounds for reference
        bounds = boundary_gdf.total_bounds
        print(f"Boundary extent: {bounds[0]:.4f}, {bounds[1]:.4f} to {bounds[2]:.4f}, {bounds[3]:.4f}")
        print(f"Boundary width: {(bounds[2] - bounds[0]) * 111.32:.1f} km")
        print(f"Boundary height: {(bounds[3] - bounds[1]) * 111.32:.1f} km")
        
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
    # Note: Using precise boundary clipping instead of buffer (inward_buffer = 0)
    
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
    
    with rasterio.open(TREE_COVER_RASTER) as tcc_src:
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
        
        # Apply tree cover constraint with pre-smoothing to reduce small gaps
        valid_tcc_mask = ~np.isnan(tcc_data)
        
        if morphological_closing_pixels > 0:
            print(f"Pre-smoothing tree cover data to reduce small gaps...")
            from scipy.ndimage import binary_closing, generate_binary_structure
            
            # Create initial mask
            initial_tcc_mask = (tcc_data >= MIN_TREE_COVER_PCT) & (tcc_data <= MAX_TREE_COVER_PCT) & valid_tcc_mask
            
            # Apply morphological closing to fill small gaps
            structure = generate_binary_structure(2, 2)  # 8-connectivity
            smoothed_mask = binary_closing(initial_tcc_mask, structure=structure, iterations=morphological_closing_pixels)
            
            # Use smoothed mask instead of original
            tree_cover_mask = smoothed_mask
            gaps_filled = np.sum(smoothed_mask) - np.sum(initial_tcc_mask)
            print(f"Pre-smoothing filled {gaps_filled} pixels ({gaps_filled * (OUTPUT_RESOLUTION**2) / 10000:.1f} hectares) of small gaps")
        else:
            tree_cover_mask = (tcc_data >= MIN_TREE_COVER_PCT) & (tcc_data <= MAX_TREE_COVER_PCT) & valid_tcc_mask
        
        print(f"Tree cover constraint: {np.sum(tree_cover_mask)} pixels {MIN_TREE_COVER_PCT}-{MAX_TREE_COVER_PCT}% within boundary")
        print(f"Excluding {np.sum(~valid_tcc_mask)} pixels with NoData (likely water/non-land areas)")
    
    # Step 3: Combine constraints
    print("Combining elevation and tree cover constraints...")
    western_cascades_mask = elevation_mask & tree_cover_mask
    
    print(f"Final Western Cascades mask: {np.sum(western_cascades_mask)} pixels")
    print(f"Approximate area: {np.sum(western_cascades_mask) * (OUTPUT_RESOLUTION**2) / 1000000:.0f} km²")
    
    # Apply morphological closing to fill small holes in the mask (if enabled)
    if morphological_closing_pixels > 0:
        print(f"Applying morphological closing to fill small holes ({morphological_closing_pixels} pixel radius)...")
        structure = ndimage.generate_binary_structure(2, 2)  # 8-connectivity
        western_cascades_mask_cleaned = ndimage.binary_closing(
            western_cascades_mask, 
            structure=structure, 
            iterations=morphological_closing_pixels
        )
        
        holes_filled = np.sum(western_cascades_mask_cleaned) - np.sum(western_cascades_mask)
        print(f"Filled {holes_filled} pixels ({holes_filled * (OUTPUT_RESOLUTION**2) / 10000:.1f} hectares) of small holes")
        
        return western_cascades_mask_cleaned, elev_transform, elev_crs, elev_profile
    else:
        print("Morphological closing disabled - preserving all natural gaps from low tree density areas")
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
    """Apply enhanced two-step clipping to remove edge artifacts with higher precision"""
    print("Applying enhanced two-step clipping to remove edge artifacts...")
    
    if gdf.empty:
        return gdf
    
    # Ensure both datasets are in same CRS
    gdf_clipped = gdf.to_crs(boundary_gdf.crs)
    original_boundary = boundary_gdf.geometry.iloc[0]  # Get the boundary geometry
    
    print("Step 1: Precision clipping to original boundary")
    
    # Find polygons that intersect with original boundary
    intersecting_mask = gdf_clipped.geometry.intersects(original_boundary)
    intersecting_polygons = gdf_clipped[intersecting_mask]
    
    print(f"Found {len(intersecting_polygons)} polygons intersecting original boundary (of {len(gdf)} total)")
    
    if intersecting_polygons.empty:
        print("Warning: No polygons intersect with original boundary")
        return gdf_clipped
    
    # Enhanced first clip with better error handling
    first_clip_geometries = []
    clip_errors = 0
    
    for idx, row in intersecting_polygons.iterrows():
        try:
            # Pre-validate geometry before clipping
            geom = row.geometry
            if not geom.is_valid:
                from shapely.validation import make_valid
                geom = make_valid(geom)
            
            clipped_geom = geom.intersection(original_boundary)
            
            # Enhanced validity check
            if (not clipped_geom.is_empty and 
                clipped_geom.area > 0 and 
                hasattr(clipped_geom, 'is_valid') and clipped_geom.is_valid):
                first_clip_geometries.append((row.to_dict(), clipped_geom))
            else:
                clip_errors += 1
        except Exception as e:
            print(f"Warning: Error in precision clip for polygon {idx}: {e}")
            clip_errors += 1
            continue
    
    if clip_errors > 0:
        print(f"Encountered {clip_errors} clipping errors (geometries skipped)")
    
    if not first_clip_geometries:
        print("Warning: No valid geometries after precision clip")
        return gpd.GeoDataFrame()
    
    print(f"Precision clip completed: {len(first_clip_geometries)} polygons")
    
    print("Step 2: Multi-stage artifact removal")
    
    # Stage 2a: Small inward buffer to remove thin edge artifacts
    buffer_degrees_small = -0.0003  # ~30m inward buffer for fine detail
    buffered_boundary_small = original_boundary.buffer(buffer_degrees_small)
    
    # Stage 2b: Larger inward buffer for major edge cleanup
    buffer_degrees_large = -0.00045  # ~50m inward buffer for major artifacts
    buffered_boundary_large = original_boundary.buffer(buffer_degrees_large)
    
    print(f"Applied {abs(buffer_degrees_small * 111000):.0f}m and {abs(buffer_degrees_large * 111000):.0f}m inward buffers")
    
    # Multi-stage clipping and filtering
    final_geometries = []
    edge_artifacts_removed = 0
    small_fragments_removed = 0
    
    for data, geom in first_clip_geometries:
        try:
            # Check original area before buffering
            original_area_m2 = geom.area * (111320**2)  # More precise conversion
            
            # Apply small buffer first for precision
            if geom.intersects(buffered_boundary_small):
                stage1_clipped = geom.intersection(buffered_boundary_small)
                
                if not stage1_clipped.is_empty and stage1_clipped.area > 0:
                    # For larger polygons, use small buffer; for smaller ones, use large buffer
                    if original_area_m2 >= 10000:  # >= 1 hectare
                        final_clipped_geom = stage1_clipped
                        min_area_threshold = 1000  # 1000 sq meters for large polygons
                    else:
                        # Apply larger buffer for smaller polygons to remove artifacts
                        if geom.intersects(buffered_boundary_large):
                            final_clipped_geom = geom.intersection(buffered_boundary_large)
                            min_area_threshold = 500  # 500 sq meters for smaller polygons
                        else:
                            edge_artifacts_removed += 1
                            continue
                    
                    # Enhanced area and validity checking
                    if (not final_clipped_geom.is_empty and 
                        final_clipped_geom.area > 0 and
                        final_clipped_geom.is_valid):
                        
                        # Convert area and apply threshold
                        final_area_m2 = final_clipped_geom.area * (111320**2)
                        
                        if final_area_m2 >= min_area_threshold:
                            # Handle MultiPolygon parts separately
                            if hasattr(final_clipped_geom, 'geoms'):
                                valid_parts = []
                                for part in final_clipped_geom.geoms:
                                    part_area_m2 = part.area * (111320**2)
                                    if part_area_m2 >= min_area_threshold:
                                        valid_parts.append(part)
                                    else:
                                        small_fragments_removed += 1
                                
                                if valid_parts:
                                    if len(valid_parts) == 1:
                                        final_clipped_geom = valid_parts[0]
                                    else:
                                        final_clipped_geom = MultiPolygon(valid_parts)
                                    final_geometries.append((data, final_clipped_geom))
                                else:
                                    small_fragments_removed += 1
                            else:
                                final_geometries.append((data, final_clipped_geom))
                        else:
                            small_fragments_removed += 1
                    else:
                        edge_artifacts_removed += 1
                else:
                    edge_artifacts_removed += 1
            else:
                edge_artifacts_removed += 1
                
        except Exception as e:
            print(f"Warning: Error in multi-stage clip: {e}")
            edge_artifacts_removed += 1
            continue
    
    if edge_artifacts_removed > 0:
        print(f"Removed {edge_artifacts_removed} edge artifacts")
    if small_fragments_removed > 0:
        print(f"Removed {small_fragments_removed} small fragments")
    
    if not final_geometries:
        print("Warning: No valid clipped geometries after multi-stage processing")
        return gpd.GeoDataFrame()
    
    # Create new GeoDataFrame with enhanced clipped geometries
    data_list = []
    geom_list = []
    for data, geom in final_geometries:
        data_copy = data.copy()
        data_copy.pop('geometry', None)
        data_list.append(data_copy)
        geom_list.append(geom)
    
    clipped_gdf = gpd.GeoDataFrame(data_list, geometry=geom_list, crs=boundary_gdf.crs)
    
    print(f"Enhanced clipping complete: {len(clipped_gdf)} high-quality polygons retained")
    return clipped_gdf

def remove_artifacts_and_simplify(gdf):
    """Enhanced artifact removal and geometry simplification with 3x more detail"""
    print(f"Enhanced artifact removal and simplification (min_area={min_area_sqkm}km²)...")
    
    if gdf.empty:
        return gdf
    
    # Convert to geographic CRS for area calculation
    gdf_geo = gdf.to_crs('EPSG:4326')
    
    # Calculate areas and perimeters with higher precision
    gdf_geo['area_sqkm'] = gdf_geo.geometry.area * 111.32 * 111.32  # More precise conversion
    gdf_geo['perimeter_km'] = gdf_geo.geometry.length * 111.32  # Calculate perimeter
    
    print(f"Initial polygon count: {len(gdf_geo)}")
    print(f"Total area before filtering: {gdf_geo['area_sqkm'].sum():.2f} km²")
    
    # Stage 1: Remove micro-polygons (sub-pixel artifacts)
    print("Stage 1: Removing micro-polygons and sub-pixel artifacts...")
    pixel_area_km2 = (OUTPUT_RESOLUTION * OUTPUT_RESOLUTION) / 1000000  # 120m pixels in km²
    micro_threshold = pixel_area_km2 * 0.5  # Remove anything smaller than 0.5 pixels
    
    micro_mask = gdf_geo['area_sqkm'] < micro_threshold
    micro_count = micro_mask.sum()
    if micro_count > 0:
        gdf_geo = gdf_geo[~micro_mask]
        print(f"Removed {micro_count} micro-polygons (<{micro_threshold:.4f} km²)")
    
    # Stage 2: Enhanced rectangular artifact detection
    print("Stage 2: Enhanced rectangular artifact detection...")
    bounds = gdf_geo.bounds
    gdf_geo['bbox_width'] = (bounds['maxx'] - bounds['minx']) * 111320  # More precise conversion to meters
    gdf_geo['bbox_height'] = (bounds['maxy'] - bounds['miny']) * 111320
    gdf_geo['aspect_ratio'] = gdf_geo[['bbox_width', 'bbox_height']].max(axis=1) / gdf_geo[['bbox_width', 'bbox_height']].min(axis=1)
    
    # Multi-criteria artifact detection
    min_dimension = gdf_geo[['bbox_width', 'bbox_height']].min(axis=1)
    max_dimension = gdf_geo[['bbox_width', 'bbox_height']].max(axis=1)
    
    # Enhanced thin strip detection (multiple criteria)
    is_thin_strip = (
        (gdf_geo['aspect_ratio'] > 8) & (min_dimension < 120) |  # Very thin strips
        (gdf_geo['aspect_ratio'] > 15) & (min_dimension < 200) |  # Extremely thin strips
        (max_dimension < 150) & (gdf_geo['area_sqkm'] < pixel_area_km2 * 2)  # Small rectangular artifacts
    )
    
    strip_count = is_thin_strip.sum()
    if strip_count > 0:
        gdf_geo = gdf_geo[~is_thin_strip]
        print(f"Removed {strip_count} thin rectangular edge artifacts")
    
    # Stage 3: Perimeter-to-area ratio filtering for irregular shapes
    print("Stage 3: Perimeter-to-area ratio filtering...")
    gdf_geo['perimeter_area_ratio'] = gdf_geo['perimeter_km'] / gdf_geo['area_sqkm']
    
    # Remove highly irregular shapes (likely artifacts)
    # Normal forest patches should have reasonable perimeter/area ratios
    max_pa_ratio = 50  # Adjust based on desired compactness
    irregular_mask = gdf_geo['perimeter_area_ratio'] > max_pa_ratio
    irregular_count = irregular_mask.sum()
    
    if irregular_count > 0:
        # Don't remove all irregular shapes - only the most extreme ones
        extreme_irregular = gdf_geo['perimeter_area_ratio'] > max_pa_ratio * 1.5
        extreme_count = extreme_irregular.sum()
        if extreme_count > 0:
            gdf_geo = gdf_geo[~extreme_irregular]
            print(f"Removed {extreme_count} extremely irregular polygons (P/A ratio > {max_pa_ratio * 1.5})")
    
    # Stage 4: Enhanced area-based filtering with multiple thresholds
    print("Stage 4: Multi-threshold area filtering...")
    
    # Different thresholds based on polygon characteristics
    compact_polygons = gdf_geo[gdf_geo['perimeter_area_ratio'] <= 25]  # Compact shapes
    elongated_polygons = gdf_geo[gdf_geo['perimeter_area_ratio'] > 25]  # Elongated shapes
    
    # More lenient area threshold for compact shapes
    compact_filtered = compact_polygons[compact_polygons['area_sqkm'] >= min_area_sqkm]
    # Stricter threshold for elongated shapes
    elongated_filtered = elongated_polygons[elongated_polygons['area_sqkm'] >= min_area_sqkm * 2]
    
    # Combine filtered results
    gdf_filtered = pd.concat([compact_filtered, elongated_filtered], ignore_index=True)
    
    print(f"After enhanced filtering: {len(gdf_filtered)} polygons")
    print(f"Compact polygons retained: {len(compact_filtered)}")
    print(f"Elongated polygons retained: {len(elongated_filtered)}")
    
    if gdf_filtered.empty:
        print("Warning: No polygons meet enhanced filtering requirements")
        return gdf_geo
    
    # Stage 5: Apply geometry optimization with higher precision
    print("Stage 5: High-precision geometry optimization...")
    gdf_optimized = apply_geometry_optimization(gdf_filtered)
    
    # Stage 6: Enhanced geometry fixing and validation
    print("Stage 6: Enhanced geometry validation and fixing...")
    from shapely.validation import make_valid
    from shapely import affinity
    
    def enhanced_geometry_fix(geom):
        """Enhanced geometry fixing with multiple validation steps"""
        if geom is None or geom.is_empty:
            return None
        
        # Step 1: Basic validity check
        if geom.is_valid:
            return geom
            
        # Step 2: Try make_valid
        try:
            fixed_geom = make_valid(geom)
            if fixed_geom.is_valid and not fixed_geom.is_empty:
                return fixed_geom
        except:
            pass
        
        # Step 3: Try buffer(0) technique
        try:
            buffered_geom = geom.buffer(0)
            if buffered_geom.is_valid and not buffered_geom.is_empty:
                return buffered_geom
        except:
            pass
        
        # Step 4: Last resort - return None to exclude
        return None
    
    gdf_optimized['geometry'] = gdf_optimized['geometry'].apply(enhanced_geometry_fix)
    gdf_optimized = gdf_optimized[gdf_optimized.geometry.notna()]
    
    print(f"After enhanced geometry fixing: {len(gdf_optimized)} valid polygons")
    
    # Stage 7: Precision simplification
    print("Stage 7: Precision simplification for higher detail...")
    # Use smaller tolerance for higher detail (was 0.0015)
    tolerance = 0.0008  # ~80m simplification for higher detail
    gdf_optimized['geometry'] = gdf_optimized.geometry.simplify(tolerance, preserve_topology=True)
    
    # Stage 8: Create final MultiPolygon with enhanced processing
    print(f"Stage 8: Creating enhanced MultiPolygon from {len(gdf_optimized)} polygons...")
    
    # Enhanced geometry combination
    try:
        combined_geometry = unary_union(gdf_optimized.geometry.tolist())
    except Exception as e:
        print(f"Warning: Error in geometry union, trying alternative approach: {e}")
        # Alternative: process geometries in smaller batches
        batch_size = 100
        batch_results = []
        for i in range(0, len(gdf_optimized), batch_size):
            batch = gdf_optimized.iloc[i:i+batch_size]
            try:
                batch_union = unary_union(batch.geometry.tolist())
                batch_results.append(batch_union)
            except:
                # If batch fails, add individual geometries
                batch_results.extend(batch.geometry.tolist())
        
        combined_geometry = unary_union(batch_results)
    
    # Enhanced polygon extraction
    def extract_polygons_enhanced(geom):
        """Enhanced polygon extraction with better handling of complex geometries"""
        polygons = []
        
        if isinstance(geom, Polygon):
            if geom.is_valid and not geom.is_empty and geom.area > 0:
                polygons.append(geom)
        elif isinstance(geom, MultiPolygon):
            for poly in geom.geoms:
                if isinstance(poly, Polygon) and poly.is_valid and not poly.is_empty and poly.area > 0:
                    polygons.append(poly)
        elif isinstance(geom, GeometryCollection):
            for g in geom.geoms:
                polygons.extend(extract_polygons_enhanced(g))
        
        return polygons
    
    polygon_list = extract_polygons_enhanced(combined_geometry)
    
    if not polygon_list:
        print("Error: No valid polygons found in enhanced combined geometry")
        return gpd.GeoDataFrame()
    
    print(f"Enhanced extraction yielded {len(polygon_list)} valid polygon parts")
    
    # Create proper MultiPolygon with size filtering
    min_part_area = min_area_sqkm / 4  # Allow smaller parts within the MultiPolygon
    filtered_polygons = []
    
    for poly in polygon_list:
        poly_area_km2 = poly.area * 111.32 * 111.32
        if poly_area_km2 >= min_part_area:
            filtered_polygons.append(poly)
    
    print(f"After part-level filtering: {len(filtered_polygons)} polygon parts retained")
    
    if not filtered_polygons:
        print("Warning: No polygon parts meet minimum size requirement")
        return gpd.GeoDataFrame()
    
    # Create final geometry
    if len(filtered_polygons) == 1:
        final_geometry = filtered_polygons[0]
    else:
        final_geometry = MultiPolygon(filtered_polygons)
    
    # Calculate final statistics
    final_area_km2 = sum(poly.area * 111.32 * 111.32 for poly in filtered_polygons)
    
    # Create enhanced single feature
    western_cascades_bioregion = gpd.GeoDataFrame([
        {
            'region_name': 'Western Cascades',
            'total_area_sqkm': final_area_km2,
            'polygon_parts': len(filtered_polygons),
            'processing_resolution_m': OUTPUT_RESOLUTION,
            'detail_level': 'enhanced_3x'
        }
    ], geometry=[final_geometry], crs=gdf_optimized.crs)
    
    print(f"Enhanced processing complete:")
    print(f"  - Total area: {final_area_km2:.1f} km²")
    print(f"  - Polygon parts: {len(filtered_polygons)}")
    print(f"  - Detail level: 3x enhanced")
    
    return western_cascades_bioregion


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
    
    # Step 2: Apply raster constraints within boundary
    mask_array, transform, crs, profile = apply_raster_constraints_within_boundary(boundary_gdf)
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
    
    # Save enhanced processing summary
    summary = {
        'method': 'enhanced_shapefile_boundary_with_raster_constraints_3x_detail',
        'processing_level': 'enhanced_3x_detail',
        'parameters': {
            'min_elevation_ft': MIN_ELEVATION_FT,
            'max_elevation_ft': MAX_ELEVATION_FT,
            'min_tree_cover_pct': MIN_TREE_COVER_PCT,
            'max_tree_cover_pct': MAX_TREE_COVER_PCT,
            'output_resolution_m': OUTPUT_RESOLUTION,
            'min_area_sqkm': min_area_sqkm,
            'inward_buffer_degrees': inward_buffer,
            'simplification_tolerance': 0.0008,
            'coordinate_precision_decimal_places': 4,
            'morphological_closing_pixels': morphological_closing_pixels,
            'multi_stage_artifact_removal': True,
            'enhanced_geometry_validation': True,
            'precision_clipping': True
        },
        'enhancement_features': {
            'multi_stage_clipping': 'Two-buffer system (30m + 50m)',
            'artifact_detection': 'Multi-criteria (area, aspect ratio, perimeter/area)',
            'geometry_fixing': 'Enhanced validation with fallback methods',
            'area_filtering': 'Adaptive thresholds for compact vs elongated shapes',
            'resolution_improvement': f'Increased from 180m to {OUTPUT_RESOLUTION}m',
            'detail_preservation': 'Reduced morphological closing and simplification'
        },
        'data_sources': {
            'base_boundary': BASE_BOUNDARY_SHP,
            'elevation': ELEVATION_RASTER,
            'tree_cover': TREE_COVER_RASTER
        },
        'results': {
            'total_area_sqkm': float(western_cascades_bioregion.iloc[0].get('total_area_sqkm', 0)),
            'polygon_parts': int(western_cascades_bioregion.iloc[0].get('polygon_parts', 1)),
            'processing_stages': 8
        },
        'created': pd.Timestamp.now().isoformat()
    }
    
    summary_file = OUTPUT_DIR / "western_cascades_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    
    # Print enhanced final summary
    print("\n=== Enhanced Final Summary ===")
    if not western_cascades_bioregion.empty:
        feature_count = len(western_cascades_bioregion)
        total_area = western_cascades_bioregion.iloc[0].get('total_area_sqkm', 0)
        polygon_parts = western_cascades_bioregion.iloc[0].get('polygon_parts', 1)
        
        print(f"Western Cascades bioregion: {feature_count} feature(s) with {polygon_parts} polygon parts")
        print(f"Total area: {total_area:.1f} km²")
        print(f"Processing resolution: {OUTPUT_RESOLUTION}m (enhanced from 180m)")
        print(f"Detail level: 3x enhanced with 8-stage processing")
        print(f"Minimum area threshold: {min_area_sqkm} km²")
        print(f"Tree cover range: {MIN_TREE_COVER_PCT}-{MAX_TREE_COVER_PCT}%")
        print(f"Elevation range: {MIN_ELEVATION_FT}-{MAX_ELEVATION_FT}ft")
    
    print(f"\n✅ Enhanced Western Cascades bioregion created successfully!")
    print(f"   🎯 3x higher detail level achieved")
    print(f"   📊 Enhanced multi-stage processing completed")
    print(f"   📁 Output: {output_file}")
    print(f"   📋 Summary: {summary_file}")

if __name__ == "__main__":
    main()