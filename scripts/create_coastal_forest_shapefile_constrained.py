#!/usr/bin/env python3
"""
Create Coastal Forest Bioregion using Shapefile Boundary + Raster Constraints

This script creates a coastal forest bioregion by:
1. Loading user-defined coastal boundary shapefile as base area
2. Applying elevation constraints within the boundary
3. Applying tree cover constraints (≥5%) within the boundary
4. Creating natural boundaries following topography and forest coverage
5. Validating against shore pine plot locations

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

# Configuration - Coastal Forest Bioregion Parameters
MAX_ELEVATION_FT = 1000       # Broader coastal forest zone (was 500)
MIN_TREE_COVER_PCT = 10       # Forest definition threshold (was 5)
OUTPUT_RESOLUTION = 120       # Match elevation data resolution

# Performance optimization settings (Latest 2025)
min_area_sqkm = 0.5           # Larger minimum for fewer polygons
skip_merging = True           # Keep individual polygons for performance
inward_buffer = 0             # No buffer - use precise boundary clipping instead

# Geometry simplification settings (Latest 2025)
simplification_tolerance = 0.001      # ~100m simplification for web performance
coordinate_precision_decimal_places = 3  # 3 decimal places = ~100m precision

# Species validation
SHORE_PINE_SPCD = 108         # Lodgepole pine (includes shore pine)

# File paths
BASE_BOUNDARY_SHP = "/Users/JJC/trees/outputs/bioregions/coastal_forest_broad.shp"
ELEVATION_RASTER = "outputs/mapbox_masks/pnw_elevation_120m_mapbox.tif"
# Use full CONUS dataset for dynamic clipping
TREE_COVER_RASTER = "data/raw/nlcd_tcc_conus_2021_v2021-4.tif"
PLOT_DATA = "outputs/plot_carbon_percentiles_latest_surveys.geojson"
OUTPUT_DIR = Path("outputs/bioregions")

def convert_elevation_to_feet(elevation_meters):
    """Convert elevation from meters to feet"""
    return elevation_meters * 3.28084

def load_base_boundary():
    """Load the user-defined coastal boundary shapefile"""
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
            print("3. Draw the coastal forest boundary you want")
            print("4. Save it as a shapefile to the path above")
            print("\nAlternatively, you can use the longitude cutoff approach with the v2 script")
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

def create_clipped_nlcd_if_needed(boundary_gdf):
    """Create clipped NLCD file for boundary if it doesn't exist"""
    clipped_nlcd_path = "outputs/mapbox_masks/nlcd_tcc_coastal_forest_extent.tif"
    
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
        
        # Apply elevation constraint - remove areas >1000ft (data already in feet)
        elev_data_ft = elev_data  # Data is already in feet  
        elevation_mask = (elev_data_ft <= MAX_ELEVATION_FT) & (~np.isnan(elev_data_ft))
        
        print(f"Elevation constraint: {np.sum(elevation_mask)} pixels ≤ {MAX_ELEVATION_FT}ft within boundary")
        
        # Debug: Check elevation values where mask is True
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
    coastal_forest_mask = elevation_mask & tree_cover_mask
    
    print(f"Final coastal forest mask: {np.sum(coastal_forest_mask)} pixels")
    print(f"Approximate area: {np.sum(coastal_forest_mask) * (OUTPUT_RESOLUTION**2) / 1000000:.0f} km²")
    
    return coastal_forest_mask, elev_transform, elev_crs, elev_profile

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
        [{'coastal_forest': 1} for _ in polygons],
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

def remove_offshore_bands_and_simplify(gdf, boundary_gdf, min_area_sqkm=0.08, max_perimeter_area_ratio=None):
    """Remove thin offshore bands and simplify geometry with latest optimizations"""
    if min_area_sqkm is None:
        min_area_sqkm = 0.5  # Use default minimum area
    
    print(f"Removing offshore bands and simplifying (min_area={min_area_sqkm}km², max_P/A_ratio={max_perimeter_area_ratio})...")
    """Remove thin offshore bands and simplify geometry"""
    print(f"Removing offshore bands and simplifying (min_area={min_area_sqkm}km², max_P/A_ratio={max_perimeter_area_ratio})...")
    
    if gdf.empty:
        return gdf
    
    # Apply two-step clipping first
    gdf = clip_to_original_boundary(gdf, boundary_gdf)
    
    if gdf.empty:
        return gdf
    
    # Convert to geographic CRS for area calculation
    gdf_geo = gdf.to_crs('EPSG:4326')
    
    # Skip water exclusion for now - was removing valid coastal forests
    print("Skipping water exclusion to preserve coastal forest areas")
    
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
        print("Warning: No polygons meet minimum area requirement - using all remaining polygons")
        large_polygons = gdf_geo
    
    # Skip perimeter-to-area ratio filter for now
    if max_perimeter_area_ratio is not None:
        compact_polygons = large_polygons[large_polygons['perimeter_area_ratio'] <= max_perimeter_area_ratio]
        print(f"After compactness filter: {len(compact_polygons)} polygons with P/A ratio ≤ {max_perimeter_area_ratio}")
        
        if compact_polygons.empty:
            print("Warning: No polygons meet compactness requirement, using area-filtered polygons")
            compact_polygons = large_polygons
    else:
        compact_polygons = large_polygons
        print("Skipping perimeter/area ratio filter - keeping all polygons above minimum area")
    
    # Simplify geometry and remove small holes for web use
    tolerance = simplification_tolerance  # Use global setting (~100m)
    hole_threshold = 0.001  # Remove holes smaller than ~100m radius (0.03 km²)
    
    print(f"Simplifying geometries (tolerance={tolerance}) and removing small holes...")
    
    def simplify_and_remove_holes(geom):
        """Simplify geometry and remove small holes"""
        if geom is None or geom.is_empty:
            return geom
        
        # First simplify
        simplified = geom.simplify(tolerance)
        
        # Remove small holes from polygons
        if hasattr(simplified, 'exterior'):  # Single Polygon
            exterior = simplified.exterior
            holes = [hole for hole in simplified.interiors 
                    if Polygon(hole).area >= hole_threshold]
            return Polygon(exterior, holes)
        elif hasattr(simplified, 'geoms'):  # MultiPolygon
            cleaned_polys = []
            for poly in simplified.geoms:
                if hasattr(poly, 'exterior'):
                    exterior = poly.exterior
                    holes = [hole for hole in poly.interiors 
                            if Polygon(hole).area >= hole_threshold]
                    cleaned_polys.append(Polygon(exterior, holes))
            return MultiPolygon(cleaned_polys) if cleaned_polys else simplified
        else:
            return simplified
    
    compact_polygons['geometry'] = compact_polygons.geometry.apply(simplify_and_remove_holes)
    
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
    coastal_bioregion = gpd.GeoDataFrame([
        {'region_name': 'Coastal Forest'}
    ], geometry=[final_geometry], crs=compact_polygons.crs)
    
    print(f"Created single MultiPolygon feature for Coastal Forest")
    
    return coastal_bioregion

def validate_with_shore_pine_plots(bioregion_gdf):
    """Validate bioregion against shore pine plot locations"""
    print("Validating bioregion against shore pine plots...")
    
    try:
        # Load plot data
        plots_gdf = gpd.read_file(PLOT_DATA)
        
        # Filter for shore pine plots (Pinus contorta)
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
        
        print(f"Found {len(coastal_pine_plots)} shore pine plots in elevation range")
        
        if len(coastal_pine_plots) == 0:
            return {'validation': 'no_plots_found'}
        
        # Check how many plots fall within bioregion
        plots_in_bioregion = gpd.sjoin(coastal_pine_plots, bioregion_gdf, how='inner', predicate='within')
        
        validation_stats = {
            'total_shore_pine_plots': len(shore_pine_plots),
            'coastal_shore_pine_plots': len(coastal_pine_plots),
            'plots_captured_by_bioregion': len(plots_in_bioregion),
            'capture_rate': len(plots_in_bioregion) / len(coastal_pine_plots) if len(coastal_pine_plots) > 0 else 0,
            'validation_success': len(plots_in_bioregion) >= len(coastal_pine_plots) * 0.7  # 70% capture rate
        }
        
        print(f"Validation results:")
        print(f"  Total shore pine plots: {validation_stats['total_shore_pine_plots']}")
        print(f"  Coastal shore pine plots: {validation_stats['coastal_shore_pine_plots']}")
        print(f"  Plots captured by bioregion: {validation_stats['plots_captured_by_bioregion']}")
        print(f"  Capture rate: {validation_stats['capture_rate']:.1%}")
        
        # Save validation plots for inspection
        if len(coastal_pine_plots) > 0:
            plot_output = OUTPUT_DIR / "shore_pine_validation_plots_shapefile.geojson"
            coastal_pine_plots.to_file(plot_output, driver='GeoJSON')
            print(f"Saved validation plots to {plot_output}")
        
        return validation_stats
        
    except Exception as e:
        print(f"Warning: Could not validate with plot data: {e}")
        return {'validation': 'error', 'error': str(e)}

def main():
    """Main processing function"""
    print("=== Creating Coastal Forest Bioregion (Shapefile + Raster Constraints) ===")
    
    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Load base boundary shapefile
    boundary_gdf = load_base_boundary()
    if boundary_gdf is None:
        print("Error: Could not load base boundary shapefile")
        return
    
    # Step 2: Create clipped NLCD first
    clipped_nlcd_path = create_clipped_nlcd_if_needed(boundary_gdf)
    
    # Apply raster constraints within boundary
    mask_array, transform, crs, profile = apply_raster_constraints_within_boundary(boundary_gdf, clipped_nlcd_path)
    if mask_array is None:
        print("Error: Could not apply raster constraints")
        return
    
    # Step 3: Convert raster to vector polygons
    coastal_polygons = raster_to_polygons_optimized(mask_array, transform, crs)
    
    # Step 4: Remove offshore bands and simplify
    coastal_bioregion = remove_offshore_bands_and_simplify(coastal_polygons, boundary_gdf)
    
    if coastal_bioregion.empty:
        print("Error: No coastal forest bioregion created")
        return
    
    
    # Step 6: Save outputs
    output_file = OUTPUT_DIR / "coastal_forest_constrained.geojson"
    print(f"\nSaving coastal forest bioregion to {output_file}")
    coastal_bioregion.to_file(output_file, driver='GeoJSON')
    
    # Save processing summary
    summary = {
        'method': 'shapefile_boundary_with_raster_constraints',
        'parameters': {
            'max_elevation_ft': MAX_ELEVATION_FT,
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
            'tree_cover': TREE_COVER_RASTER,
            'validation': PLOT_DATA
        },
        'created': pd.Timestamp.now().isoformat()
    }
    
    summary_file = OUTPUT_DIR / "coastal_forest_shapefile_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    
    # Print final summary
    print("\n=== Final Summary ===")
    if not coastal_bioregion.empty:
        area_sqkm = coastal_bioregion.iloc[0].get('total_area_sqkm', coastal_bioregion.geometry.area.sum() * 111 * 111)
        print(f"Coastal forest bioregion area: {area_sqkm:.0f} km²")
        
    
    print(f"\n✅ Coastal forest bioregion created successfully!")
    print(f"   Output: {output_file}")
    print(f"   Summary: {summary_file}")

if __name__ == "__main__":
    main()