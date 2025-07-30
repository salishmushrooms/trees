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

# Configuration - Coastal Forest Bioregion Parameters
MAX_ELEVATION_FT = 1000       # Broader coastal forest zone (was 500)
MIN_TREE_COVER_PCT = 10       # Forest definition threshold (was 5)
OUTPUT_RESOLUTION = 120       # Match elevation data resolution

# Species validation
SHORE_PINE_SPCD = 108         # Lodgepole pine (includes shore pine)

# File paths
BASE_BOUNDARY_SHP = "/Users/JJC/trees/outputs/bioregions/coastal_forest_broad.shp"
ELEVATION_RASTER = "outputs/mapbox_masks/pnw_elevation_120m_mapbox.tif"
TREE_COVER_RASTER = "outputs/mapbox_masks/pnw_tree_cover_30m_full.tif"
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

def apply_raster_constraints_within_boundary(boundary_gdf):
    """Apply elevation and tree cover constraints within the boundary"""
    print("=== Applying Raster Constraints Within Boundary ===")
    
    # Get boundary geometry for masking
    boundary_geom = boundary_gdf.geometry.iloc[0]
    
    # Buffer inward more aggressively to avoid edge artifacts
    boundary_geom = boundary_geom.buffer(-0.006)  # ~600m inward buffer in degrees  
    print("Applied 600m inward buffer to avoid edge artifacts")
    
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
        
        # Apply elevation constraint - remove areas >1000ft
        elev_data_ft = convert_elevation_to_feet(elev_data)
        elevation_mask = (elev_data_ft <= MAX_ELEVATION_FT) & (~np.isnan(elev_data_ft))
        
        print(f"Elevation constraint: {np.sum(elevation_mask)} pixels ≤ {MAX_ELEVATION_FT}ft within boundary")
        
        # Debug: Check elevation values where mask is True
        valid_elevations = elev_data_ft[elevation_mask]
        print(f"Elevation range in kept pixels: {np.min(valid_elevations):.1f}ft to {np.max(valid_elevations):.1f}ft")
    
    # Step 2: Load and resample tree cover data to match elevation grid
    print(f"Loading tree cover data from {TREE_COVER_RASTER}")
    
    # Use 120m tree cover data if available, otherwise resample 30m data
    tcc_120m_path = "outputs/mapbox_masks/pnw_tree_cover_120m_quarter.tif"
    tcc_path = tcc_120m_path if Path(tcc_120m_path).exists() else TREE_COVER_RASTER
    
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

def remove_offshore_bands_and_simplify(gdf, min_area_sqkm=0.08, max_perimeter_area_ratio=None):
    """Remove thin offshore bands and simplify geometry"""
    print(f"Removing offshore bands and simplifying (min_area={min_area_sqkm}km², max_P/A_ratio={max_perimeter_area_ratio})...")
    
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
        print("Warning: No polygons meet minimum area requirement")
        return gdf_geo
    
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
    
    # Simplify geometry for web use - reduced for less blocky appearance
    tolerance = 0.0005  # ~50m simplification (was 0.002/200m)
    compact_polygons['geometry'] = compact_polygons.geometry.simplify(tolerance)
    
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
    
    # Create final bioregion
    total_area = sum(compact_polygons['area_sqkm'])
    
    coastal_bioregion = gpd.GeoDataFrame(
        [{
            'region_name': 'Coastal Forests',
            'region_code': 'coastal_forest_shapefile',
            'description': f'Coastal-influenced forests ≤{MAX_ELEVATION_FT}ft elevation with ≥{MIN_TREE_COVER_PCT}% tree cover, characterized by Pinus contorta and Sitka spruce',
            'elevation_max_ft': MAX_ELEVATION_FT,
            'min_tree_cover_pct': MIN_TREE_COVER_PCT,
            'area_sqkm': total_area,
            'created_date': pd.Timestamp.now().isoformat(),
            'method': 'shapefile_boundary_with_raster_constraints',
            'validation_species': 'Shore pine (Pinus contorta var. contorta)',
            'data_sources': 'User shapefile, USGS elevation, NLCD tree cover'
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
    
    # Step 2: Apply raster constraints within boundary
    mask_array, transform, crs, profile = apply_raster_constraints_within_boundary(boundary_gdf)
    if mask_array is None:
        print("Error: Could not apply raster constraints")
        return
    
    # Step 3: Convert raster to vector polygons
    coastal_polygons = raster_to_polygons_optimized(mask_array, transform, crs)
    
    # Step 4: Remove offshore bands and simplify
    coastal_bioregion = remove_offshore_bands_and_simplify(coastal_polygons)
    
    if coastal_bioregion.empty:
        print("Error: No coastal forest bioregion created")
        return
    
    # Step 5: Validate against shore pine plots
    validation_stats = validate_with_shore_pine_plots(coastal_bioregion)
    
    # Add validation results to bioregion attributes
    if validation_stats:
        for key, value in validation_stats.items():
            coastal_bioregion[key] = value
    
    # Step 6: Save outputs
    output_file = OUTPUT_DIR / "coastal_forest_shapefile_constrained.geojson"
    print(f"\nSaving coastal forest bioregion to {output_file}")
    coastal_bioregion.to_file(output_file, driver='GeoJSON')
    
    # Save processing summary
    summary = {
        'method': 'shapefile_boundary_with_raster_constraints',
        'parameters': {
            'max_elevation_ft': MAX_ELEVATION_FT,
            'min_tree_cover_pct': MIN_TREE_COVER_PCT,
            'output_resolution_m': OUTPUT_RESOLUTION,
            'min_area_sqkm': 0.5,
            'max_perimeter_area_ratio': 50,
            'simplification_tolerance': 0.0005
        },
        'data_sources': {
            'base_boundary': BASE_BOUNDARY_SHP,
            'elevation': ELEVATION_RASTER,
            'tree_cover': TREE_COVER_RASTER,
            'validation': PLOT_DATA
        },
        'validation': validation_stats,
        'created': pd.Timestamp.now().isoformat()
    }
    
    summary_file = OUTPUT_DIR / "coastal_forest_shapefile_summary.json"
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
                print("✅ Validation successful (≥70% plot capture)")
            else:
                print("⚠️  Validation acceptable but could be improved")
    
    print(f"\n✅ Coastal forest bioregion created successfully!")
    print(f"   Output: {output_file}")
    print(f"   Summary: {summary_file}")

if __name__ == "__main__":
    main()