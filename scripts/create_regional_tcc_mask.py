#!/usr/bin/env python3
"""
Regional Tree Canopy Cover (TCC) Mask Creator

Creates comprehensive tree canopy cover masks for the entire Pacific Northwest region
(Washington, Oregon, Idaho) that can be reused across multiple species distribution mapping projects.

Features:
- Full regional coverage (not species-specific)
- Incremental processing to fill gaps in existing masks
- Multiple threshold support (1%, 5%, 10%, 20%, etc.)
- High-resolution processing with configurable parameters
- Caching and resumable processing
- Geographic bounds validation and clipping
- Batch creation of commonly needed mask variants

Usage:
  python create_regional_tcc_mask.py --threshold 5 --resolution 30 --min-area 2.0
  python create_regional_tcc_mask.py --threshold 20 --resolution 60 --min-area 0.5
  python create_regional_tcc_mask.py --fill-gaps --existing-mask cache/species_masks/pnw_tcc_mask_5pct_30m_2.0sqkm.geojson
  python create_regional_tcc_mask.py --create-batch  # Create commonly needed variants
"""

import argparse
import sqlite3
import pandas as pd
import geopandas as gpd
import numpy as np
from pathlib import Path
from shapely.geometry import Point, box
from shapely.ops import unary_union
from tqdm import tqdm
import warnings
import json
from datetime import datetime
import time
warnings.filterwarnings('ignore')

# Try to import required libraries
try:
    import rasterio
    import rasterio.features
    from rasterio.mask import mask
    from rasterio.windows import from_bounds
    from rasterio.enums import Resampling
    from rasterio.warp import transform_bounds
    from rasterio.transform import Affine
    HAS_RASTERIO = True
except ImportError:
    HAS_RASTERIO = False
    print("❌ Rasterio not available - this script requires rasterio")
    exit(1)

# Try to import scipy for morphological operations
try:
    from scipy import ndimage
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    print("⚠️  Scipy not available - morphological processing will be disabled")

# Configuration
PACIFIC_NW_BOUNDS = [-125.0, 42.0, -110.0, 49.0]  # W, S, E, N (WGS84) - Full PNW region
US_CANADA_BORDER = 49.0  # Latitude to clip at northern border

# Tree canopy cover configuration
TCC_FILE = "data/raw/nlcd_tcc_conus_2021_v2021-4.tif"

# Cache directory
CACHE_DIR = Path("cache/species_masks")
CACHE_DIR.mkdir(parents=True, exist_ok=True)

# Processing defaults
DEFAULT_THRESHOLD_PCT = 5           # Default tree cover threshold
DEFAULT_RESOLUTION_M = 30           # Default processing resolution (meters)
DEFAULT_MIN_AREA_SQKM = 2.0        # Default minimum polygon area (sq km)
DEFAULT_MORPHOLOGICAL = True       # Default morphological processing
DEFAULT_HIGH_RES_MODE = True       # Default high-resolution mode

# Simplified versions to create for each mask
SIMPLIFIED_VERSIONS = {
    'qgis': {'tolerance': 0.01, 'desc': 'Ultra-simplified for QGIS display (~1km tolerance)'},
    'display': {'tolerance': 0.001, 'desc': 'For QGIS display (~100m tolerance)'},
    'analysis': {'tolerance': 0.0005, 'desc': 'For analysis work (~50m tolerance)'},
    'processing': {'tolerance': 0.0001, 'desc': 'For processing (higher detail, ~10m tolerance)'}
}

def get_cache_filename(threshold_pct, resolution_m, min_area_sqkm):
    """Generate standardized cache filename"""
    return CACHE_DIR / f"pnw_tcc_mask_{threshold_pct}pct_{resolution_m}m_{min_area_sqkm}sqkm.geojson"

def load_existing_mask(cache_file):
    """Load existing mask file if it exists"""
    if not cache_file.exists():
        return None
    
    print(f"📦 Loading existing mask: {cache_file.name}")
    file_size_mb = cache_file.stat().st_size / (1024 * 1024)
    
    # Skip very large files that might cause memory issues
    if file_size_mb > 500:
        print(f"⚠️  Cache file very large ({file_size_mb:.1f} MB) - consider manual review")
        response = input("Continue loading? (y/n): ")
        if response.lower() != 'y':
            return None
    
    try:
        gdf = gpd.read_file(cache_file)
        if len(gdf) == 0:
            print(f"⚠️  Cache file is empty - will rebuild")
            return None
        
        geometry = gdf.iloc[0].geometry
        if geometry is None or geometry.is_empty:
            print(f"⚠️  Cache file has invalid geometry - will rebuild")
            return None
        
        print(f"✅ Loaded existing mask ({file_size_mb:.1f} MB)")
        return geometry
        
    except Exception as e:
        print(f"❌ Error loading existing mask: {e}")
        return None

def save_mask(mask_geometry, cache_file, metadata):
    """Save mask to cache with metadata"""
    if mask_geometry is None:
        print("⚠️  No geometry to save")
        return
    
    print(f"💾 Saving mask to: {cache_file.name}")
    
    # Create GeoDataFrame with metadata
    gdf = gpd.GeoDataFrame([metadata], geometry=[mask_geometry], crs='EPSG:4326')
    
    try:
        gdf.to_file(cache_file, driver='GeoJSON')
        file_size_mb = cache_file.stat().st_size / (1024 * 1024)
        print(f"✅ Saved mask ({file_size_mb:.1f} MB)")
        
    except Exception as e:
        print(f"❌ Error saving mask: {e}")
        if cache_file.exists():
            cache_file.unlink()

def create_simplified_versions(mask_geometry, base_cache_file, metadata):
    """Create simplified versions of the mask for different use cases"""
    print(f"\n🎨 Creating simplified versions for different use cases...")
    
    from shapely.validation import make_valid
    
    # Count original vertices
    def count_vertices(geom):
        """Count vertices in a geometry"""
        if hasattr(geom, 'geoms'):
            return sum(len(g.exterior.coords) for g in geom.geoms if hasattr(g, 'exterior'))
        elif hasattr(geom, 'exterior'):
            return len(geom.exterior.coords)
        else:
            return 0
    
    original_vertices = count_vertices(mask_geometry)
    print(f"  📊 Original vertices: {original_vertices:,}")
    
    simplified_files = {}
    
    for version_name, version_config in SIMPLIFIED_VERSIONS.items():
        tolerance = version_config['tolerance']
        description = version_config['desc']
        
        print(f"\n  🔄 Creating {version_name} version...")
        print(f"     {description}")
        
        try:
            # Create simplified geometry
            simplified_geom = mask_geometry.simplify(tolerance, preserve_topology=True)
            
            # Validate geometry
            if not simplified_geom.is_valid:
                simplified_geom = make_valid(simplified_geom)
            
            # Count vertices in simplified version
            simplified_vertices = count_vertices(simplified_geom)
            reduction_pct = ((original_vertices - simplified_vertices) / original_vertices * 100) if original_vertices > 0 else 0
            
            print(f"     📉 Vertices: {simplified_vertices:,} ({reduction_pct:.1f}% reduction)")
            
            # Create filename for this version
            base_name = base_cache_file.stem  # Remove .geojson
            simplified_file = base_cache_file.parent / f"{base_name}_{version_name}.geojson"
            
            # Update metadata for this version
            simplified_metadata = metadata.copy()
            simplified_metadata.update({
                'version_type': version_name,
                'simplification_tolerance': tolerance,
                'original_vertices': original_vertices,
                'simplified_vertices': simplified_vertices,
                'vertex_reduction_pct': round(reduction_pct, 1),
                'description': description
            })
            
            # Save simplified version
            gdf = gpd.GeoDataFrame([simplified_metadata], geometry=[simplified_geom], crs='EPSG:4326')
            gdf.to_file(simplified_file, driver='GeoJSON')
            
            file_size_mb = simplified_file.stat().st_size / (1024 * 1024)
            print(f"     ✅ Saved: {simplified_file.name} ({file_size_mb:.1f} MB)")
            
            simplified_files[version_name] = simplified_file
            
        except Exception as e:
            print(f"     ❌ Error creating {version_name} version: {e}")
            continue
    
    print(f"\n📊 SIMPLIFIED VERSIONS SUMMARY:")
    print(f"  📁 Original: {base_cache_file.name} ({base_cache_file.stat().st_size / (1024 * 1024):.1f} MB)")
    for version_name, version_file in simplified_files.items():
        file_size_mb = version_file.stat().st_size / (1024 * 1024)
        reduction = ((base_cache_file.stat().st_size - version_file.stat().st_size) / base_cache_file.stat().st_size * 100)
        print(f"  📁 {version_name.capitalize()}: {version_file.name} ({file_size_mb:.1f} MB, {reduction:.1f}% smaller)")
    
    return simplified_files

def get_unprocessed_bounds(existing_mask_geometry, full_bounds):
    """Calculate bounds for areas not covered by existing mask"""
    if existing_mask_geometry is None:
        return [full_bounds]  # Process entire region
    
    print("🔍 Calculating unprocessed areas...")
    
    # Create full region polygon
    full_region = box(*full_bounds)
    
    # Find areas not covered by existing mask
    unprocessed_areas = full_region.difference(existing_mask_geometry)
    
    if unprocessed_areas.is_empty:
        print("✅ Entire region already processed")
        return []
    
    # Convert to list of bounds for processing
    bounds_list = []
    if hasattr(unprocessed_areas, 'geoms'):
        for geom in unprocessed_areas.geoms:
            bounds_list.append(list(geom.bounds))
    else:
        bounds_list.append(list(unprocessed_areas.bounds))
    
    total_unprocessed_area = sum([
        (bounds[2] - bounds[0]) * (bounds[3] - bounds[1]) 
        for bounds in bounds_list
    ])
    full_area = (full_bounds[2] - full_bounds[0]) * (full_bounds[3] - full_bounds[1])
    unprocessed_pct = (total_unprocessed_area / full_area) * 100
    
    print(f"📊 Unprocessed areas: {len(bounds_list)} regions ({unprocessed_pct:.1f}% of total)")
    
    return bounds_list

def create_regional_tcc_mask(
    threshold_pct=DEFAULT_THRESHOLD_PCT,
    resolution_m=DEFAULT_RESOLUTION_M,
    min_area_sqkm=DEFAULT_MIN_AREA_SQKM,
    morphological_processing=DEFAULT_MORPHOLOGICAL,
    high_res_mode=DEFAULT_HIGH_RES_MODE,
    fill_gaps=False,
    existing_mask_file=None,
    create_simplified=True
):
    """Create comprehensive TCC mask for Pacific Northwest region"""
    
    print(f"🌲 CREATING REGIONAL TREE CANOPY COVER MASK")
    print("="*60)
    print(f"Region: Pacific Northwest (WA, OR, ID)")
    print(f"Bounds: {PACIFIC_NW_BOUNDS}")
    print(f"Threshold: <{threshold_pct}% tree cover (areas to exclude)")
    print(f"Resolution: {resolution_m}m")
    print(f"Min area: {min_area_sqkm} sq km")
    print(f"Morphological processing: {'Enabled' if morphological_processing else 'Disabled'}")
    print(f"High-res mode: {'Enabled' if high_res_mode else 'Disabled'}")
    print()
    
    # Check if TCC file exists
    if not Path(TCC_FILE).exists():
        print(f"❌ TCC file not found: {TCC_FILE}")
        return None
    
    # Generate cache filename
    cache_file = get_cache_filename(threshold_pct, resolution_m, min_area_sqkm)
    
    # Load existing mask if filling gaps
    existing_mask = None
    if fill_gaps:
        if existing_mask_file:
            existing_mask = load_existing_mask(Path(existing_mask_file))
        else:
            existing_mask = load_existing_mask(cache_file)
    
    # Determine processing bounds
    if fill_gaps and existing_mask is not None:
        processing_bounds_list = get_unprocessed_bounds(existing_mask, PACIFIC_NW_BOUNDS)
        if not processing_bounds_list:
            print("✅ No additional processing needed")
            return existing_mask
    else:
        processing_bounds_list = [PACIFIC_NW_BOUNDS]
        print("🗺️  Processing entire Pacific Northwest region")
    
    # Process each bounds region
    all_low_cover_geometries = []
    if existing_mask is not None:
        all_low_cover_geometries.append(existing_mask)
    
    total_start_time = time.time()
    
    for i, bounds in enumerate(processing_bounds_list):
        print(f"\n🔄 Processing region {i+1}/{len(processing_bounds_list)}")
        print(f"📦 Bounds: {[round(x, 3) for x in bounds]}")
        
        region_start_time = time.time()
        
        try:
            with rasterio.open(TCC_FILE) as tcc_src:
                print(f"📍 TCC CRS: {tcc_src.crs}")
                print(f"📏 TCC original resolution: {tcc_src.res[0]:.0f}m")
                
                # Reproject bounds to TCC CRS
                bounds_tcc = transform_bounds("EPSG:4326", tcc_src.crs, *bounds)
                print(f"📦 Bounds in TCC CRS: {[round(x) for x in bounds_tcc]}")
                
                # Calculate downsampling
                original_res = tcc_src.res[0]
                if high_res_mode and resolution_m <= 30:
                    downsample_factor = max(1, int(resolution_m / original_res))
                    print(f"🔍 HIGH-RES MODE: {original_res:.0f}m → {original_res * downsample_factor:.0f}m")
                else:
                    downsample_factor = max(1, int(resolution_m / original_res))
                    print(f"📉 Standard processing: {original_res:.0f}m → {original_res * downsample_factor:.0f}m")
                
                # Read data
                window = from_bounds(*bounds_tcc, tcc_src.transform)
                window = window.intersection(
                    rasterio.windows.Window(0, 0, tcc_src.width, tcc_src.height)
                )
                
                if window.width <= 0 or window.height <= 0:
                    print(f"⚠️  No TCC data overlap for this region")
                    continue
                
                # Read with downsampling
                out_height = max(1, int(window.height // downsample_factor))
                out_width = max(1, int(window.width // downsample_factor))
                
                tcc_data = tcc_src.read(
                    1,
                    window=window,
                    out_shape=(out_height, out_width),
                    resampling=Resampling.average
                )
                
                # Update transform
                window_transform = tcc_src.window_transform(window)
                window_transform = Affine(
                    window_transform.a * downsample_factor,
                    window_transform.b,
                    window_transform.c,
                    window_transform.d,
                    window_transform.e * downsample_factor,
                    window_transform.f
                )
                
                data_mb = tcc_data.nbytes / (1024 * 1024)
                print(f"🗺️  Read TCC data: {tcc_data.shape} ({data_mb:.1f} MB)")
                
                # Create low tree cover mask
                valid_data_mask = (tcc_data >= 0) & (tcc_data <= 100)
                low_tree_cover_mask = (tcc_data < threshold_pct) & valid_data_mask
                
                low_cover_pixels = np.sum(low_tree_cover_mask)
                total_pixels = np.sum(valid_data_mask)
                low_cover_pct = (low_cover_pixels / total_pixels * 100) if total_pixels > 0 else 0
                
                print(f"📊 Low tree cover: {low_cover_pixels:,}/{total_pixels:,} pixels ({low_cover_pct:.1f}%)")
                
                if low_cover_pixels == 0:
                    print(f"✅ All areas have ≥{threshold_pct}% tree cover")
                    continue
                
                # Apply morphological processing
                if HAS_SCIPY and morphological_processing:
                    print(f"🔄 Applying morphological processing...")
                    structure = np.ones((5, 5))
                    
                    if high_res_mode:
                        # Aggressive processing for high-res
                        low_tree_cover_mask = ndimage.binary_closing(
                            low_tree_cover_mask, structure=structure, iterations=2
                        )
                        low_tree_cover_mask = ndimage.binary_opening(
                            low_tree_cover_mask, structure=structure, iterations=2
                        )
                        low_tree_cover_mask = ndimage.binary_closing(
                            low_tree_cover_mask, structure=structure, iterations=1
                        )
                        print(f"✅ High-res morphological processing complete")
                    else:
                        # Standard processing
                        low_tree_cover_mask = ndimage.binary_closing(
                            low_tree_cover_mask, structure=structure, iterations=3
                        )
                        low_tree_cover_mask = ndimage.binary_opening(
                            low_tree_cover_mask, structure=structure, iterations=2
                        )
                        print(f"✅ Standard morphological processing complete")
                
                # Convert to polygons with size filtering
                print(f"🔄 Converting to polygons...")
                
                pixel_area_sqm = abs(window_transform.a * window_transform.e)
                min_pixels = int((min_area_sqkm * 1_000_000) / pixel_area_sqm)
                
                print(f"📏 Min area: {min_area_sqkm} sq km = {min_pixels} pixels")
                
                region_shapes = []
                shape_count = 0
                processed_count = 0
                
                start_time = time.time()
                last_progress = start_time
                
                for geom, value in rasterio.features.shapes(
                    low_tree_cover_mask.astype(np.uint8),
                    mask=low_tree_cover_mask,
                    transform=window_transform
                ):
                    processed_count += 1
                    
                    # Progress reporting
                    current_time = time.time()
                    if current_time - last_progress >= 10.0:  # Every 10 seconds
                        elapsed = current_time - start_time
                        print(f"    ⏱️  Progress: {processed_count:,} processed, {shape_count:,} kept ({elapsed:.1f}s)")
                        last_progress = current_time
                    
                    if value == 1:
                        # Size filtering
                        from shapely.geometry import shape as shapely_shape
                        try:
                            polygon = shapely_shape(geom)
                            area_sqkm = polygon.area / 1_000_000
                            
                            if area_sqkm >= min_area_sqkm:
                                region_shapes.append(geom)
                                shape_count += 1
                        except:
                            continue
                    
                    # Safety limit
                    if shape_count >= 100000:
                        print(f"⚠️  Reached safety limit of 100,000 polygons")
                        break
                
                total_elapsed = time.time() - start_time
                print(f"✅ Polygon extraction: {processed_count:,} processed, {shape_count:,} kept ({total_elapsed:.1f}s)")
                
                if region_shapes:
                    # Create geometries
                    from shapely.geometry import shape
                    geometries = []
                    
                    print(f"🔄 Creating {len(region_shapes):,} geometry objects...")
                    for geom_dict in tqdm(region_shapes, desc="Creating geometries"):
                        try:
                            geom = shape(geom_dict)
                            if geom.is_valid and not geom.is_empty:
                                geometries.append(geom)
                        except:
                            continue
                    
                    if geometries:
                        # Create GeoDataFrame and reproject
                        region_gdf = gpd.GeoDataFrame(
                            [{'source_region': i+1}] * len(geometries),
                            geometry=geometries,
                            crs=tcc_src.crs
                        )
                        
                        if region_gdf.crs != "EPSG:4326":
                            print(f"🔄 Reprojecting to WGS84...")
                            region_gdf = region_gdf.to_crs("EPSG:4326")
                        
                        # Merge region polygons
                        print(f"🔗 Merging {len(region_gdf):,} polygons...")
                        try:
                            region_merged = unary_union(region_gdf.geometry.tolist())
                            all_low_cover_geometries.append(region_merged)
                            print(f"✅ Region {i+1} processed successfully")
                        except Exception as e:
                            print(f"❌ Error merging region {i+1}: {e}")
                            continue
                
                region_elapsed = time.time() - region_start_time
                print(f"⏱️  Region {i+1} completed in {region_elapsed:.1f}s")
                
        except Exception as e:
            print(f"❌ Error processing region {i+1}: {e}")
            continue
    
    # Combine all regions
    if not all_low_cover_geometries:
        print("❌ No low tree cover areas found in any region")
        return None
    
    print(f"\n🔗 FINAL MERGE: Combining {len(all_low_cover_geometries)} region(s)...")
    final_start_time = time.time()
    
    try:
        if len(all_low_cover_geometries) == 1:
            final_geometry = all_low_cover_geometries[0]
        else:
            final_geometry = unary_union(all_low_cover_geometries)
        
        final_elapsed = time.time() - final_start_time
        print(f"✅ Final merge completed in {final_elapsed:.1f}s")
        
        # Clip to Pacific Northwest bounds
        print(f"✂️  Clipping to Pacific Northwest bounds...")
        pnw_box = box(*PACIFIC_NW_BOUNDS)
        final_geometry = final_geometry.intersection(pnw_box)
        
        # Save result
        metadata = {
            'mask_type': 'regional_tcc',
            'threshold_pct': threshold_pct,
            'resolution_m': resolution_m,
            'min_area_sqkm': min_area_sqkm,
            'region': 'Pacific Northwest (WA, OR, ID)',
            'bounds': PACIFIC_NW_BOUNDS,
            'created_date': datetime.now().isoformat(),
            'processing_time_seconds': time.time() - total_start_time,
            'morphological_processing': morphological_processing,
            'high_res_mode': high_res_mode,
            'fill_gaps_mode': fill_gaps
        }
        
        save_mask(final_geometry, cache_file, metadata)
        
        # Create simplified versions for different use cases (unless disabled)
        simplified_files = {}
        if create_simplified and not fill_gaps:  # Don't create simplified versions when filling gaps
            simplified_files = create_simplified_versions(final_geometry, cache_file, metadata)
        
        total_elapsed = time.time() - total_start_time
        print(f"\n✅ REGIONAL TCC MASK COMPLETE!")
        print(f"⏱️  Total processing time: {total_elapsed:.1f}s ({total_elapsed/60:.1f} min)")
        print(f"📁 Primary output: {cache_file}")
        print(f"📁 Simplified versions: {len(simplified_files)} files created")
        
        print(f"\n💡 USAGE RECOMMENDATIONS:")
        print(f"  🖥️  For QGIS display: {cache_file.parent / (cache_file.stem + '_qgis.geojson')}")
        print(f"  📊 For analysis: {cache_file.parent / (cache_file.stem + '_analysis.geojson')}")
        print(f"  ⚙️  For processing: {cache_file.parent / (cache_file.stem + '_processing.geojson')}")
        
        return final_geometry
        
    except Exception as e:
        print(f"❌ Error in final merge: {e}")
        return None

def list_existing_masks():
    """List existing TCC mask files with their versions"""
    print("📋 EXISTING TCC MASK FILES")
    print("="*50)
    
    # Group masks by base configuration
    mask_groups = {}
    mask_files = list(CACHE_DIR.glob("pnw_tcc_mask_*.geojson"))
    
    for mask_file in mask_files:
        try:
            # Parse filename: pnw_tcc_mask_5pct_30m_2.0sqkm.geojson or pnw_tcc_mask_5pct_30m_2.0sqkm_processing.geojson
            name_parts = mask_file.stem.split('_')
            
            if len(name_parts) >= 6:
                threshold = name_parts[3].replace('pct', '')
                resolution = name_parts[4].replace('m', '')
                min_area = name_parts[5].replace('sqkm', '')
                
                # Check if this is a simplified version
                version_type = 'original'
                if len(name_parts) > 6:
                    version_type = name_parts[6]
                
                base_config = f"{threshold}pct_{resolution}m_{min_area}sqkm"
                
                if base_config not in mask_groups:
                    mask_groups[base_config] = {}
                
                file_size_mb = mask_file.stat().st_size / (1024 * 1024)
                mask_groups[base_config][version_type] = {
                    'file': mask_file,
                    'size_mb': file_size_mb
                }
        except:
            continue  # Skip files with unexpected naming
    
    if not mask_groups:
        print("  (No TCC mask files found)")
        return
    
    for config, versions in sorted(mask_groups.items()):
        print(f"\n📁 Configuration: {config}")
        
        # Show versions in order of detail level
        version_order = ['original', 'processing', 'analysis', 'display', 'qgis']
        for version_type in version_order:
            if version_type in versions:
                version_info = versions[version_type]
                print(f"  └─ {version_type}: {version_info['file'].name} ({version_info['size_mb']:.1f} MB)")
    
    print(f"\n💡 USAGE TIPS:")
    print(f"  🖥️  Use 'qgis' version for QGIS display (smallest, fastest)")
    print(f"  📊 Use 'analysis' version for analysis work (balanced)")
    print(f"  ⚙️  Use 'processing' version for species processing (detailed)")
    print(f"  📄 Use 'original' version for maximum detail (largest)")

def get_processing_version_path(threshold_pct, resolution_m, min_area_sqkm):
    """Get the processing version path for a TCC mask configuration"""
    base_filename = f"pnw_tcc_mask_{threshold_pct}pct_{resolution_m}m_{min_area_sqkm}sqkm"
    processing_file = CACHE_DIR / f"{base_filename}_processing.geojson"
    return processing_file

def main():
    parser = argparse.ArgumentParser(description='Create regional tree canopy cover masks with simplified versions')
    parser.add_argument('--threshold', type=int, default=DEFAULT_THRESHOLD_PCT,
                       help=f'Tree cover threshold percentage (default: {DEFAULT_THRESHOLD_PCT})')
    parser.add_argument('--resolution', type=int, default=DEFAULT_RESOLUTION_M,
                       help=f'Processing resolution in meters (default: {DEFAULT_RESOLUTION_M})')
    parser.add_argument('--min-area', type=float, default=DEFAULT_MIN_AREA_SQKM,
                       help=f'Minimum polygon area in sq km (default: {DEFAULT_MIN_AREA_SQKM})')
    parser.add_argument('--no-morphological', action='store_true',
                       help='Disable morphological processing')
    parser.add_argument('--no-high-res', action='store_true',
                       help='Disable high-resolution processing mode')
    parser.add_argument('--fill-gaps', action='store_true',
                       help='Fill gaps in existing mask instead of recreating')
    parser.add_argument('--existing-mask', type=str,
                       help='Path to existing mask file for gap filling')
    parser.add_argument('--list-existing', action='store_true',
                       help='List existing mask files with versions')
    parser.add_argument('--no-simplified-versions', action='store_true',
                       help='Skip creating simplified versions (qgis, display, analysis, processing)')
    
    args = parser.parse_args()
    
    if args.list_existing:
        list_existing_masks()
        return
    
    # Validate TCC file
    if not Path(TCC_FILE).exists():
        print(f"❌ TCC file not found: {TCC_FILE}")
        print("Please ensure the NLCD tree canopy cover file is available")
        return
    
    # Create mask
    result = create_regional_tcc_mask(
        threshold_pct=args.threshold,
        resolution_m=args.resolution,
        min_area_sqkm=args.min_area,
        morphological_processing=not args.no_morphological,
        high_res_mode=not args.no_high_res,
        fill_gaps=args.fill_gaps,
        existing_mask_file=args.existing_mask,
        create_simplified=not args.no_simplified_versions
    )
    
    if result is not None:
        print("\n🎉 SUCCESS! Regional TCC mask created.")
        print("\n📖 Usage in other scripts:")
        cache_file = get_cache_filename(args.threshold, args.resolution, args.min_area)
        print(f"  TCC_MASK_FILE = '{cache_file}'")
        print(f"  # Load with: gdf = gpd.read_file(TCC_MASK_FILE)")
        print(f"  # Use geometry: mask = gdf.geometry.iloc[0]")
        
        print(f"\n💡 Species distribution script usage:")
        print(f"  python create_species_distribution_maps_cached.py --regional-tcc-mask {cache_file}")
    else:
        print("❌ Failed to create regional TCC mask")

if __name__ == "__main__":
    main() 