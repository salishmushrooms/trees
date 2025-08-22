#!/usr/bin/env python3
"""
Create high-resolution PRISM precipitation raster for Pacific Northwest states.

This script creates a seamless, high-resolution precipitation raster covering
Oregon, Washington, Idaho, and Northern California for use in zonal statistics
calculations. The output is optimized for efficient analysis with proper
coordinate system, resolution, and data format.

Key features:
- Native PRISM resolution (800m) preserved
- Covers OR, WA, ID, and Northern CA (north of 40°N)
- Albers Equal Area projection for accurate area calculations
- Optimized for zonal statistics with proper NoData handling
- Includes both annual and monthly precipitation data

Usage:
    python scripts/create_prism_precipitation_raster.py [--monthly]
"""

import os
import sys
import numpy as np
import rasterio
from rasterio.merge import merge
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.crs import CRS
import geopandas as gpd
from shapely.geometry import box
import json
from datetime import datetime
import argparse

# Define the geographic extent for OR, WA, ID, and Northern CA
STUDY_EXTENT = {
    'west': -125.0,  # Western edge of OR/WA coast
    'east': -110.0,  # Eastern edge of Idaho
    'south': 40.0,   # Northern California boundary
    'north': 49.5    # Canadian border
}

# Target CRS for analysis - Albers Equal Area for accurate area calculations
TARGET_CRS = CRS.from_epsg(5070)  # NAD83 / Conus Albers

def get_state_boundaries():
    """Load state boundaries for extent definition."""
    # Try to load from common locations
    state_paths = [
        "data/raw/cb_2018_us_state_20m/cb_2018_us_state_20m.shp",
        "data/raw/states/cb_2018_us_state_20m.shp",
        "data/raw/states.shp"
    ]
    
    states_gdf = None
    for path in state_paths:
        if os.path.exists(path):
            states_gdf = gpd.read_file(path)
            break
    
    if states_gdf is None:
        print("Warning: Could not find state boundaries shapefile.")
        print("Using bounding box extent instead.")
        return None
    
    # Filter to our states of interest
    target_states = ['Oregon', 'Washington', 'Idaho', 'California']
    states_gdf = states_gdf[states_gdf['NAME'].isin(target_states)]
    
    # For California, clip to northern portion
    ca_mask = states_gdf['NAME'] == 'California'
    if ca_mask.any():
        ca_geom = states_gdf.loc[ca_mask, 'geometry'].iloc[0]
        north_ca_box = box(STUDY_EXTENT['west'], STUDY_EXTENT['south'], 
                          STUDY_EXTENT['east'], 42.0)  # CA northern border
        ca_clipped = ca_geom.intersection(north_ca_box)
        states_gdf.loc[ca_mask, 'geometry'] = ca_clipped
    
    return states_gdf

def create_study_area_mask(bounds, crs):
    """Create a mask for the study area."""
    # Create bounding box
    study_box = box(bounds['west'], bounds['south'], 
                    bounds['east'], bounds['north'])
    
    # Create GeoDataFrame
    study_gdf = gpd.GeoDataFrame([1], geometry=[study_box], crs='EPSG:4326')
    
    # Load state boundaries if available
    states_gdf = get_state_boundaries()
    if states_gdf is not None:
        # Reproject to WGS84 for intersection
        states_gdf = states_gdf.to_crs('EPSG:4326')
        # Get union of state geometries
        study_area = states_gdf.unary_union
        study_gdf = gpd.GeoDataFrame([1], geometry=[study_area], crs='EPSG:4326')
    
    # Reproject to target CRS
    study_gdf = study_gdf.to_crs(crs)
    
    return study_gdf

def load_prism_data(data_type='annual'):
    """Load PRISM precipitation data files."""
    base_path = "data/raw/PRISM_ppt_30yr_normal_800mM4_all_bil"
    
    if data_type == 'annual':
        file_path = os.path.join(base_path, "PRISM_ppt_30yr_normal_800mM4_annual_bil.bil")
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"PRISM annual precipitation file not found: {file_path}")
        return [file_path]
    
    elif data_type == 'monthly':
        # Load all monthly files
        monthly_files = []
        months = ['01', '02', '03', '04', '05', '06', 
                 '07', '08', '09', '10', '11', '12']
        
        for month in months:
            file_path = os.path.join(base_path, f"PRISM_ppt_30yr_normal_800mM4_{month}_bil.bil")
            if os.path.exists(file_path):
                monthly_files.append(file_path)
            else:
                print(f"Warning: Missing month {month} file")
        
        if not monthly_files:
            raise FileNotFoundError("No PRISM monthly precipitation files found")
        
        return monthly_files

def process_precipitation_raster(input_files, output_path, study_area_gdf):
    """Process and merge PRISM precipitation data."""
    
    # For single file (annual), no merging needed
    if len(input_files) == 1:
        with rasterio.open(input_files[0]) as src:
            # Get metadata
            src_crs = src.crs
            src_transform = src.transform
            src_bounds = src.bounds
            
            print(f"Source PRISM data:")
            print(f"  CRS: {src_crs}")
            print(f"  Resolution: {src.res[0]:.1f}m x {src.res[1]:.1f}m")
            print(f"  Bounds: {src_bounds}")
            
            # Calculate transform for target CRS and extent
            study_bounds = study_area_gdf.total_bounds
            dst_transform, dst_width, dst_height = calculate_default_transform(
                src_crs, TARGET_CRS, src.width, src.height,
                left=study_bounds[0], bottom=study_bounds[1],
                right=study_bounds[2], top=study_bounds[3],
                resolution=800  # Maintain 800m resolution
            )
            
            # Create output array
            output_data = np.zeros((dst_height, dst_width), dtype=np.float32)
            
            # Reproject
            reproject(
                source=rasterio.band(src, 1),
                destination=output_data,
                src_transform=src.transform,
                src_crs=src.crs,
                dst_transform=dst_transform,
                dst_crs=TARGET_CRS,
                resampling=Resampling.bilinear,
                dst_nodata=-9999
            )
            
            # Write output
            with rasterio.open(
                output_path,
                'w',
                driver='GTiff',
                height=dst_height,
                width=dst_width,
                count=1,
                dtype='float32',
                crs=TARGET_CRS,
                transform=dst_transform,
                nodata=-9999,
                compress='lzw',
                tiled=True,
                blockxsize=512,
                blockysize=512
            ) as dst:
                dst.write(output_data, 1)
                
                # Add metadata
                dst.update_tags(
                    source='PRISM Climate Group',
                    description='Annual precipitation normals (1981-2010)',
                    units='millimeters',
                    coverage='OR, WA, ID, Northern CA',
                    processing_date=datetime.now().isoformat(),
                    resolution='800m'
                )
    
    else:
        # For monthly data, create multi-band raster
        # First, get dimensions from first file
        with rasterio.open(input_files[0]) as src:
            src_crs = src.crs
            study_bounds = study_area_gdf.total_bounds
            
            dst_transform, dst_width, dst_height = calculate_default_transform(
                src_crs, TARGET_CRS, src.width, src.height,
                left=study_bounds[0], bottom=study_bounds[1],
                right=study_bounds[2], top=study_bounds[3],
                resolution=800
            )
        
        # Create output with 12 bands
        with rasterio.open(
            output_path,
            'w',
            driver='GTiff',
            height=dst_height,
            width=dst_width,
            count=12,
            dtype='float32',
            crs=TARGET_CRS,
            transform=dst_transform,
            nodata=-9999,
            compress='lzw',
            tiled=True,
            blockxsize=512,
            blockysize=512
        ) as dst:
            
            # Process each month
            for i, file_path in enumerate(input_files, 1):
                print(f"  Processing month {i}...")
                
                with rasterio.open(file_path) as src:
                    # Create output array for this month
                    output_data = np.zeros((dst_height, dst_width), dtype=np.float32)
                    
                    # Reproject
                    reproject(
                        source=rasterio.band(src, 1),
                        destination=output_data,
                        src_transform=src.transform,
                        src_crs=src.crs,
                        dst_transform=dst_transform,
                        dst_crs=TARGET_CRS,
                        resampling=Resampling.bilinear,
                        dst_nodata=-9999
                    )
                    
                    # Write band
                    dst.write(output_data, i)
                    dst.set_band_description(i, f"Month_{i:02d}_precipitation")
            
            # Add metadata
            dst.update_tags(
                source='PRISM Climate Group',
                description='Monthly precipitation normals (1981-2010)',
                units='millimeters',
                coverage='OR, WA, ID, Northern CA',
                processing_date=datetime.now().isoformat(),
                resolution='800m',
                band_info='Bands 1-12 represent January-December'
            )
    
    print(f"\nOutput raster saved: {output_path}")
    print(f"  CRS: {TARGET_CRS}")
    print(f"  Dimensions: {dst_width} x {dst_height}")
    print(f"  Resolution: 800m x 800m")
    print(f"  Bands: {12 if len(input_files) > 1 else 1}")

def create_summary_statistics(raster_path, output_path):
    """Generate summary statistics for the precipitation raster."""
    
    with rasterio.open(raster_path) as src:
        # Initialize statistics
        stats = {
            'metadata': {
                'source': 'PRISM Climate Group - 30-year normals (1981-2010)',
                'processing_date': datetime.now().isoformat(),
                'coverage': 'Oregon, Washington, Idaho, Northern California',
                'crs': str(src.crs),
                'resolution_m': src.res[0],
                'dimensions': f"{src.width} x {src.height}",
                'bands': src.count,
                'units': 'millimeters',
                'nodata_value': src.nodata
            },
            'extent': {
                'west': src.bounds.left,
                'east': src.bounds.right,
                'south': src.bounds.bottom,
                'north': src.bounds.top
            }
        }
        
        if src.count == 1:
            # Annual data
            data = src.read(1)
            valid_data = data[data != src.nodata]
            
            stats['annual_statistics'] = {
                'min_mm': float(np.min(valid_data)),
                'max_mm': float(np.max(valid_data)),
                'mean_mm': float(np.mean(valid_data)),
                'median_mm': float(np.median(valid_data)),
                'std_mm': float(np.std(valid_data)),
                'min_inches': float(np.min(valid_data) / 25.4),
                'max_inches': float(np.max(valid_data) / 25.4),
                'mean_inches': float(np.mean(valid_data) / 25.4),
                'valid_pixels': int(len(valid_data)),
                'total_pixels': int(data.size),
                'coverage_percent': float(len(valid_data) / data.size * 100)
            }
        
        else:
            # Monthly data
            stats['monthly_statistics'] = {}
            month_names = ['January', 'February', 'March', 'April', 'May', 'June',
                          'July', 'August', 'September', 'October', 'November', 'December']
            
            for i in range(1, src.count + 1):
                data = src.read(i)
                valid_data = data[data != src.nodata]
                
                stats['monthly_statistics'][month_names[i-1]] = {
                    'band': i,
                    'min_mm': float(np.min(valid_data)),
                    'max_mm': float(np.max(valid_data)),
                    'mean_mm': float(np.mean(valid_data)),
                    'median_mm': float(np.median(valid_data)),
                    'mean_inches': float(np.mean(valid_data) / 25.4)
                }
        
        # Add usage notes
        stats['usage_notes'] = {
            'zonal_statistics': 'Use rasterio or rasterstats for efficient zonal calculations',
            'coordinate_system': 'Albers Equal Area (EPSG:5070) for accurate area calculations',
            'resolution': 'Native PRISM 800m resolution preserved',
            'optimization': 'Tiled GeoTIFF with LZW compression for fast access'
        }
    
    # Save statistics
    with open(output_path, 'w') as f:
        json.dump(stats, f, indent=2)
    
    print(f"\nSummary statistics saved: {output_path}")
    return stats

def main():
    """Main processing function."""
    parser = argparse.ArgumentParser(description='Create PRISM precipitation raster for PNW states')
    parser.add_argument('--monthly', action='store_true', 
                       help='Create monthly precipitation raster (12 bands) instead of annual')
    args = parser.parse_args()
    
    print("Creating high-resolution PRISM precipitation raster...")
    print(f"Coverage: Oregon, Washington, Idaho, Northern California")
    print(f"Data type: {'Monthly (12 bands)' if args.monthly else 'Annual'}")
    print("=" * 70)
    
    # Create output directory
    os.makedirs("outputs/climate/zonal_stats", exist_ok=True)
    
    # Step 1: Define study area
    print("\n1. Defining study area...")
    study_area_gdf = create_study_area_mask(STUDY_EXTENT, TARGET_CRS)
    print(f"   Study area bounds: {study_area_gdf.total_bounds}")
    
    # Step 2: Load PRISM data
    print("\n2. Loading PRISM precipitation data...")
    data_type = 'monthly' if args.monthly else 'annual'
    input_files = load_prism_data(data_type)
    print(f"   Found {len(input_files)} input file(s)")
    
    # Step 3: Process and reproject data
    print("\n3. Processing precipitation raster...")
    output_filename = f"prism_precipitation_pnw_{'monthly' if args.monthly else 'annual'}_800m.tif"
    output_path = os.path.join("outputs/climate/zonal_stats", output_filename)
    
    process_precipitation_raster(input_files, output_path, study_area_gdf)
    
    # Step 4: Generate summary statistics
    print("\n4. Generating summary statistics...")
    stats_filename = f"prism_precipitation_pnw_{'monthly' if args.monthly else 'annual'}_stats.json"
    stats_path = os.path.join("outputs/climate/zonal_stats", stats_filename)
    
    stats = create_summary_statistics(output_path, stats_path)
    
    # Print summary
    print("\n" + "=" * 70)
    print("Processing completed successfully!")
    print(f"\nOutput files:")
    print(f"  Raster: {output_path}")
    print(f"  Statistics: {stats_path}")
    
    if 'annual_statistics' in stats:
        print(f"\nPrecipitation range:")
        print(f"  {stats['annual_statistics']['min_inches']:.1f} - {stats['annual_statistics']['max_inches']:.1f} inches")
        print(f"  Mean: {stats['annual_statistics']['mean_inches']:.1f} inches")
    
    print(f"\nOptimized for zonal statistics with:")
    print(f"  - Albers Equal Area projection (EPSG:5070)")
    print(f"  - 800m native resolution")
    print(f"  - Tiled GeoTIFF format")
    print(f"  - Coverage: {stats['annual_statistics']['coverage_percent']:.1f}% of region")

if __name__ == "__main__":
    main()