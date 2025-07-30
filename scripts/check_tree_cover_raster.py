#!/usr/bin/env python3
"""
Check tree cover raster properties and sample coastal areas
"""

import rasterio
import numpy as np
from pathlib import Path

def check_raster(raster_path):
    """Check raster properties and sample data"""
    print(f"\n{'='*60}")
    print(f"Checking: {raster_path}")
    print('='*60)
    
    try:
        with rasterio.open(raster_path) as src:
            print(f"Shape: {src.shape} (height x width)")
            print(f"CRS: {src.crs}")
            print(f"Resolution: {src.res} (x, y)")
            print(f"Bounds: {src.bounds}")
            print(f"Data type: {src.dtypes[0]}")
            print(f"NoData value: {src.nodata}")
            
            # Sample different coastal areas
            samples = [
                ("Olympic Peninsula", -124.5, 47.5, -123.5, 48.5),
                ("WA Coast", -124.5, 46.0, -123.5, 47.0),
                ("OR North Coast", -124.5, 45.0, -123.5, 46.0),
                ("OR Central Coast", -124.5, 43.5, -123.5, 44.5),
                ("Port Angeles area", -124.0, 48.0, -123.0, 48.3)
            ]
            
            print("\nCoastal Area Samples:")
            print("-" * 50)
            
            for name, west, south, east, north in samples:
                try:
                    window = rasterio.windows.from_bounds(west, south, east, north, src.transform)
                    data = src.read(1, window=window)
                    
                    if data.size > 0:
                        # Filter out NoData values if specified
                        if src.nodata is not None:
                            valid_data = data[data != src.nodata]
                            if valid_data.size > 0:
                                print(f"{name:20} - Min: {valid_data.min():3d}, Max: {valid_data.max():3d}, Mean: {valid_data.mean():5.1f}%")
                            else:
                                print(f"{name:20} - All NoData values")
                        else:
                            print(f"{name:20} - Min: {data.min():3d}, Max: {data.max():3d}, Mean: {data.mean():5.1f}%")
                    else:
                        print(f"{name:20} - No data in bounds")
                        
                except Exception as e:
                    print(f"{name:20} - Error: {str(e)}")
                    
    except Exception as e:
        print(f"Error opening raster: {e}")

def main():
    """Check all available tree cover rasters"""
    print("Tree Cover Raster Analysis")
    print("==========================")
    
    # List of rasters to check
    rasters = [
        "outputs/mapbox_masks/pnw_tree_cover_90m_mapbox.tif",
        "outputs/mapbox_masks/pnw_tree_cover_30m_full.tif",
        "outputs/mapbox_masks/pnw_tree_cover_60m_half.tif",
        "outputs/mapbox_masks/pnw_tree_cover_120m_quarter.tif",
        "data/raw/nlcd_tcc_conus_2021_v2021-4.tif",
        "data/raw/nlcd_tcc_waor.tif"
    ]
    
    for raster_path in rasters:
        if Path(raster_path).exists():
            check_raster(raster_path)
        else:
            print(f"\nSkipping {raster_path} - file not found")
    
    print("\nAnalysis complete!")

if __name__ == "__main__":
    main()