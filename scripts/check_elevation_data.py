#!/usr/bin/env python3
"""
Check for elevation/DEM data in the project
"""

import os
from pathlib import Path
import rasterio

def find_elevation_files():
    """Find all potential elevation data files"""
    print("Searching for elevation/DEM data...")
    print("="*60)
    
    # Search patterns
    patterns = ['*elevation*', '*elev*', '*dem*', '*height*', '*terrain*']
    extensions = ['.tif', '.tiff', '.img', '.bil', '.hgt']
    
    # Directories to search
    search_dirs = [
        'data/raw',
        'cache',
        'outputs',
        '.'
    ]
    
    found_files = set()
    
    for search_dir in search_dirs:
        if not Path(search_dir).exists():
            continue
            
        print(f"\nSearching in {search_dir}...")
        
        for pattern in patterns:
            # Search for files matching pattern
            for ext in extensions:
                try:
                    for file in Path(search_dir).rglob(f"*{pattern}*{ext}"):
                        found_files.add(file)
                except ValueError:
                    # Skip invalid patterns
                    continue
    
    # Also check specific mapbox masks directory
    mapbox_dir = Path("outputs/mapbox_masks")
    if mapbox_dir.exists():
        print(f"\nChecking {mapbox_dir}...")
        for file in mapbox_dir.glob("*elevation*"):
            found_files.add(file)
            
    # Display results
    print("\n" + "="*60)
    print("FOUND ELEVATION/DEM FILES:")
    print("="*60)
    
    if found_files:
        for file in sorted(found_files):
            print(f"\n{file}")
            # Get file info
            size_mb = file.stat().st_size / 1024 / 1024
            print(f"  Size: {size_mb:.1f} MB")
            
            # If it's a raster, try to get properties
            if file.suffix.lower() in ['.tif', '.tiff']:
                try:
                    with rasterio.open(file) as src:
                        print(f"  Shape: {src.shape}")
                        print(f"  CRS: {src.crs}")
                        print(f"  Resolution: {src.res}")
                except:
                    print("  (Unable to read raster properties)")
    else:
        print("No elevation/DEM files found")
    
    # Check for elevation masks
    print("\n" + "="*60)
    print("ELEVATION MASKS:")
    print("="*60)
    
    mask_dir = Path("cache/species_masks")
    if mask_dir.exists():
        elevation_masks = sorted(mask_dir.glob("*elevation*"))
        if elevation_masks:
            for mask in elevation_masks:
                print(f"  {mask.name}")
        else:
            print("  No elevation masks found")
    
    # Look for DEM in unexpected places
    print("\n" + "="*60)
    print("CHECKING FOR DEM IN NLCD FILES:")
    print("="*60)
    
    nlcd_files = [
        "data/raw/nlcd_tcc_conus_2021_v2021-4.tif",
        "data/raw/nlcd_tcc_waor.tif"
    ]
    
    for nlcd_file in nlcd_files:
        if Path(nlcd_file).exists():
            print(f"\nChecking {nlcd_file}...")
            try:
                with rasterio.open(nlcd_file) as src:
                    print(f"  Bands: {src.count}")
                    print(f"  Band descriptions: {src.descriptions}")
                    if src.count > 1:
                        print("  Note: Multiple bands detected - might include elevation")
            except:
                print("  Unable to read file")

if __name__ == "__main__":
    find_elevation_files()