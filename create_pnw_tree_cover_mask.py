#!/usr/bin/env python3
"""
Create optimized PNW Tree Cover Mask for Mapbox
Extracts tree cover data from NLCD for Oregon, Washington, Idaho
Optimizes resolution and compression to stay under 300MB
"""

import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.windows import from_bounds
import numpy as np
import os
from pathlib import Path

def create_pnw_tree_cover_mask():
    """Create optimized tree cover mask for PNW region"""
    
    # Input and output paths
    input_tif = "data/raw/nlcd_tcc_conus_2021_v2021-4.tif"
    output_dir = Path("outputs/mapbox_masks")
    output_dir.mkdir(exist_ok=True)
    
    # PNW bounding box (Oregon, Washington, Idaho)
    # Expanded slightly to ensure full coverage
    pnw_bounds = {
        'west': -125.0,   # Oregon coast
        'east': -110.0,   # Eastern Idaho  
        'south': 41.5,    # Southern Oregon
        'north': 49.5     # Northern Washington
    }
    
    print("🌲 Creating PNW Tree Cover Mask for Mapbox")
    print(f"📁 Input: {input_tif}")
    print(f"🎯 Target region: {pnw_bounds}")
    
    with rasterio.open(input_tif) as src:
        print(f"📊 Source resolution: {src.res[0]}m")
        print(f"📐 Source CRS: {src.crs}")
        print(f"📏 Source size: {src.width} x {src.height}")
        
        # Transform bounds to source CRS
        from rasterio.warp import transform_bounds
        pnw_bounds_src = transform_bounds(
            'EPSG:4326', src.crs,
            pnw_bounds['west'], pnw_bounds['south'],
            pnw_bounds['east'], pnw_bounds['north']
        )
        
        print(f"🗺️  PNW bounds in source CRS: {pnw_bounds_src}")
        
        # Create window for PNW region
        window = from_bounds(*pnw_bounds_src, src.transform)
        print(f"🔲 Extraction window: {window}")
        
        # Read the PNW subset
        print("📖 Reading PNW subset...")
        pnw_data = src.read(1, window=window)
        pnw_transform = src.window_transform(window)
        
        print(f"✅ Extracted region: {pnw_data.shape}")
        print(f"📊 Data range: {pnw_data.min()} - {pnw_data.max()}")
        print(f"🌳 Tree cover stats: mean={pnw_data[pnw_data < 255].mean():.1f}%")
        
        # Calculate different resolution options
        resolutions = [
            (30, "30m_full"),      # Original resolution
            (60, "60m_half"),      # Half resolution  
            (90, "90m_third"),     # Third resolution
            (120, "120m_quarter")  # Quarter resolution
        ]
        
        for target_res, suffix in resolutions:
            output_file = output_dir / f"pnw_tree_cover_{suffix}.tif"
            
            # Calculate new dimensions
            scale_factor = target_res / 30.0  # Original is 30m
            new_height = int(pnw_data.shape[0] / scale_factor)
            new_width = int(pnw_data.shape[1] / scale_factor)
            
            print(f"\n🔄 Creating {target_res}m resolution version...")
            print(f"📐 New size: {new_width} x {new_height}")
            
            # Calculate new transform
            new_transform = rasterio.transform.from_bounds(
                *pnw_bounds_src, new_width, new_height
            )
            
            # Resample data
            resampled_data = np.empty((new_height, new_width), dtype=np.uint8)
            
            reproject(
                pnw_data,
                resampled_data,
                src_transform=pnw_transform,
                src_crs=src.crs,
                dst_transform=new_transform,
                dst_crs=src.crs,
                resampling=Resampling.average
            )
            
            # Optimize for Mapbox:
            # - Convert no-data (255) to 0 for transparency
            # - Keep tree cover values 1-100
            mapbox_data = resampled_data.copy()
            mapbox_data[resampled_data == 255] = 0  # No data = transparent
            mapbox_data[resampled_data == 0] = 1    # 0% cover = very light
            
            # Write optimized TIF
            profile = {
                'driver': 'GTiff',
                'dtype': 'uint8',
                'nodata': 0,
                'width': new_width,
                'height': new_height,
                'count': 1,
                'crs': 'EPSG:4326',  # Web Mercator for Mapbox
                'transform': rasterio.transform.from_bounds(
                    pnw_bounds['west'], pnw_bounds['south'],
                    pnw_bounds['east'], pnw_bounds['north'],
                    new_width, new_height
                ),
                'compress': 'lzw',
                'tiled': True,
                'blockxsize': 512,
                'blockysize': 512,
                'BIGTIFF': 'IF_NEEDED'
            }
            
            # Reproject to WGS84 for Mapbox
            final_data = np.empty((new_height, new_width), dtype=np.uint8)
            
            reproject(
                mapbox_data,
                final_data,
                src_transform=new_transform,
                src_crs=src.crs,
                dst_transform=profile['transform'],
                dst_crs=profile['crs'],
                resampling=Resampling.average,
                dst_nodata=0
            )
            
            # Write the final data
            with rasterio.open(output_file, 'w', **profile) as dst:
                dst.write(final_data, 1)
            
            # Check file size
            file_size_mb = os.path.getsize(output_file) / 1024 / 1024
            print(f"💾 File size: {file_size_mb:.1f} MB")
            
            if file_size_mb <= 300:
                print(f"✅ Under 300MB limit!")
            else:
                print(f"⚠️  Over 300MB limit")
        
        print(f"\n🎯 Recommended for Mapbox:")
        print(f"📁 Files saved to: {output_dir}")
        print(f"🔧 Next steps:")
        print(f"   1. Upload the smallest file under 300MB to Mapbox Studio")
        print(f"   2. Create a tileset from the raster")
        print(f"   3. Use in map style with expressions like:")
        print(f"      ['>=', ['raster-value'], 20] // Show areas with >20% tree cover")

def create_mapbox_style_examples():
    """Create example Mapbox style expressions for tree cover masking"""
    
    examples = {
        "basic_tree_cover_mask": {
            "description": "Show species data only in areas with >20% tree cover",
            "filter": [">=", ["raster-value"], 20],
            "paint": {
                "circle-opacity": [
                    "case",
                    [">=", ["raster-value"], 20], 0.8,
                    0.2
                ]
            }
        },
        
        "graduated_tree_cover": {
            "description": "Opacity based on tree cover percentage",
            "paint": {
                "circle-opacity": [
                    "interpolate",
                    ["linear"],
                    ["raster-value"],
                    0, 0.1,    # 0% tree cover = very transparent
                    20, 0.4,   # 20% tree cover = semi-transparent  
                    50, 0.7,   # 50% tree cover = mostly opaque
                    80, 1.0    # 80%+ tree cover = fully opaque
                ]
            }
        },
        
        "forest_type_masking": {
            "description": "Different styles for different forest densities",
            "paint": {
                "circle-color": [
                    "case",
                    [">=", ["raster-value"], 60], "#006d2c",  # Dense forest
                    [">=", ["raster-value"], 30], "#41ab5d",  # Medium forest
                    [">=", ["raster-value"], 10], "#a1d99b",  # Light forest
                    "#fee5d9"  # Non-forest
                ]
            }
        },
        
        "species_specific_masking": {
            "description": "Different thresholds for different species",
            "paint": {
                "circle-opacity": [
                    "case",
                    ["==", ["get", "species"], "western_juniper"],
                    [
                        "case", 
                        [">=", ["raster-value"], 5], 0.8,  # Juniper needs less tree cover
                        0.1
                    ],
                    [
                        "case",
                        [">=", ["raster-value"], 30], 0.8,  # Other species need more
                        0.1
                    ]
                ]
            }
        }
    }
    
    output_file = Path("outputs/mapbox_masks/tree_cover_style_examples.json")
    import json
    with open(output_file, 'w') as f:
        json.dump(examples, f, indent=2)
    
    print(f"📝 Style examples saved to: {output_file}")

if __name__ == "__main__":
    create_pnw_tree_cover_mask()
    create_mapbox_style_examples() 