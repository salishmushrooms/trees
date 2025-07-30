#!/usr/bin/env python3
"""
Create optimized PNW Elevation Mask for Mapbox
Converts elevation data to feet and optimizes for dynamic elevation masking
"""

import rasterio
from rasterio.warp import reproject, Resampling, transform_bounds
from rasterio.windows import from_bounds
import numpy as np
import os
from pathlib import Path

def create_elevation_mask():
    """Create optimized elevation mask for PNW region"""
    
    input_dem = "/Users/JJC/morel-maps/QGIS/static-layers/DEM-western-us-combined.vrt"
    output_dir = Path("outputs/mapbox_masks")
    output_dir.mkdir(exist_ok=True)
    
    # PNW bounding box (same as tree cover)
    pnw_bounds = {
        'west': -125.0,   # Oregon coast
        'east': -110.0,   # Eastern Idaho  
        'south': 41.5,    # Southern Oregon
        'north': 49.5     # Northern Washington
    }
    
    target_res = 120  # 120m resolution for reasonable file size
    
    print("🏔️ Creating PNW Elevation Mask for Mapbox")
    print(f"🎯 Resolution: {target_res}m")
    print(f"📍 Region: Oregon, Washington, Idaho")
    print(f"📁 Source: {input_dem}")
    
    with rasterio.open(input_dem) as src:
        print(f"📊 Source resolution: ~{abs(src.res[0]) * 111000:.0f}m")
        print(f"📐 Source size: {src.width} x {src.height}")
        print(f"🗺️  Source CRS: {src.crs}")
        
        # Transform bounds to source CRS (already WGS84, so should be same)
        pnw_bounds_src = transform_bounds(
            'EPSG:4326', src.crs,
            pnw_bounds['west'], pnw_bounds['south'],
            pnw_bounds['east'], pnw_bounds['north']
        )
        
        print(f"🎯 PNW bounds: {pnw_bounds_src}")
        
        # Create window for PNW region
        window = from_bounds(*pnw_bounds_src, src.transform)
        print(f"🔲 Extraction window: {window}")
        
        # Read the PNW subset
        print("📖 Reading PNW elevation data...")
        pnw_data = src.read(1, window=window)
        pnw_transform = src.window_transform(window)
        
        print(f"✅ Loaded: {pnw_data.shape}")
        print(f"📊 Elevation range: {pnw_data.min():.0f} - {pnw_data.max():.0f} meters")
        
        # Calculate elevation statistics (excluding no-data)
        valid_elevations = pnw_data[pnw_data != -32768]
        if len(valid_elevations) > 0:
            print(f"🏔️ Elevation stats:")
            print(f"   Mean: {valid_elevations.mean():.0f}m ({valid_elevations.mean() * 3.28084:.0f}ft)")
            print(f"   Median: {np.median(valid_elevations):.0f}m ({np.median(valid_elevations) * 3.28084:.0f}ft)")
            print(f"   Std: {valid_elevations.std():.0f}m")
            
            # Show elevation zones
            print(f"🌊 Sea level (0-500m): {(valid_elevations <= 500).sum() / len(valid_elevations) * 100:.1f}%")
            print(f"🏞️ Lowland (500-1000m): {((valid_elevations > 500) & (valid_elevations <= 1000)).sum() / len(valid_elevations) * 100:.1f}%")
            print(f"⛰️  Montane (1000-2000m): {((valid_elevations > 1000) & (valid_elevations <= 2000)).sum() / len(valid_elevations) * 100:.1f}%")
            print(f"🏔️ Alpine (>2000m): {(valid_elevations > 2000).sum() / len(valid_elevations) * 100:.1f}%")
        
        # Calculate new dimensions for target resolution
        # Current resolution is ~30m, so scale accordingly
        current_res = abs(src.res[0]) * 111000  # Convert degrees to meters approximately
        scale_factor = current_res / target_res  # Fix: divide current by target for downsampling
        new_height = int(pnw_data.shape[0] / scale_factor)
        new_width = int(pnw_data.shape[1] / scale_factor)
        
        print(f"\n🔄 Resampling to {target_res}m resolution...")
        print(f"📐 New size: {new_width} x {new_height}")
        
        # Resample data using average (good for elevation)
        resampled_data = np.empty((new_height, new_width), dtype=np.float32)
        
        new_transform = rasterio.transform.from_bounds(
            *pnw_bounds_src, new_width, new_height
        )
        
        reproject(
            pnw_data,
            resampled_data,
            src_transform=pnw_transform,
            src_crs=src.crs,
            dst_transform=new_transform,
            dst_crs=src.crs,
            resampling=Resampling.average,
            src_nodata=-32768,
            dst_nodata=-32768
        )
        
        # Optimize for Mapbox elevation masking
        print("🎨 Optimizing for Mapbox...")
        
        # Convert to feet and optimize data range
        mapbox_data = resampled_data.copy()
        
        # Convert meters to feet
        mapbox_data[resampled_data != -32768] *= 3.28084
        
        # Handle no-data and negative elevations
        mapbox_data[resampled_data == -32768] = 0  # No data = 0 (sea level)
        mapbox_data[mapbox_data < 0] = 0  # Negative elevations = 0
        
        # Cap at reasonable maximum (15,000 feet)
        mapbox_data[mapbox_data > 15000] = 15000
        
        # Convert to uint16 for smaller file size (0-65535 range)
        mapbox_data_int = mapbox_data.astype(np.uint16)
        
        print(f"✅ Optimized elevation range: {mapbox_data_int.min()} - {mapbox_data_int.max()} feet")
        
        # Create final output file
        output_file = output_dir / "pnw_elevation_120m_mapbox.tif"
        
        # Profile optimized for Mapbox
        profile = {
            'driver': 'GTiff',
            'dtype': 'uint16',
            'nodata': 0,
            'width': new_width,
            'height': new_height,
            'count': 1,
            'crs': 'EPSG:4326',
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
        
        print("💾 Writing final elevation file...")
        with rasterio.open(output_file, 'w', **profile) as dst:
            dst.write(mapbox_data_int, 1)
        
        # Verify results
        file_size_mb = os.path.getsize(output_file) / 1024 / 1024
        print(f"\n🎉 SUCCESS!")
        print(f"📁 File: {output_file}")
        print(f"💾 Size: {file_size_mb:.1f} MB")
        print(f"📐 Resolution: {target_res}m")
        print(f"🗺️  Coverage: Oregon, Washington, Idaho")
        
        # Final verification
        with rasterio.open(output_file) as verify:
            verify_data = verify.read(1)
            elevation_pixels = verify_data[verify_data > 0]
            
            print(f"\n✅ Verification:")
            print(f"   Elevation range: {verify_data.min()} - {verify_data.max()} feet")
            print(f"   Valid pixels: {len(elevation_pixels):,} ({len(elevation_pixels)/verify_data.size*100:.1f}%)")
            if len(elevation_pixels) > 0:
                print(f"   Average elevation: {elevation_pixels.mean():.0f} feet")
        
        print(f"\n🚀 Ready for Mapbox!")
        print(f"📋 Upload steps:")
        print(f"   1. Go to Mapbox Studio > Tilesets")
        print(f"   2. Upload: {output_file}")
        print(f"   3. Create tileset from raster")
        print(f"   4. Use in map style with elevation expressions")
        
        # Create usage examples
        create_elevation_examples(output_dir)

def create_elevation_examples(output_dir):
    """Create comprehensive usage examples for elevation masking"""
    
    examples = {
        "mapbox_gl_js_elevation_masking": {
            "basic_elevation_ranges": {
                "description": "Basic elevation range masking",
                "code": {
                    "filter": [
                        "all",
                        [">=", ["raster-value"], 1000],  # Above 1000ft
                        ["<=", ["raster-value"], 7000]   # Below 7000ft
                    ],
                    "paint": {
                        "circle-opacity": [
                            "case",
                            [
                                "all",
                                [">=", ["raster-value"], 1000],
                                ["<=", ["raster-value"], 7000]
                            ], 0.8,
                            0.1
                        ]
                    }
                }
            },
            
            "species_elevation_zones": {
                "description": "Different elevation zones for different species",
                "code": {
                    "paint": {
                        "circle-opacity": [
                            "case",
                            ["==", ["get", "species"], "douglas_fir"],
                            [
                                "case",
                                [
                                    "all",
                                    [">=", ["raster-value"], 0],
                                    ["<=", ["raster-value"], 6000]
                                ], 0.8, 0.1
                            ],
                            ["==", ["get", "species"], "western_juniper"],
                            [
                                "case",
                                [
                                    "all", 
                                    [">=", ["raster-value"], 2000],
                                    ["<=", ["raster-value"], 8000]
                                ], 0.8, 0.1
                            ],
                            ["==", ["get", "species"], "bigleaf_maple"],
                            [
                                "case",
                                [
                                    "all",
                                    [">=", ["raster-value"], 0],
                                    ["<=", ["raster-value"], 4000]
                                ], 0.8, 0.1
                            ],
                            0.5
                        ]
                    }
                }
            },
            
            "elevation_gradient_coloring": {
                "description": "Color points based on elevation zones",
                "code": {
                    "paint": {
                        "circle-color": [
                            "interpolate", ["linear"], ["raster-value"],
                            0, "#2c7fb8",      # Sea level - blue
                            1000, "#7fcdbb",   # Lowland - light blue
                            3000, "#41b6c4",   # Foothill - teal
                            6000, "#253494",   # Montane - dark blue
                            9000, "#081d58"    # Alpine - very dark blue
                        ]
                    }
                }
            },
            
            "combined_tree_elevation_masking": {
                "description": "Combine tree cover and elevation masking",
                "code": {
                    "paint": {
                        "circle-opacity": [
                            "case",
                            [
                                "all",
                                [">=", ["raster-value"], 1000],      # Elevation filter
                                ["<=", ["raster-value"], 7000],
                                [">=", ["get", "tree_cover"], 20]    # Tree cover filter (from separate raster)
                            ], 0.8,
                            0.1
                        ]
                    }
                }
            }
        },
        
        "mapbox_studio_elevation_expressions": {
            "coastal_filter": "['<=', ['raster-value'], 1000]",
            "montane_filter": "[['>=', ['raster-value'], 3000], ['<=', ['raster-value'], 8000]]",
            "alpine_filter": "['>=', ['raster-value'], 8000]",
            "elevation_opacity": "['interpolate', ['linear'], ['raster-value'], 0, 0.3, 5000, 1.0]"
        },
        
        "common_elevation_zones": {
            "coastal": "0-1000 feet",
            "valley_foothill": "1000-3000 feet", 
            "montane": "3000-6000 feet",
            "subalpine": "6000-9000 feet",
            "alpine": ">9000 feet"
        },
        
        "species_elevation_preferences": {
            "douglas_fir": "0-6000 feet (prefers 1000-4000)",
            "bigleaf_maple": "0-4000 feet (prefers 500-2500)",
            "western_juniper": "2000-8000 feet (prefers 3000-6000)",
            "ponderosa_pine": "1000-7000 feet (prefers 2000-5000)",
            "western_hemlock": "0-5000 feet (prefers coastal to 3000)"
        }
    }
    
    import json
    output_file = output_dir / "elevation_masking_examples.json"
    with open(output_file, 'w') as f:
        json.dump(examples, f, indent=2)
    
    print(f"📝 Elevation examples: {output_file}")

if __name__ == "__main__":
    create_elevation_mask() 