#!/usr/bin/env python3
"""
Create PNW Elevation Mask using GDAL command line tools
Simpler approach to avoid memory issues
"""

import subprocess
import os
from pathlib import Path
import rasterio
import numpy as np

def create_elevation_mask_gdal():
    """Create elevation mask using GDAL tools"""
    
    input_dem = "/Users/JJC/morel-maps/QGIS/static-layers/DEM-western-us-combined.vrt"
    output_dir = Path("outputs/mapbox_masks")
    output_dir.mkdir(exist_ok=True)
    
    # PNW bounding box
    pnw_bounds = [-125.0, 41.5, -110.0, 49.5]  # [west, south, east, north]
    
    print("🏔️ Creating PNW Elevation Mask using GDAL")
    print(f"📍 Region: Oregon, Washington, Idaho")
    print(f"📁 Source: {input_dem}")
    
    # Step 1: Clip to PNW region and resample to 120m (0.001 degrees ≈ 111m)
    temp_clipped = output_dir / "temp_pnw_elevation_clipped.tif"
    temp_resampled = output_dir / "temp_pnw_elevation_resampled.tif"
    final_output = output_dir / "pnw_elevation_120m_mapbox.tif"
    
    print("🔄 Step 1: Clipping to PNW region...")
    clip_cmd = [
        "gdalwarp",
        "-te", str(pnw_bounds[0]), str(pnw_bounds[1]), str(pnw_bounds[2]), str(pnw_bounds[3]),  # target extent
        "-tr", "0.001", "0.001",  # target resolution (≈120m)
        "-r", "average",  # resampling method
        "-co", "COMPRESS=LZW",
        "-co", "TILED=YES",
        "-co", "BLOCKXSIZE=512",
        "-co", "BLOCKYSIZE=512",
        input_dem,
        str(temp_clipped)
    ]
    
    try:
        result = subprocess.run(clip_cmd, capture_output=True, text=True, check=True)
        print("✅ Clipping completed")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error clipping: {e}")
        print(f"STDOUT: {e.stdout}")
        print(f"STDERR: {e.stderr}")
        return
    
    # Step 2: Convert to feet and optimize for Mapbox
    print("🔄 Step 2: Converting to feet and optimizing...")
    
    with rasterio.open(temp_clipped) as src:
        print(f"📊 Clipped size: {src.width} x {src.height}")
        print(f"📊 Elevation range: {src.read(1).min():.0f} - {src.read(1).max():.0f} meters")
        
        # Read data
        elevation_data = src.read(1)
        
        # Convert meters to feet and handle no-data
        elevation_feet = elevation_data.copy().astype(np.float32)
        
        # Handle no-data values
        no_data_mask = (elevation_data == src.nodata) | (elevation_data < -1000)
        
        # Convert to feet
        elevation_feet[~no_data_mask] *= 3.28084
        
        # Set no-data and negative values to 0
        elevation_feet[no_data_mask] = 0
        elevation_feet[elevation_feet < 0] = 0
        
        # Cap at 15,000 feet
        elevation_feet[elevation_feet > 15000] = 15000
        
        # Convert to uint16 for smaller file size
        elevation_uint16 = elevation_feet.astype(np.uint16)
        
        print(f"✅ Converted to feet: {elevation_uint16.min()} - {elevation_uint16.max()} feet")
        
        # Write optimized file
        profile = src.profile.copy()
        profile.update({
            'dtype': 'uint16',
            'nodata': 0,
            'compress': 'lzw',
            'tiled': True,
            'blockxsize': 512,
            'blockysize': 512
        })
        
        with rasterio.open(final_output, 'w', **profile) as dst:
            dst.write(elevation_uint16, 1)
    
    # Clean up temp files
    if temp_clipped.exists():
        temp_clipped.unlink()
    if temp_resampled.exists():
        temp_resampled.unlink()
    
    # Verify results
    file_size_mb = os.path.getsize(final_output) / 1024 / 1024
    print(f"\n🎉 SUCCESS!")
    print(f"📁 File: {final_output}")
    print(f"💾 Size: {file_size_mb:.1f} MB")
    print(f"📐 Resolution: ~120m")
    
    # Final verification
    with rasterio.open(final_output) as verify:
        verify_data = verify.read(1)
        elevation_pixels = verify_data[verify_data > 0]
        
        print(f"\n✅ Verification:")
        print(f"   Size: {verify.width} x {verify.height}")
        print(f"   Elevation range: {verify_data.min()} - {verify_data.max()} feet")
        print(f"   Valid pixels: {len(elevation_pixels):,} ({len(elevation_pixels)/verify_data.size*100:.1f}%)")
        if len(elevation_pixels) > 0:
            print(f"   Average elevation: {elevation_pixels.mean():.0f} feet")
    
    print(f"\n🚀 Ready for Mapbox!")
    print(f"📋 Next steps:")
    print(f"   1. Upload to Mapbox Studio > Tilesets")
    print(f"   2. Use elevation expressions like:")
    print(f"      ['>=', ['raster-value'], 2000] for above 2000ft")
    print(f"      ['<=', ['raster-value'], 6000] for below 6000ft")
    
    # Create combined usage examples
    create_combined_examples(output_dir)

def create_combined_examples(output_dir):
    """Create examples combining tree cover and elevation masking"""
    
    examples = {
        "combined_masking_approach": {
            "description": "Using both tree cover and elevation rasters together",
            "setup": {
                "tree_cover_raster": "pnw_tree_cover_90m_mapbox.tif (0-100% values)",
                "elevation_raster": "pnw_elevation_120m_mapbox.tif (0-15000ft values)"
            }
        },
        
        "mapbox_gl_js_dual_masking": {
            "species_habitat_masking": {
                "description": "Combine elevation and tree cover for species habitat",
                "code": """
// Add both raster sources
map.addSource('tree-cover', {
    type: 'raster',
    url: 'mapbox://your-username.tree-cover-tileset'
});

map.addSource('elevation', {
    type: 'raster', 
    url: 'mapbox://your-username.elevation-tileset'
});

// Species layer with dual masking
map.addLayer({
    id: 'douglas-fir-habitat',
    type: 'circle',
    source: 'species-plots',
    'source-layer': 'douglas_fir_plots',
    paint: {
        'circle-opacity': [
            'case',
            [
                'all',
                ['>=', ['raster-value', 'tree-cover'], 20],  // >20% tree cover
                ['>=', ['raster-value', 'elevation'], 500],   // >500ft elevation
                ['<=', ['raster-value', 'elevation'], 6000]   // <6000ft elevation
            ], 0.8,
            0.1
        ],
        'circle-color': [
            'interpolate', ['linear'], ['raster-value', 'elevation'],
            0, '#2c7fb8',      // Low elevation - blue
            3000, '#41b6c4',   // Mid elevation - teal  
            6000, '#253494'    // High elevation - dark blue
        ]
    }
});
"""
            },
            
            "interactive_controls": {
                "description": "Dynamic controls for both masks",
                "code": """
// HTML controls
<div class="controls">
    <label>Min Tree Cover: <span id="treeCoverValue">20</span>%</label>
    <input type="range" id="treeCoverSlider" min="0" max="100" value="20">
    
    <label>Min Elevation: <span id="elevationMinValue">0</span>ft</label>
    <input type="range" id="elevationMinSlider" min="0" max="10000" value="0">
    
    <label>Max Elevation: <span id="elevationMaxValue">8000</span>ft</label>
    <input type="range" id="elevationMaxSlider" min="0" max="15000" value="8000">
</div>

// JavaScript event handlers
document.getElementById('treeCoverSlider').addEventListener('input', function(e) {
    const minCover = parseInt(e.target.value);
    const minElev = parseInt(document.getElementById('elevationMinSlider').value);
    const maxElev = parseInt(document.getElementById('elevationMaxSlider').value);
    
    updateSpeciesFilter(minCover, minElev, maxElev);
});

function updateSpeciesFilter(minCover, minElev, maxElev) {
    const filter = [
        'all',
        ['>=', ['raster-value', 'tree-cover'], minCover],
        ['>=', ['raster-value', 'elevation'], minElev],
        ['<=', ['raster-value', 'elevation'], maxElev]
    ];
    
    map.setFilter('species-habitat', filter);
}
"""
            }
        },
        
        "species_specific_combinations": {
            "douglas_fir": {
                "tree_cover": ">=20%",
                "elevation": "0-6000ft",
                "description": "Coastal to montane forests"
            },
            "western_juniper": {
                "tree_cover": ">=5%", 
                "elevation": "2000-8000ft",
                "description": "High desert and montane zones"
            },
            "bigleaf_maple": {
                "tree_cover": ">=30%",
                "elevation": "0-4000ft", 
                "description": "Moist lowland and foothill forests"
            }
        }
    }
    
    import json
    output_file = output_dir / "combined_masking_examples.json"
    with open(output_file, 'w') as f:
        json.dump(examples, f, indent=2)
    
    print(f"📝 Combined masking examples: {output_file}")

if __name__ == "__main__":
    create_elevation_mask_gdal() 