#!/usr/bin/env python3
"""
Create Mapbox mbtiles from precipitation zones for web mapping.

This script converts the precipitation zone GeoJSON to mbtiles format
using Tippecanoe with appropriate zoom levels and layer names.

Usage:
    python scripts/create_precipitation_mbtiles.py
"""

import os
import sys
import subprocess
import json
from datetime import datetime

def check_tippecanoe():
    """Check if tippecanoe is installed."""
    try:
        result = subprocess.run(['tippecanoe', '--version'], 
                              capture_output=True, text=True)
        print(f"Using Tippecanoe: {result.stdout.strip()}")
        return True
    except FileNotFoundError:
        print("ERROR: Tippecanoe is not installed!")
        print("Install with: brew install tippecanoe")
        return False

def create_precipitation_mbtiles():
    """Create mbtiles from precipitation zones GeoJSON."""
    
    # Input and output paths
    input_geojson = "outputs/climate/pnw_precipitation_zones.geojson"
    output_mbtiles = "outputs/climate/pnw_precipitation_zones.mbtiles"
    
    # Check if input exists
    if not os.path.exists(input_geojson):
        print(f"ERROR: Input file not found: {input_geojson}")
        print("Please run process_prism_precipitation_pnw.py first")
        return False
    
    # Create mbtiles with Tippecanoe
    cmd = [
        'tippecanoe',
        '-o', output_mbtiles,
        '-f',  # Force overwrite
        '-z', '14',  # Max zoom level 14
        '-Z', '0',   # Min zoom level 0
        '-l', 'precipitation_zones',  # Layer name
        '--detect-shared-borders',  # Better handling of adjacent polygons
        '--coalesce-densest-as-needed',  # Simplify at lower zooms
        '--extend-zooms-if-still-dropping',  # Keep features at all zooms
        '--attribute-type=zone_code:int',  # Type hints
        '--attribute-type=area_km2:float',
        '--no-feature-limit',  # Include all features
        '--no-tile-size-limit',  # Allow larger tiles if needed
        input_geojson
    ]
    
    print(f"Creating mbtiles with max zoom 14...")
    print(f"Command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"ERROR: Tippecanoe failed!")
            print(f"STDERR: {result.stderr}")
            return False
            
        print(f"Successfully created: {output_mbtiles}")
        
        # Get file size
        size_mb = os.path.getsize(output_mbtiles) / (1024 * 1024)
        print(f"File size: {size_mb:.2f} MB")
        
        # Create metadata file
        metadata = {
            "name": "PNW Precipitation Zones",
            "description": "Annual precipitation zones for Pacific Northwest bioregion analysis",
            "source": "PRISM Climate Group - 30-year normals (1981-2010)",
            "created": datetime.now().isoformat(),
            "layer": "precipitation_zones",
            "zoom_range": "0-14",
            "attributes": {
                "zone_name": "Precipitation zone category",
                "zone_code": "Numeric zone identifier",
                "precip_range": "Precipitation range in inches",
                "area_km2": "Zone area in square kilometers"
            },
            "zones": {
                "very_dry": "0-15 inches (Eastern WA/OR arid)",
                "dry": "15-25 inches (Eastern Cascades)",
                "moderate": "25-50 inches (Western valleys)",
                "wet": "50-80 inches (Western Cascades)",
                "very_wet": "80-150 inches (Coastal forests)",
                "rainforest": ">150 inches (Olympic rainforest)"
            }
        }
        
        metadata_path = output_mbtiles.replace('.mbtiles', '_metadata.json')
        with open(metadata_path, 'w') as f:
            json.dump(metadata, f, indent=2)
        print(f"Created metadata: {metadata_path}")
        
        return True
        
    except Exception as e:
        print(f"ERROR: {str(e)}")
        return False

def main():
    """Main function."""
    print("Creating precipitation zones mbtiles for web mapping...")
    print("=" * 70)
    
    # Check dependencies
    if not check_tippecanoe():
        sys.exit(1)
    
    # Create mbtiles
    if create_precipitation_mbtiles():
        print("\n" + "=" * 70)
        print("Precipitation mbtiles created successfully!")
        print("\nTo use in Mapbox:")
        print("1. Upload outputs/climate/pnw_precipitation_zones.mbtiles to Mapbox Studio")
        print("2. Use layer name 'precipitation_zones' in your style")
        print("3. Style by 'zone_name' or 'zone_code' properties")
        print("4. Max zoom level is 14 for optimal performance")
    else:
        print("\nFailed to create mbtiles")
        sys.exit(1)

if __name__ == "__main__":
    main()