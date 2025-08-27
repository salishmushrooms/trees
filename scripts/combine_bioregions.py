#!/usr/bin/env python3
"""
Simplified Bioregion Combination Script

Combines bioregion files containing only region_name for Mapbox styling.
Each bioregion is now a single MultiPolygon feature for optimal performance.

Usage:
    python scripts/combine_bioregions.py

Output Files:
    outputs/bioregions/all_bioregions_combined.geojson    # Combined bioregions for Mapbox

Author: Pacific Northwest Forest Habitat Mapping Project
Updated: 2025-08-27 - Simplified for Mapbox styling only
"""

import geopandas as gpd
import pandas as pd
from pathlib import Path
import json
from datetime import datetime
import warnings
import subprocess
import shutil

warnings.filterwarnings('ignore')

# Bioregion files to combine (filename only - region_name comes from the file)
BIOREGION_FILES = [
    'eastern_cascades_constrained.geojson',
    'western_cascades_constrained.geojson', 
    'urban_range_constrained.geojson',
    'high_cascades_constrained.geojson',
    'olympic_mountains_constrained.geojson',
    'coastal_forest_constrained.geojson',
    'coast_range_constrained.geojson'
]


def combine_all_bioregions():
    """
    Combine all bioregion files into a single GeoJSON for Mapbox
    Each file should contain a single MultiPolygon feature with region_name property
    """
    
    print("=== Combining Pacific Northwest Forest Bioregions ===")
    
    bioregions_dir = Path('outputs/bioregions')
    
    # Find existing bioregion files
    existing_files = []
    missing_files = []
    
    for filename in BIOREGION_FILES:
        file_path = bioregions_dir / filename
        if file_path.exists():
            existing_files.append(file_path)
        else:
            missing_files.append(filename)
    
    if missing_files:
        print("⚠️  Warning: Some configured files not found:")
        for missing in missing_files:
            print(f"   • {missing}")
        print()
    
    if not existing_files:
        print("❌ No bioregion files found")
        return
    
    print(f"📁 Found {len(existing_files)} bioregion files:")
    
    combined_regions = []
    
    for file_path in existing_files:
        print(f"   📄 Loading {file_path.name}...")
        
        try:
            # Load bioregion (should be single MultiPolygon with region_name)
            gdf = gpd.read_file(file_path)
            
            if gdf.empty:
                print(f"      ⚠️  Warning: {file_path.name} contains no features")
                continue
            
            # Ensure consistent CRS
            if gdf.crs != 'EPSG:4326':
                gdf = gdf.to_crs('EPSG:4326')
            
            # Validate required columns
            if 'region_name' not in gdf.columns:
                print(f"      ❌ Warning: {file_path.name} missing region_name column")
                continue
            
            combined_regions.append(gdf)
            
            region_name = gdf['region_name'].iloc[0]
            print(f"      ✅ Loaded {region_name}")
            
        except Exception as e:
            print(f"      ❌ Error loading {file_path.name}: {e}")
            continue
    
    if not combined_regions:
        print("❌ No valid bioregion files could be loaded")
        return
    
    print(f"\n🔄 Combining {len(combined_regions)} bioregions...")
    
    # Combine all regions
    all_bioregions = pd.concat(combined_regions, ignore_index=True)
    
    # Add master dataset metadata
    all_bioregions['dataset'] = 'PNW_Forest_Bioregions'
    all_bioregions['dataset_version'] = '1.0'
    all_bioregions['created_date'] = datetime.now().isoformat()
    all_bioregions['coordinate_system'] = 'EPSG:4326'
    
    # Ensure geometry column is properly set
    if 'geometry' not in all_bioregions.columns:
        print("❌ Error: No geometry column found in combined data")
        return
    
    # Save combined file
    output_path = bioregions_dir / 'all_bioregions_combined.geojson'
    try:
        all_bioregions.to_file(output_path, driver='GeoJSON')
        print(f"✅ Saved combined bioregions: {output_path}")
    except Exception as e:
        print(f"❌ Error saving combined file: {e}")
        return
    
    # Final summary output
    print(f"\n🎉 Successfully combined bioregions!")
    print(f"   📊 Regions: {len(combined_regions)}")
    print(f"   🗺️  Features: {len(all_bioregions)} (1 per region)")
    print(f"   📄 Output: {output_path.name}")
    
    # List individual regions
    print(f"\n📋 Included regions:")
    for region in all_bioregions['region_name']:
        print(f"   • {region}")
    
    return output_path

def validate_combined_bioregions(combined_file=None):
    """
    Validate the combined bioregions file for common issues
    
    Args:
        combined_file: Path to combined bioregions file (optional)
    """
    
    if combined_file is None:
        combined_file = Path('outputs/bioregions/all_bioregions_combined.geojson')
    
    if not combined_file.exists():
        print(f"❌ Combined file not found: {combined_file}")
        return False
    
    print(f"\n🔍 Validating combined bioregions: {combined_file.name}")
    
    try:
        gdf = gpd.read_file(combined_file)
        
        # Basic validation checks
        print(f"   ✅ File loads successfully")
        print(f"   ✅ Contains {len(gdf)} features")
        print(f"   ✅ CRS: {gdf.crs}")
        
        # Check required columns
        required_columns = ['region_name', 'geometry']
        missing_columns = [col for col in required_columns if col not in gdf.columns]
        
        if missing_columns:
            print(f"   ⚠️  Missing columns: {missing_columns}")
        else:
            print(f"   ✅ All required columns present")
        
        # Check for empty geometries
        empty_geoms = gdf.geometry.is_empty.sum()
        if empty_geoms > 0:
            print(f"   ⚠️  {empty_geoms} empty geometries found")
        else:
            print(f"   ✅ No empty geometries")
        
        # Check for valid geometries
        invalid_geoms = (~gdf.geometry.is_valid).sum()
        if invalid_geoms > 0:
            print(f"   ⚠️  {invalid_geoms} invalid geometries found")
        else:
            print(f"   ✅ All geometries valid")
        
        # Region summary
        regions = gdf['region_name'].value_counts()
        print(f"   📊 Regions breakdown:")
        for region, count in regions.items():
            print(f"      • {region}: {count} features")
        
        return True
        
    except Exception as e:
        print(f"   ❌ Validation error: {e}")
        return False

def create_bioregions_mbtiles(geojson_file=None, max_zoom=13):
    """
    Create Mapbox mbtiles from combined bioregions using Tippecanoe
    
    Args:
        geojson_file: Path to combined bioregions GeoJSON (optional)
        max_zoom: Maximum zoom level for tiles (default: 14)
    """
    
    if geojson_file is None:
        geojson_file = Path('outputs/bioregions/all_bioregions_combined.geojson')
    
    if not geojson_file.exists():
        print(f"❌ GeoJSON file not found: {geojson_file}")
        return None
    
    print(f"\n🗺️  Creating Mapbox tiles from bioregions...")
    
    # Check if tippecanoe is available
    if not shutil.which('tippecanoe'):
        print("❌ Tippecanoe not found. Install with:")
        print("   brew install tippecanoe")
        return None
    
    # Output path for mbtiles
    mbtiles_dir = Path('outputs/mapbox_mbtiles')
    mbtiles_dir.mkdir(exist_ok=True)
    mbtiles_file = mbtiles_dir / 'pnw_forest_bioregions.mbtiles'
    
    # Remove existing mbtiles file if it exists
    if mbtiles_file.exists():
        mbtiles_file.unlink()
        print(f"   🗑️  Removed existing {mbtiles_file.name}")
    
    # Tippecanoe command
    tippecanoe_cmd = [
        'tippecanoe',
        '--output', str(mbtiles_file),
        '--maximum-zoom', str(max_zoom),
        '--minimum-zoom', '0',
        '--base-zoom', str(max_zoom),
        '--force',  # Overwrite existing files
        '--drop-densest-as-needed',  # Automatically drop features at high zoom levels
        '--extend-zooms-if-still-dropping',  # Extend zoom levels if needed
        '--layer', 'bioregions',  # Layer name in the mbtiles
        '--attribution', 'Pacific Northwest Forest Habitat Mapping Project',
        '--name', 'PNW Forest Bioregions',
        '--description', 'Pacific Northwest Forest Bioregions for habitat mapping and species distribution analysis',
        # Feature properties to include
        '--include', 'region_name',
        str(geojson_file)
    ]
    
    print(f"   ⚙️  Running tippecanoe...")
    print(f"      Max zoom: {max_zoom}")
    print(f"      Output: {mbtiles_file.name}")
    
    try:
        # Run tippecanoe
        result = subprocess.run(
            tippecanoe_cmd,
            capture_output=True,
            text=True,
            check=True
        )
        
        print(f"   ✅ Successfully created mbtiles file")
        
        # Get file size
        file_size_mb = mbtiles_file.stat().st_size / (1024 * 1024)
        print(f"      File size: {file_size_mb:.1f} MB")
        
        # Update summary with mbtiles info
        summary_file = geojson_file.parent / 'bioregions_combination_summary.json'
        if summary_file.exists():
            with open(summary_file, 'r') as f:
                summary = json.load(f)
            
            summary['output_files']['mbtiles'] = str(mbtiles_file)
            summary['mbtiles_info'] = {
                'max_zoom': max_zoom,
                'layer_name': 'bioregions',
                'file_size_mb': round(file_size_mb, 1),
                'created_date': datetime.now().isoformat(),
                'tippecanoe_version': get_tippecanoe_version()
            }
            
            with open(summary_file, 'w') as f:
                json.dump(summary, f, indent=2)
            
            print(f"   📄 Updated summary with mbtiles info")
        
        return mbtiles_file
        
    except subprocess.CalledProcessError as e:
        print(f"   ❌ Tippecanoe failed:")
        print(f"      Error: {e.stderr}")
        return None
    except Exception as e:
        print(f"   ❌ Error creating mbtiles: {e}")
        return None

def get_tippecanoe_version():
    """Get tippecanoe version for metadata"""
    try:
        result = subprocess.run(['tippecanoe', '--version'], capture_output=True, text=True)
        return result.stdout.strip()
    except:
        return "unknown"

def main():
    """Main function to run bioregion combination and mbtiles creation process"""
    
    print("=== Pacific Northwest Forest Bioregions Processing ===")
    
    # Run combination
    result = combine_all_bioregions()
    
    if result:
        output_path = result
        
        # Run validation
        is_valid = validate_combined_bioregions(output_path)
        
        if is_valid:
            # Create mbtiles
            mbtiles_path = create_bioregions_mbtiles(output_path, max_zoom=14)
            
            print(f"\n🎉 Processing complete!")
            print(f"   📄 GeoJSON: {output_path.name}")
            if mbtiles_path:
                print(f"   🗺️  MBTiles: {mbtiles_path.name}")
                print(f"   📤 Upload to Mapbox Studio: {mbtiles_path}")
            
            print(f"\n🚀 Ready for use!")
            print(f"   Load in QGIS: {output_path}")
            print(f"   Web mapping: Use region_name field for styling")
            
            if mbtiles_path:
                print(f"\n📱 Mapbox Integration:")
                print(f"   1. Upload {mbtiles_path.name} to Mapbox Studio")
                print(f"   2. Create style with 'bioregions' layer")
                print(f"   3. Use region_name property for styling")
                print(f"   4. Each bioregion is a single MultiPolygon for optimal performance")
        else:
            print(f"\n⚠️  Validation failed - skipping mbtiles creation")
    
    else:
        print(f"\n❌ Combination process failed")

if __name__ == "__main__":
    main()