#!/usr/bin/env python3
"""
Bioregion Combination Script

Automatically combines all individual bioregion files into a single comprehensive 
GeoJSON file for analysis and web mapping. Designed to work with the standardized
bioregion creation workflow.

Usage:
    python scripts/combine_bioregions.py

Input Files (configured in BIOREGION_FILES list):
    outputs/bioregions/coastal_forest_shapefile_constrained.geojson
    outputs/bioregions/eastern_cascades_species_based.geojson
    outputs/bioregions/high_cascades_species_based.geojson
    outputs/bioregions/olympic_mountains_constrained.geojson

Output Files:
    outputs/bioregions/all_bioregions_combined.geojson    # Combined bioregions
    outputs/bioregions/bioregions_combination_summary.json    # Summary statistics

Author: Pacific Northwest Forest Habitat Mapping Project
Created: Based on coastal forest bioregion implementation
Updated: 2025-07-29
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

# =============================================================================
# BIOREGION FILES CONFIGURATION
# =============================================================================
# Add or remove bioregion files to include in the combination process
# Files should be located in outputs/bioregions/
# Format: (filename, desired_name) - if desired_name is None, it will be auto-generated
BIOREGION_FILES = [
    ('coastal_forest_shapefile_constrained.geojson', 'Coastal Forest'),
    ('eastern_cascades_species_based.geojson', 'Eastern Cascades'),
    ('high_cascades_species_based.geojson', 'High Cascades'),
    ('olympic_mountains_constrained.geojson', 'Olympic Mountains'),
    # Add new bioregion files here as needed:
    # ('new_bioregion_constrained.geojson', 'New Bioregion Name'),
    # ('another_bioregion.geojson', None),  # Will use auto-generated name
]
# =============================================================================

def create_display_name(file_path):
    """Create human-readable display name from bioregion file path"""
    
    # Extract base name and remove suffixes
    base_name = file_path.stem
    
    # Remove common suffixes that shouldn't be in the display name
    suffixes_to_remove = ['_constrained', '_species_based', '_shapefile', '_broad']
    for suffix in suffixes_to_remove:
        base_name = base_name.replace(suffix, '')
    
    # Define name mappings for better display names
    name_mappings = {
        'coastal_forest': 'Coastal Forest',
        'olympic_mountains': 'Olympic Mountains',
        'eastern_cascades': 'Eastern Cascades', 
        'high_cascades': 'High Cascades',
        'riparian_forest': 'Riparian Forest',
        'alpine_forest': 'Alpine Forest'
    }
    
    # Return mapped name or formatted version
    if base_name in name_mappings:
        return name_mappings[base_name]
    else:
        # Fallback: title case with underscores as spaces
        return base_name.replace('_', ' ').title()

def combine_all_bioregions():
    """
    Automatically combine all bioregion files in the bioregions directory
    
    This function:
    1. Scans for all files matching '*_constrained.geojson' pattern
    2. Loads each bioregion and adds standardized metadata
    3. Combines into single GeoDataFrame with region identifiers
    4. Calculates summary statistics for each region
    5. Saves combined file and summary report
    """
    
    print("=== Combining Pacific Northwest Forest Bioregions ===")
    
    bioregions_dir = Path('outputs/bioregions')
    
    # Build file paths from hardcoded configuration
    bioregion_files = []
    missing_files = []
    
    for entry in BIOREGION_FILES:
        # Handle both tuple and string formats for backwards compatibility
        if isinstance(entry, tuple):
            filename, desired_name = entry
        else:
            filename = entry
            desired_name = None
            
        file_path = bioregions_dir / filename
        if file_path.exists():
            bioregion_files.append((file_path, desired_name))
        else:
            missing_files.append(filename)
    
    if missing_files:
        print("⚠️  Warning: Some configured files not found:")
        for missing in missing_files:
            print(f"   • {missing}")
        print()
    
    if not bioregion_files:
        print("❌ No bioregion files found from configuration:")
        for entry in BIOREGION_FILES:
            if isinstance(entry, tuple):
                filename, _ = entry
            else:
                filename = entry
            print(f"   • {filename}")
        print(f"   Searched in: {bioregions_dir.absolute()}")
        return
    
    print(f"📁 Found {len(bioregion_files)} bioregion files:")
    
    combined_regions = []
    region_summaries = []
    total_area = 0
    
    for file_entry in bioregion_files:
        # Unpack file path and desired name
        if isinstance(file_entry, tuple):
            file_path, desired_name = file_entry
        else:
            file_path = file_entry
            desired_name = None
            
        # Use hardcoded name if provided, otherwise generate from filename
        if desired_name:
            display_name = desired_name
        else:
            display_name = create_display_name(file_path)
            
        region_code = file_path.stem.replace('_constrained', '').replace('_species_based', '').replace('_shapefile', '')
        
        print(f"   📄 Loading {display_name} ({region_code})...")
        
        try:
            # Load bioregion
            gdf = gpd.read_file(file_path)
            
            if gdf.empty:
                print(f"      ⚠️  Warning: {display_name} contains no features")
                continue
            
            # Ensure consistent CRS
            if gdf.crs != 'EPSG:4326':
                gdf = gdf.to_crs('EPSG:4326')
            
            # Reduce precision for web mapping (~10m precision)
            print(f"      🔧 Reducing precision to ~10m...")
            
            # Simplify geometry (remove unnecessary vertices)
            tolerance_degrees = 0.0001  # ~10m at mid-latitudes
            gdf['geometry'] = gdf.geometry.simplify(tolerance_degrees, preserve_topology=True)
            
            # Round coordinates to ~10m precision (5 decimal places)
            from shapely.ops import transform
            import math
            
            def round_coordinates(geom):
                """Round coordinates to 5 decimal places (~1m precision)"""
                def round_coords(x, y, z=None):
                    rounded_x = round(x, 5)
                    rounded_y = round(y, 5)
                    return (rounded_x, rounded_y) if z is None else (rounded_x, rounded_y, round(z, 5))
                return transform(round_coords, geom)
            
            gdf['geometry'] = gdf['geometry'].apply(round_coordinates)
            
            # Add standardized metadata
            gdf['region_name'] = display_name  # Human-readable name for Mapbox
            gdf['region_code'] = region_code.upper().replace('_', '')[:8]  # 8-char code
            gdf['source_file'] = file_path.name
            gdf['bioregion_type'] = 'forest_bioregion'
            
            # Calculate area statistics
            area_km2 = gdf.geometry.area.sum() * 111 * 111  # Approximate conversion
            gdf['area_km2'] = area_km2
            total_area += area_km2
            
            # Count features
            feature_count = len(gdf)
            
            # Add processing timestamp
            gdf['combined_date'] = datetime.now().isoformat()
            
            combined_regions.append(gdf)
            
            # Store region summary
            region_summary = {
                'region_name': display_name,
                'region_code': gdf['region_code'].iloc[0],
                'source_file': file_path.name,
                'feature_count': feature_count,
                'area_km2': float(area_km2),
                'precision_reduced': True,
                'simplification_tolerance_degrees': tolerance_degrees,
                'coordinate_precision_decimal_places': 5,
                'bounds': {
                    'minx': float(gdf.bounds.minx.min()),
                    'miny': float(gdf.bounds.miny.min()),
                    'maxx': float(gdf.bounds.maxx.max()),
                    'maxy': float(gdf.bounds.maxy.max())
                }
            }
            region_summaries.append(region_summary)
            
            print(f"      ✅ Loaded {feature_count} features, {area_km2:,.0f} km²")
            
        except Exception as e:
            print(f"      ❌ Error loading {display_name}: {e}")
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
    
    # Create comprehensive summary
    summary = {
        'dataset_info': {
            'name': 'Pacific Northwest Forest Bioregions',
            'version': '1.0',
            'description': 'Combined forest bioregions for habitat mapping and species distribution analysis',
            'coordinate_system': 'EPSG:4326 (WGS84)',
            'created_date': datetime.now().isoformat()
        },
        'summary_statistics': {
            'total_regions': len(combined_regions),
            'total_features': len(all_bioregions),
            'total_area_km2': float(total_area),
            'combined_bounds': {
                'minx': float(all_bioregions.bounds.minx.min()),
                'miny': float(all_bioregions.bounds.miny.min()),
                'maxx': float(all_bioregions.bounds.maxx.max()),
                'maxy': float(all_bioregions.bounds.maxy.max())
            }
        },
        'individual_regions': region_summaries,
        'output_files': {
            'combined_bioregions': str(output_path),
            'summary_report': 'bioregions_combination_summary.json'
        },
        'usage_notes': {
            'web_mapping': 'Use region_name field for layer styling and filtering',
            'analysis': 'Use region_code for joining with species data',
            'validation': 'Check area_km2 field for bioregion size statistics'
        }
    }
    
    # Save summary report
    summary_path = bioregions_dir / 'bioregions_combination_summary.json'
    try:
        with open(summary_path, 'w') as f:
            json.dump(summary, f, indent=2)
        print(f"✅ Saved summary report: {summary_path}")
    except Exception as e:
        print(f"❌ Error saving summary: {e}")
    
    # Final summary output
    print(f"\n🎉 Successfully combined bioregions!")
    print(f"   📊 Regions: {len(combined_regions)}")
    print(f"   🗺️  Features: {len(all_bioregions):,}")
    print(f"   📐 Total area: {total_area:,.0f} km²")
    print(f"   📄 Output: {output_path.name}")
    
    # List individual regions
    print(f"\n📋 Included regions:")
    for summary in region_summaries:
        print(f"   • {summary['region_name']}: {summary['area_km2']:,.0f} km² ({summary['feature_count']} features)")
    
    return output_path, summary_path

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
        required_columns = ['region_name', 'region_code', 'area_km2', 'geometry']
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
            area = gdf[gdf['region_name'] == region]['area_km2'].sum()
            print(f"      • {region}: {count} features, {area:,.0f} km²")
        
        return True
        
    except Exception as e:
        print(f"   ❌ Validation error: {e}")
        return False

def create_bioregions_mbtiles(geojson_file=None, max_zoom=14):
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
        '--include', 'region_code', 
        '--include', 'area_km2',
        '--include', 'bioregion_type',
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
        output_path, summary_path = result
        
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
            print(f"   📊 Summary: {summary_path.name}")
            
            print(f"\n🚀 Ready for use!")
            print(f"   Load in QGIS: {output_path}")
            print(f"   Web mapping: Use region_name field for styling")
            print(f"   Analysis: Join with species data using region_code")
            
            if mbtiles_path:
                print(f"\n📱 Mapbox Integration:")
                print(f"   1. Upload {mbtiles_path.name} to Mapbox Studio")
                print(f"   2. Create style with 'bioregions' layer")
                print(f"   3. Use region_name property for styling:")
                print(f"      - 'Coastal Forest', 'Olympic Mountains', 'Eastern Cascades', 'High Cascades'")
                print(f"   4. Geometry simplified to ~10m precision for web performance")
        else:
            print(f"\n⚠️  Validation failed - skipping mbtiles creation")
    
    else:
        print(f"\n❌ Combination process failed")

if __name__ == "__main__":
    main()