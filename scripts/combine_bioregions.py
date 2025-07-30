#!/usr/bin/env python3
"""
Bioregion Combination Script

Automatically combines all individual bioregion files into a single comprehensive 
GeoJSON file for analysis and web mapping. Designed to work with the standardized
bioregion creation workflow.

Usage:
    python scripts/combine_bioregions.py

Input Files:
    outputs/bioregions/*_constrained.geojson    # Individual bioregion files

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
    bioregion_files = list(bioregions_dir.glob('*_constrained.geojson'))
    
    if not bioregion_files:
        print("❌ No bioregion files found matching pattern '*_constrained.geojson'")
        print(f"   Searched in: {bioregions_dir.absolute()}")
        return
    
    print(f"📁 Found {len(bioregion_files)} bioregion files:")
    
    combined_regions = []
    region_summaries = []
    total_area = 0
    
    for file_path in bioregion_files:
        # Extract region name from filename
        region_name = file_path.stem.replace('_constrained', '')
        
        print(f"   📄 Loading {region_name}...")
        
        try:
            # Load bioregion
            gdf = gpd.read_file(file_path)
            
            if gdf.empty:
                print(f"      ⚠️  Warning: {region_name} contains no features")
                continue
            
            # Ensure consistent CRS
            if gdf.crs != 'EPSG:4326':
                gdf = gdf.to_crs('EPSG:4326')
            
            # Add standardized metadata
            gdf['region_name'] = region_name
            gdf['region_code'] = region_name.upper().replace('_', '')[:8]  # 8-char code
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
                'region_name': region_name,
                'region_code': gdf['region_code'].iloc[0],
                'source_file': file_path.name,
                'feature_count': feature_count,
                'area_km2': float(area_km2),
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
            print(f"      ❌ Error loading {region_name}: {e}")
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

def main():
    """Main function to run bioregion combination process"""
    
    # Run combination
    result = combine_all_bioregions()
    
    if result:
        output_path, summary_path = result
        
        # Run validation
        validate_combined_bioregions(output_path)
        
        print(f"\n🚀 Ready for use!")
        print(f"   Load in QGIS: {output_path}")
        print(f"   Web mapping: Use region_name field for styling")
        print(f"   Analysis: Join with species data using region_code")
    
    else:
        print(f"\n❌ Combination process failed")

if __name__ == "__main__":
    main()