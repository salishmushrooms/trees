#!/usr/bin/env python3
"""
Fix Invalid Geometries in Bioregion Files

This script checks for and repairs invalid geometries in bioregion geojson files.
Common issues include self-intersections, duplicate vertices, and topology errors.

Usage:
    python scripts/fix_bioregion_geometries.py
"""

import geopandas as gpd
import pandas as pd
from pathlib import Path
from shapely.geometry import Point, Polygon, MultiPolygon
from shapely.ops import unary_union
from shapely.validation import make_valid
import warnings
import json

warnings.filterwarnings('ignore')

# Files to check and fix
BIOREGION_FILES = [
    'eastern_cascades_constrained.geojson',
    'western_cascades_constrained.geojson', 
    'urban_range_constrained.geojson',
    'high_cascades_constrained.geojson',
    'olympic_mountains_constrained.geojson',
    'coastal_forest_shapefile_constrained.geojson',
    'coast_range_constrained.geojson'
]

def fix_geometry(geom):
    """Fix invalid geometries using shapely's make_valid function"""
    if geom is None or geom.is_empty:
        return None
    
    # If geometry is already valid, return as-is
    if geom.is_valid:
        return geom
    
    try:
        # Use make_valid to fix common issues
        fixed_geom = make_valid(geom)
        
        # If make_valid resulted in a geometry collection, try to extract largest polygon
        if hasattr(fixed_geom, 'geom_type') and fixed_geom.geom_type == 'GeometryCollection':
            polygons = []
            for g in fixed_geom.geoms:
                if g.geom_type in ['Polygon', 'MultiPolygon']:
                    polygons.append(g)
            
            if polygons:
                # Return the largest polygon/multipolygon
                largest = max(polygons, key=lambda x: x.area)
                return largest
            else:
                print(f"Warning: Could not extract polygon from geometry collection")
                return None
        
        return fixed_geom
        
    except Exception as e:
        print(f"Error fixing geometry: {e}")
        return None

def check_and_fix_bioregion_file(file_path):
    """Check a single bioregion file for invalid geometries and fix them"""
    
    if not file_path.exists():
        print(f"⚠️  File not found: {file_path}")
        return {'status': 'not_found', 'file': file_path.name}
    
    print(f"\n📍 Checking {file_path.name}...")
    
    try:
        # Load the bioregion file
        gdf = gpd.read_file(file_path)
        print(f"   Loaded {len(gdf)} features")
        
        # Check for invalid geometries
        invalid_mask = ~gdf.geometry.is_valid
        invalid_count = invalid_mask.sum()
        
        if invalid_count == 0:
            print(f"   ✅ All geometries are valid")
            return {'status': 'valid', 'file': file_path.name, 'features': len(gdf)}
        
        print(f"   ⚠️  Found {invalid_count} invalid geometries")
        
        # Get details about invalid geometries
        invalid_reasons = []
        for idx, geom in gdf[invalid_mask].geometry.items():
            try:
                reason = geom.is_valid_reason if hasattr(geom, 'is_valid_reason') else "Unknown"
                invalid_reasons.append(f"Feature {idx}: {reason}")
            except:
                invalid_reasons.append(f"Feature {idx}: Could not determine reason")
        
        print(f"   Invalid geometry reasons:")
        for reason in invalid_reasons[:5]:  # Show first 5
            print(f"     - {reason}")
        if len(invalid_reasons) > 5:
            print(f"     ... and {len(invalid_reasons) - 5} more")
        
        # Fix invalid geometries
        print(f"   🔧 Fixing invalid geometries...")
        original_crs = gdf.crs
        
        fixed_geometries = []
        fixed_count = 0
        removed_count = 0
        
        for idx, row in gdf.iterrows():
            geom = row.geometry
            
            if not geom.is_valid:
                fixed_geom = fix_geometry(geom)
                if fixed_geom is not None and not fixed_geom.is_empty:
                    # Create new row with fixed geometry
                    new_row = row.copy()
                    new_row.geometry = fixed_geom
                    fixed_geometries.append(new_row)
                    fixed_count += 1
                else:
                    print(f"     ⚠️  Removed feature {idx} - could not fix geometry")
                    removed_count += 1
            else:
                # Keep valid geometries as-is
                fixed_geometries.append(row)
        
        # Create new GeoDataFrame with fixed geometries
        if fixed_geometries:
            fixed_gdf = gpd.GeoDataFrame(fixed_geometries, crs=original_crs)
            
            # Verify all geometries are now valid
            final_invalid = (~fixed_gdf.geometry.is_valid).sum()
            
            if final_invalid == 0:
                # Create backup of original file
                backup_path = file_path.with_suffix('.backup.geojson')
                if not backup_path.exists():
                    print(f"   💾 Creating backup: {backup_path.name}")
                    gdf.to_file(backup_path, driver='GeoJSON')
                
                # Save fixed file
                print(f"   💾 Saving fixed geometries...")
                fixed_gdf.to_file(file_path, driver='GeoJSON')
                
                print(f"   ✅ Fixed {fixed_count} geometries, removed {removed_count} features")
                print(f"   📊 Final file: {len(fixed_gdf)} features, all valid")
                
                return {
                    'status': 'fixed', 
                    'file': file_path.name,
                    'original_features': len(gdf),
                    'final_features': len(fixed_gdf),
                    'fixed_count': fixed_count,
                    'removed_count': removed_count
                }
            else:
                print(f"   ❌ Still have {final_invalid} invalid geometries after fixing")
                return {
                    'status': 'failed', 
                    'file': file_path.name,
                    'remaining_invalid': final_invalid
                }
        else:
            print(f"   ❌ No valid geometries remaining after fixes")
            return {'status': 'failed', 'file': file_path.name, 'error': 'No valid geometries'}
            
    except Exception as e:
        print(f"   ❌ Error processing file: {e}")
        return {'status': 'error', 'file': file_path.name, 'error': str(e)}

def main():
    """Check and fix all bioregion files"""
    
    print("=== Bioregion Geometry Repair Tool ===")
    print("Checking for invalid geometries in bioregion files...")
    
    bioregions_dir = Path("outputs/bioregions")
    results = []
    
    for filename in BIOREGION_FILES:
        file_path = bioregions_dir / filename
        result = check_and_fix_bioregion_file(file_path)
        results.append(result)
    
    # Summary report
    print("\n" + "="*60)
    print("GEOMETRY REPAIR SUMMARY")
    print("="*60)
    
    valid_files = [r for r in results if r['status'] == 'valid']
    fixed_files = [r for r in results if r['status'] == 'fixed']
    failed_files = [r for r in results if r['status'] in ['failed', 'error']]
    not_found_files = [r for r in results if r['status'] == 'not_found']
    
    print(f"✅ Files already valid: {len(valid_files)}")
    for r in valid_files:
        print(f"   - {r['file']} ({r['features']} features)")
    
    print(f"🔧 Files repaired: {len(fixed_files)}")
    for r in fixed_files:
        print(f"   - {r['file']}: {r['original_features']} → {r['final_features']} features "
              f"(fixed {r['fixed_count']}, removed {r['removed_count']})")
    
    if failed_files:
        print(f"❌ Files with issues: {len(failed_files)}")
        for r in failed_files:
            print(f"   - {r['file']}: {r.get('error', 'Failed to fix')}")
    
    if not_found_files:
        print(f"⚠️  Files not found: {len(not_found_files)}")
        for r in not_found_files:
            print(f"   - {r['file']}")
    
    # Save summary
    summary_path = bioregions_dir / "geometry_repair_summary.json"
    with open(summary_path, 'w') as f:
        json.dump({
            'timestamp': pd.Timestamp.now().isoformat(),
            'results': results,
            'summary': {
                'valid_files': len(valid_files),
                'fixed_files': len(fixed_files),  
                'failed_files': len(failed_files),
                'not_found_files': len(not_found_files)
            }
        }, f, indent=2)
    
    print(f"\n💾 Detailed summary saved to: {summary_path}")
    
    if len(fixed_files) > 0 or len(valid_files) == len(BIOREGION_FILES):
        print("\n🎉 Geometry repair completed successfully!")
        print("You can now try running combine_bioregions.py again.")
    else:
        print(f"\n⚠️  Some files still have issues. Check the summary above.")
    
    return results

if __name__ == "__main__":
    main()