#!/usr/bin/env python3
"""
Pre-process large bioregion files to reduce size for web mapping

This script applies the same simplification and precision reduction that
should have been applied during creation but may have been missed.
"""

import geopandas as gpd
import pandas as pd
from pathlib import Path
from shapely.ops import transform
from shapely.validation import make_valid
import warnings

warnings.filterwarnings('ignore')

# Target files (over 5MB)
LARGE_FILES = [
    'high_cascades_constrained.geojson',
    'urban_range_constrained.geojson', 
    'eastern_cascades_constrained.geojson',
    'western_cascades_constrained.geojson'
]

def preprocess_bioregion_file(file_path):
    """Pre-process a single bioregion file"""
    
    print(f"\n📍 Processing {file_path.name}...")
    
    # Get original file size
    original_size_mb = file_path.stat().st_size / (1024 * 1024)
    print(f"   📏 Original size: {original_size_mb:.1f} MB")
    
    # Load file
    print(f"   📖 Loading file...")
    gdf = gpd.read_file(file_path)
    print(f"   ✅ Loaded {len(gdf)} features")
    
    # Create backup
    backup_path = file_path.with_suffix('.backup.geojson')
    if not backup_path.exists():
        print(f"   💾 Creating backup...")
        gdf.to_file(backup_path, driver='GeoJSON')
    
    # Apply simplification and precision reduction
    print(f"   🔧 Applying geometry simplification...")
    
    # Simplify geometry (remove unnecessary vertices) - more aggressive for large files
    tolerance_degrees = 0.001  # ~100m at mid-latitudes
    gdf['geometry'] = gdf.geometry.simplify(tolerance_degrees, preserve_topology=True)
    
    # Round coordinates to 3 decimal places (~100m precision)
    def round_coordinates(geom):
        """Round coordinates to 3 decimal places (~100m precision)"""
        if geom is None or geom.is_empty:
            return None
        def round_coords(x, y, z=None):
            rounded_x = round(x, 3)
            rounded_y = round(y, 3)
            return (rounded_x, rounded_y) if z is None else (rounded_x, rounded_y, round(z, 3))
        return transform(round_coords, geom)
    
    gdf['geometry'] = gdf['geometry'].apply(round_coordinates)
    
    # Fix any invalid geometries
    invalid_count = (~gdf.geometry.is_valid).sum()
    if invalid_count > 0:
        print(f"   🔧 Fixing {invalid_count} invalid geometries...")
        
        def fix_geometry(geom):
            if geom is None or geom.is_empty:
                return None
            if geom.is_valid:
                return geom
            try:
                return make_valid(geom)
            except:
                return None
        
        gdf['geometry'] = gdf['geometry'].apply(fix_geometry)
        gdf = gdf[gdf.geometry.notna() & ~gdf.geometry.is_empty]
        
        fixed_invalid = (~gdf.geometry.is_valid).sum()
        print(f"   ✅ {invalid_count - fixed_invalid} geometries fixed")
    
    # Remove any duplicate or unnecessary columns for size reduction
    essential_columns = ['geometry', 'region_name', 'region_code', 'area_sqkm', 'total_area_sqkm']
    columns_to_keep = [col for col in essential_columns if col in gdf.columns] + ['geometry']
    
    # Keep some original attributes but remove verbose ones
    for col in gdf.columns:
        if col not in columns_to_keep and col not in ['geometry']:
            # Keep important metadata but drop verbose fields
            if any(keyword in col.lower() for keyword in ['name', 'type', 'date', 'elevation', 'tree_cover']):
                if col not in columns_to_keep:
                    columns_to_keep.append(col)
    
    gdf = gdf[list(set(columns_to_keep))]
    print(f"   📊 Reduced to {len(gdf.columns)} columns: {list(gdf.columns)}")
    
    # Save processed file
    print(f"   💾 Saving processed file...")
    gdf.to_file(file_path, driver='GeoJSON')
    
    # Check new file size
    new_size_mb = file_path.stat().st_size / (1024 * 1024)
    reduction_pct = (original_size_mb - new_size_mb) / original_size_mb * 100
    
    print(f"   ✅ New size: {new_size_mb:.1f} MB ({reduction_pct:.1f}% reduction)")
    print(f"   📊 Final: {len(gdf)} features")
    
    return {
        'file': file_path.name,
        'original_size_mb': original_size_mb,
        'new_size_mb': new_size_mb,
        'reduction_pct': reduction_pct,
        'features': len(gdf)
    }

def main():
    """Process all large bioregion files"""
    
    print("=== Large Bioregion File Preprocessor ===")
    print("Applying geometry simplification and precision reduction...")
    
    bioregions_dir = Path('outputs/bioregions')
    results = []
    
    for filename in LARGE_FILES:
        file_path = bioregions_dir / filename
        
        if not file_path.exists():
            print(f"⚠️  {filename} not found, skipping...")
            continue
            
        try:
            result = preprocess_bioregion_file(file_path)
            results.append(result)
            
        except Exception as e:
            print(f"❌ Error processing {filename}: {e}")
    
    # Summary
    print(f"\n" + "="*60)
    print("PREPROCESSING SUMMARY")
    print("="*60)
    
    total_original = sum(r['original_size_mb'] for r in results)
    total_new = sum(r['new_size_mb'] for r in results)
    total_reduction = (total_original - total_new) / total_original * 100 if total_original > 0 else 0
    
    for result in results:
        print(f"✅ {result['file']}: {result['original_size_mb']:.1f}MB → {result['new_size_mb']:.1f}MB "
              f"({result['reduction_pct']:.1f}% reduction, {result['features']} features)")
    
    print(f"\n📊 Total: {total_original:.1f}MB → {total_new:.1f}MB ({total_reduction:.1f}% reduction)")
    print(f"\n🎉 Large file preprocessing completed!")
    print("You can now try running combine_bioregions.py again.")

if __name__ == "__main__":
    main()