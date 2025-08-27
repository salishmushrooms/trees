#!/usr/bin/env python3
"""
Quick test for bioregion combination to identify issues
"""

import geopandas as gpd
from pathlib import Path

# Test files
TEST_FILES = [
    'eastern_cascades_constrained.geojson',
    'western_cascades_constrained.geojson', 
    'urban_range_constrained.geojson'
]

def test_file_loading():
    print("=== Testing Bioregion File Loading ===")
    
    bioregions_dir = Path('outputs/bioregions')
    
    for filename in TEST_FILES:
        file_path = bioregions_dir / filename
        
        if not file_path.exists():
            print(f"❌ {filename} - NOT FOUND")
            continue
            
        try:
            print(f"📍 Testing {filename}...")
            gdf = gpd.read_file(file_path)
            print(f"   ✅ Loaded: {len(gdf)} features")
            
            # Check for invalid geometries
            invalid_count = (~gdf.geometry.is_valid).sum()
            if invalid_count > 0:
                print(f"   ⚠️  {invalid_count} invalid geometries")
            else:
                print(f"   ✅ All geometries valid")
                
            # Test simplification
            print(f"   🔧 Testing simplification...")
            simplified = gdf.geometry.simplify(0.001, preserve_topology=True)
            invalid_after_simplify = (~simplified.is_valid).sum()
            
            if invalid_after_simplify > 0:
                print(f"   ⚠️  {invalid_after_simplify} invalid after simplification")
            else:
                print(f"   ✅ Simplification OK")
                
        except Exception as e:
            print(f"   ❌ Error loading {filename}: {e}")

if __name__ == "__main__":
    test_file_loading()