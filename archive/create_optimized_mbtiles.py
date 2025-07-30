#!/usr/bin/env python3
"""
Create optimized MBTiles directly from source GeoJSON files
Combines optimization and MBTiles creation in one step
"""

import subprocess
import tempfile
from pathlib import Path
import argparse
import json
import geopandas as gpd
from shapely.geometry import mapping

def round_coordinates(geom, precision=6):
    """Round coordinates to specified decimal places"""
    def round_coords(coords):
        return [round(float(coord), precision) for coord in coords]
    
    def round_nested_coords(coords_list):
        if not coords_list:
            return coords_list
        
        if isinstance(coords_list[0], (list, tuple)):
            return [round_nested_coords(sublist) for sublist in coords_list]
        else:
            return round_coords(coords_list)
    
    geom_dict = mapping(geom)
    coords = geom_dict['coordinates']
    
    if geom_dict['type'] == 'Point':
        geom_dict['coordinates'] = round_coords(coords)
    elif geom_dict['type'] in ['LineString', 'MultiPoint']:
        geom_dict['coordinates'] = round_nested_coords(coords)
    elif geom_dict['type'] in ['Polygon', 'MultiLineString']:
        geom_dict['coordinates'] = [round_nested_coords(ring) for ring in coords]
    elif geom_dict['type'] == 'MultiPolygon':
        geom_dict['coordinates'] = [[round_nested_coords(ring) for ring in polygon] for polygon in coords]
    
    return geom_dict

def get_essential_properties(properties):
    """Keep only essential properties for web mapping"""
    return {
        'species': properties.get('species'),
        'species_code': properties.get('species_code'),
        'total_plots': properties.get('total_plots')
    }

def optimize_geojson_for_mbtiles(input_file, precision=6):
    """Optimize GeoJSON file and return as temporary file"""
    
    # Read and optimize
    gdf = gpd.read_file(input_file)
    
    optimized_features = []
    for idx, row in gdf.iterrows():
        feature = {
            "type": "Feature", 
            "geometry": round_coordinates(row.geometry, precision),
            "properties": get_essential_properties(dict(row.drop('geometry')))
        }
        optimized_features.append(feature)
    
    optimized_geojson = {
        "type": "FeatureCollection",
        "features": optimized_features
    }
    
    # Save to temporary file
    temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.geojson', delete=False)
    json.dump(optimized_geojson, temp_file, separators=(',', ':'))
    temp_file.close()
    
    return temp_file.name

def extract_species_name(filename):
    """Extract species name from filename"""
    name = filename.replace('_continuous_distribution_pnw_elevation_tcc_constrained_from_edited_cached', '')
    name = name.replace('.geojson', '')
    return name.replace('_', ' ').replace('-', ' ').title()

def create_optimized_mbtiles(input_dir, output_file, pattern="*continuous_distribution*tcc_constrained*cached.geojson"):
    """Create MBTiles with on-the-fly optimization"""
    
    input_path = Path(input_dir)
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Find matching files
    geojson_files = list(input_path.glob(pattern))
    if not geojson_files:
        print(f"❌ No files found matching: {pattern}")
        return False
    
    print(f"🔍 Found {len(geojson_files)} files to process")
    
    # Optimize files to temporary directory
    temp_files = []
    original_size = 0
    
    try:
        for geojson_file in sorted(geojson_files):
            print(f"🔄 Optimizing: {geojson_file.name}")
            
            original_size += geojson_file.stat().st_size / (1024 * 1024)
            temp_file = optimize_geojson_for_mbtiles(geojson_file)
            temp_files.append((temp_file, geojson_file.name))
        
        print(f"📁 Total original size: {original_size:.1f} MB")
        
        # Build tippecanoe command
        cmd = [
            'tippecanoe',
            '--output', str(output_path),
            '--minimum-zoom', '4',
            '--maximum-zoom', '12',
            '--simplification', '2.0',
            '--force',
            '--no-feature-limit',
            '--no-tile-size-limit',
            '--drop-densest-as-needed',
            '--extend-zooms-if-still-dropping',
            '--name', 'PNW Species Distribution Optimized',
            '--description', 'Optimized tree species distribution data for Pacific Northwest',
            '--quiet'
        ]
        
        # Add each optimized file as a layer
        for temp_file, original_name in temp_files:
            species_name = extract_species_name(Path(original_name).stem)
            layer_name = f"{species_name.lower().replace(' ', '_').replace('-', '_')}_range"
            cmd.extend(['--layer', f"{layer_name}:{temp_file}"])
        
        print(f"🔄 Creating MBTiles with tippecanoe...")
        result = subprocess.run(cmd, check=True)
        
        if output_path.exists():
            output_size = output_path.stat().st_size / (1024 * 1024)
            print(f"✅ MBTiles created: {output_size:.1f} MB")
            print(f"📉 Size reduction: {((original_size - output_size) / original_size * 100):.1f}%")
            return True
        else:
            print("❌ MBTiles creation failed")
            return False
    
    finally:
        # Clean up temporary files
        for temp_file, _ in temp_files:
            try:
                Path(temp_file).unlink()
            except:
                pass

def main():
    parser = argparse.ArgumentParser(description='Create optimized MBTiles from GeoJSON files')
    parser.add_argument('input_dir', help='Directory containing GeoJSON files')
    parser.add_argument('-o', '--output', default='outputs/mapbox_mbtiles/pnw_species_optimized.mbtiles',
                       help='Output MBTiles file')
    parser.add_argument('--pattern', default='*continuous_distribution*tcc_constrained*cached.geojson',
                       help='File pattern to match')
    
    args = parser.parse_args()
    
    success = create_optimized_mbtiles(args.input_dir, args.output, args.pattern)
    
    if success:
        print(f"\n💡 Next steps:")
        print(f"1. Upload to Mapbox Studio: {args.output}")
        print(f"2. All species available as separate layers")
    else:
        exit(1)

if __name__ == "__main__":
    main()