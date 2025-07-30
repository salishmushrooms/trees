#!/usr/bin/env python3
"""
Optimize existing GeoJSON files for MBTiles creation by:
1. Reducing coordinate precision to ~10m accuracy (6 decimal places)
2. Removing unnecessary metadata properties
3. Optional geometry simplification
"""

import json
import os
from pathlib import Path
import argparse
import geopandas as gpd
from shapely.geometry import mapping
import numpy as np

def round_coordinates(geom, precision=6):
    """Round coordinates to specified decimal places"""
    def round_coords(coords):
        return [round(float(coord), precision) for coord in coords]
    
    def round_nested_coords(coords_list):
        if not coords_list:
            return coords_list
        
        # Check if it's a nested list (like polygon with holes)
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
    essential_props = {
        'species': properties.get('species'),
        'species_code': properties.get('species_code'),
        'total_plots': properties.get('total_plots')
    }
    
    # Remove None values
    return {k: v for k, v in essential_props.items() if v is not None}

def optimize_geojson_file(input_file, output_file=None, precision=6, simplify_tolerance=None, 
                         keep_all_properties=False, dry_run=False):
    """Optimize a single GeoJSON file"""
    
    input_path = Path(input_file)
    if output_file is None:
        output_path = input_path.parent / f"{input_path.stem}_optimized.geojson"
    else:
        output_path = Path(output_file)
    
    print(f"\n📁 Processing: {input_path.name}")
    
    # Get original file size
    original_size_mb = input_path.stat().st_size / (1024 * 1024)
    print(f"   Original size: {original_size_mb:.1f} MB")
    
    # Read the GeoJSON
    gdf = gpd.read_file(input_path)
    
    # Optional geometry simplification
    if simplify_tolerance:
        print(f"   🎨 Simplifying geometry (tolerance: {simplify_tolerance})")
        gdf.geometry = gdf.geometry.simplify(simplify_tolerance, preserve_topology=True)
    
    # Process each feature
    optimized_features = []
    for idx, row in gdf.iterrows():
        # Round coordinates
        optimized_geom = round_coordinates(row.geometry, precision)
        
        # Filter properties
        if keep_all_properties:
            properties = dict(row.drop('geometry'))
        else:
            properties = get_essential_properties(dict(row.drop('geometry')))
        
        feature = {
            "type": "Feature",
            "geometry": optimized_geom,
            "properties": properties
        }
        optimized_features.append(feature)
    
    # Create optimized GeoJSON structure
    optimized_geojson = {
        "type": "FeatureCollection",
        "features": optimized_features
    }
    
    if dry_run:
        print(f"   🔍 Would save to: {output_path}")
        return True
    
    # Save optimized file
    with open(output_path, 'w') as f:
        json.dump(optimized_geojson, f, separators=(',', ':'))  # Compact JSON
    
    # Compare file sizes
    optimized_size_mb = output_path.stat().st_size / (1024 * 1024)
    reduction_pct = ((original_size_mb - optimized_size_mb) / original_size_mb) * 100
    
    print(f"   ✅ Optimized size: {optimized_size_mb:.1f} MB")
    print(f"   📉 Size reduction: {reduction_pct:.1f}%")
    
    return True

def batch_optimize_directory(input_dir, output_dir=None, pattern="*continuous_distribution*tcc_constrained*cached.geojson", 
                            precision=6, simplify_tolerance=None, keep_all_properties=False, dry_run=False):
    """Optimize all matching GeoJSON files in a directory"""
    
    input_path = Path(input_dir)
    if output_dir is None:
        output_path = input_path / "optimized"
    else:
        output_path = Path(output_dir)
    
    if not dry_run:
        output_path.mkdir(exist_ok=True)
    
    # Find matching files
    geojson_files = list(input_path.glob(pattern))
    
    if not geojson_files:
        print(f"❌ No files found matching pattern: {pattern}")
        return False
    
    print(f"🔍 Found {len(geojson_files)} files matching pattern")
    print(f"📁 Output directory: {output_path}")
    
    total_original_size = 0
    total_optimized_size = 0
    successful_files = 0
    
    for geojson_file in sorted(geojson_files):
        try:
            # Create output filename
            output_file = output_path / f"{geojson_file.stem}_optimized.geojson"
            
            # Get original size
            original_size = geojson_file.stat().st_size / (1024 * 1024)
            total_original_size += original_size
            
            # Optimize file
            success = optimize_geojson_file(
                geojson_file, 
                output_file, 
                precision=precision,
                simplify_tolerance=simplify_tolerance,
                keep_all_properties=keep_all_properties,
                dry_run=dry_run
            )
            
            if success:
                successful_files += 1
                if not dry_run:
                    optimized_size = output_file.stat().st_size / (1024 * 1024)
                    total_optimized_size += optimized_size
        
        except Exception as e:
            print(f"   ❌ Error processing {geojson_file.name}: {e}")
    
    # Summary
    print(f"\n📊 Batch Processing Summary:")
    print(f"   ✅ Successfully processed: {successful_files}/{len(geojson_files)} files")
    
    if not dry_run and successful_files > 0:
        total_reduction_pct = ((total_original_size - total_optimized_size) / total_original_size) * 100
        print(f"   📁 Total original size: {total_original_size:.1f} MB")
        print(f"   📁 Total optimized size: {total_optimized_size:.1f} MB")
        print(f"   📉 Total size reduction: {total_reduction_pct:.1f}%")
    
    return successful_files > 0

def main():
    parser = argparse.ArgumentParser(description='Optimize GeoJSON files for MBTiles creation')
    parser.add_argument('input', help='Input GeoJSON file or directory')
    parser.add_argument('-o', '--output', help='Output file or directory (default: add _optimized suffix)')
    parser.add_argument('--precision', type=int, default=6, 
                       help='Coordinate precision in decimal places (default: 6 for ~10m accuracy)')
    parser.add_argument('--simplify', type=float, 
                       help='Geometry simplification tolerance in degrees (e.g., 0.0001)')
    parser.add_argument('--keep-all-props', action='store_true',
                       help='Keep all properties instead of filtering to essentials')
    parser.add_argument('--pattern', default='*continuous_distribution*tcc_constrained*cached.geojson',
                       help='File pattern to match when processing directory')
    parser.add_argument('--dry-run', action='store_true',
                       help='Show what would be done without creating files')
    
    args = parser.parse_args()
    
    input_path = Path(args.input)
    
    if input_path.is_file():
        # Single file processing
        success = optimize_geojson_file(
            input_path,
            args.output,
            precision=args.precision,
            simplify_tolerance=args.simplify,
            keep_all_properties=args.keep_all_props,
            dry_run=args.dry_run
        )
        
        if not success:
            print("❌ Failed to optimize file")
            exit(1)
    
    elif input_path.is_dir():
        # Directory processing
        success = batch_optimize_directory(
            input_path,
            args.output,
            pattern=args.pattern,
            precision=args.precision,
            simplify_tolerance=args.simplify,
            keep_all_properties=args.keep_all_props,
            dry_run=args.dry_run
        )
        
        if not success:
            print("❌ Failed to optimize directory")
            exit(1)
    
    else:
        print(f"❌ Input path does not exist: {input_path}")
        exit(1)
    
    if not args.dry_run:
        print(f"\n💡 Next steps:")
        print(f"   1. Test the optimized files with your MBTiles script")
        print(f"   2. If satisfied, replace the original files")
        print(f"   3. Consider updating the source script for future runs")

if __name__ == "__main__":
    main()