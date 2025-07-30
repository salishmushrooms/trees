#!/usr/bin/env python3
"""
Convert iNaturalist observations CSV to GeoJSON format.
Preserves all attributes including geoprivacy field.
"""

import csv
import json
import sys
from pathlib import Path

def csv_to_geojson(csv_file_path, output_file_path):
    """
    Convert CSV file to GeoJSON format.
    
    Args:
        csv_file_path (str): Path to input CSV file
        output_file_path (str): Path to output GeoJSON file
    """
    features = []
    skipped_count = 0
    total_count = 0
    
    print(f"Reading CSV file: {csv_file_path}")
    
    with open(csv_file_path, 'r', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile)
        
        for row in reader:
            total_count += 1
            
            # Check if latitude and longitude are present and valid
            lat = row.get('latitude', '').strip()
            lon = row.get('longitude', '').strip()
            
            if not lat or not lon or lat == '' or lon == '':
                skipped_count += 1
                continue
                
            try:
                latitude = float(lat)
                longitude = float(lon)
            except (ValueError, TypeError):
                skipped_count += 1
                continue
            
            # Create properties from all CSV columns except lat/lon
            properties = {}
            for key, value in row.items():
                if key not in ['latitude', 'longitude']:
                    # Convert empty strings to None for cleaner JSON
                    properties[key] = value if value.strip() != '' else None
            
            # Create GeoJSON feature
            feature = {
                "type": "Feature",
                "geometry": {
                    "type": "Point",
                    "coordinates": [longitude, latitude]  # GeoJSON uses [lon, lat] order
                },
                "properties": properties
            }
            
            features.append(feature)
    
    # Create GeoJSON structure
    geojson = {
        "type": "FeatureCollection",
        "features": features
    }
    
    # Write to output file
    print(f"Writing GeoJSON file: {output_file_path}")
    with open(output_file_path, 'w', encoding='utf-8') as outfile:
        json.dump(geojson, outfile, indent=2, ensure_ascii=False)
    
    print(f"Conversion complete!")
    print(f"Total records processed: {total_count}")
    print(f"Records with valid coordinates: {len(features)}")
    print(f"Records skipped (missing/invalid coordinates): {skipped_count}")
    
    return geojson

def main():
    """Main function to handle command line arguments."""
    if len(sys.argv) != 3:
        print("Usage: python csv_to_geojson_converter.py <input_csv> <output_geojson>")
        print("Example: python csv_to_geojson_converter.py observations.csv observations.geojson")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    # Check if input file exists
    if not Path(input_file).exists():
        print(f"Error: Input file '{input_file}' not found.")
        sys.exit(1)
    
    # Create output directory if it doesn't exist
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Convert CSV to GeoJSON
    csv_to_geojson(input_file, output_file)

if __name__ == "__main__":
    main() 