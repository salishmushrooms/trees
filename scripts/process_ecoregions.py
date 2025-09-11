#!/usr/bin/env python3
"""
Process Oregon and Washington Level 4 ecoregions for Mapbox upload.
Merges shapefiles, simplifies geometry, and exports as GeoJSON.
"""

import geopandas as gpd
import pandas as pd
from pathlib import Path

def process_ecoregions():
    # Load ecoregion shapefiles
    print("Loading ecoregion shapefiles...")
    or_eco = gpd.read_file('data/ecoregions/or_eco_l4/or_eco_l4.shp')
    wa_eco = gpd.read_file('data/ecoregions/wa_eco_l4/wa_eco_l4.shp')
    
    print(f"Oregon: {len(or_eco)} polygons")
    print(f"Washington: {len(wa_eco)} polygons")
    
    # Merge datasets
    print("Merging datasets...")
    combined = pd.concat([or_eco, wa_eco], ignore_index=True)
    combined = gpd.GeoDataFrame(combined)
    
    # Select essential columns
    keep_columns = [
        'US_L4CODE', 'US_L4NAME', 'US_L3CODE', 'US_L3NAME', 
        'STATE_NAME', 'EPA_REGION', 'geometry'
    ]
    combined = combined[keep_columns]
    
    # Transform to WGS84 for web compatibility
    print("Transforming to WGS84...")
    combined = combined.to_crs('EPSG:4326')
    
    # Simplify geometry for web performance (10m tolerance)
    print("Simplifying geometry...")
    # Convert tolerance from meters to degrees (rough approximation)
    tolerance_deg = 10 / 111000  # ~10m in degrees
    combined['geometry'] = combined['geometry'].simplify(tolerance_deg, preserve_topology=True)
    
    # Create output directory
    output_dir = Path('outputs/ecoregions')
    output_dir.mkdir(exist_ok=True)
    
    # Export as GeoJSON
    output_file = output_dir / 'pnw_ecoregions_l4.geojson'
    print(f"Exporting to {output_file}...")
    combined.to_file(output_file, driver='GeoJSON')
    
    print(f"Processing complete!")
    print(f"Total polygons: {len(combined)}")
    print(f"Output file: {output_file}")
    print(f"File size: {output_file.stat().st_size / 1024 / 1024:.1f} MB")
    
    return str(output_file)

if __name__ == "__main__":
    process_ecoregions()