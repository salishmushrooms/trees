#!/usr/bin/env python3
"""
Optimize large GeoJSON files for better QGIS performance by converting to GeoPackage
with spatial indexing and optional simplification.
"""

import geopandas as gpd
import pandas as pd
from pathlib import Path
import time
import sys

def optimize_geojson_to_geopackage(input_path, output_path=None):
    """Convert GeoJSON to optimized GeoPackage with spatial index."""
    
    if output_path is None:
        output_path = Path(input_path).with_suffix('.gpkg')
    
    print(f"Reading GeoJSON file: {input_path}")
    start_time = time.time()
    
    # Read the GeoJSON file
    gdf = gpd.read_file(input_path)
    
    print(f"Loaded {len(gdf)} features in {time.time() - start_time:.2f} seconds")
    print(f"CRS: {gdf.crs}")
    
    # Save as GeoPackage
    print(f"\nWriting to GeoPackage: {output_path}")
    start_time = time.time()
    
    gdf.to_file(output_path, driver='GPKG', layer='plot_carbon_percentiles')
    
    print(f"Saved in {time.time() - start_time:.2f} seconds")
    
    # Compare file sizes
    input_size = Path(input_path).stat().st_size / (1024 * 1024)
    output_size = Path(output_path).stat().st_size / (1024 * 1024)
    
    print(f"\nFile size comparison:")
    print(f"  GeoJSON: {input_size:.1f} MB")
    print(f"  GeoPackage: {output_size:.1f} MB")
    print(f"  Size reduction: {(1 - output_size/input_size) * 100:.1f}%")
    
    return output_path

def create_filtered_version(input_path, min_carbon_percentile=10):
    """Create a filtered version with only significant carbon plots."""
    
    print(f"\nCreating filtered version (min percentile: {min_carbon_percentile})")
    
    gdf = gpd.read_file(input_path)
    
    # Find all carbon percentile columns
    percentile_cols = [col for col in gdf.columns if col.endswith('_PERCENTILE') and not col.endswith('_PERCENTILE_x') and not col.endswith('_PERCENTILE_y')]
    
    # Create mask for any species with significant carbon
    mask = pd.Series(False, index=gdf.index)
    for col in percentile_cols:
        mask |= (gdf[col] >= min_carbon_percentile)
    
    filtered_gdf = gdf[mask].copy()
    
    output_path = Path(input_path).parent / f"{Path(input_path).stem}_filtered_{min_carbon_percentile}pct.gpkg"
    filtered_gdf.to_file(output_path, driver='GPKG', layer='filtered_plots')
    
    print(f"Filtered from {len(gdf)} to {len(filtered_gdf)} features ({len(filtered_gdf)/len(gdf)*100:.1f}%)")
    print(f"Saved to: {output_path}")
    
    return output_path

def create_species_specific_layers(input_path):
    """Create separate layers for each species to improve query performance."""
    
    print("\nCreating species-specific layers...")
    
    gdf = gpd.read_file(input_path)
    
    # Find all species columns
    carbon_cols = [col for col in gdf.columns if col.endswith('_AG_CARBON') and not col.endswith('_AG_CARBON_x') and not col.endswith('_AG_CARBON_y')]
    species_names = [col.replace('_AG_CARBON', '') for col in carbon_cols]
    
    output_path = Path(input_path).parent / f"{Path(input_path).stem}_by_species.gpkg"
    
    # Save base layer with common attributes
    base_cols = ['PLOT_CN', 'STATE', 'LAT', 'LON', 'ELEVATION_FT', 'TOTAL_TREES', 'TOTAL_CARBON_AG', 'geometry']
    gdf[base_cols].to_file(output_path, driver='GPKG', layer='base_plots')
    
    # Save species-specific layers
    for species in species_names[:5]:  # Limit to top 5 species for demo
        carbon_col = f"{species}_AG_CARBON"
        percentile_col = f"{species}_AG_CARBON_PERCENTILE"
        
        # Filter to plots with this species
        species_gdf = gdf[gdf[carbon_col] > 0][[*base_cols, carbon_col, percentile_col]].copy()
        
        if len(species_gdf) > 0:
            layer_name = species.lower().replace('_', '-')
            species_gdf.to_file(output_path, driver='GPKG', layer=layer_name, mode='a')
            print(f"  {layer_name}: {len(species_gdf)} plots")
    
    print(f"\nSaved multi-layer GeoPackage to: {output_path}")
    return output_path

if __name__ == "__main__":
    input_file = "outputs/plot_carbon_percentiles_latest_surveys.geojson"
    
    # Convert to GeoPackage
    gpkg_file = optimize_geojson_to_geopackage(input_file)
    
    # Create filtered version
    filtered_file = create_filtered_version(input_file, min_carbon_percentile=20)
    
    # Create species-specific layers
    species_file = create_species_specific_layers(input_file)
    
    print("\n" + "="*60)
    print("QGIS Performance Tips:")
    print("1. Use the .gpkg files instead of .geojson")
    print("2. Enable spatial indexing in QGIS layer properties")
    print("3. Use the filtered version for quick overview")
    print("4. Load only specific species layers as needed")
    print("5. Set scale-dependent rendering (hide details when zoomed out)")