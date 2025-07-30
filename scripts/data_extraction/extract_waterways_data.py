#!/usr/bin/env python3
"""
Extract Waterways Data for Mushroom Habitat Modeling

This script downloads and processes waterways data from the National Hydrography Dataset (NHD)
for Washington State, calculates distance-to-water metrics for FIA plots, 
and analyzes mushroom observations in relation to water proximity (this is not needed).

Key outputs:
- Waterways shapefile for Washington State
- Distance-to-water calculations for all FIA plots
- Water proximity analysis for mushroom observations (not needed)
- Riparian zone classifications

"""

import geopandas as gpd
import pandas as pd
import requests
import os
from pathlib import Path
import numpy as np
from shapely.geometry import Point
from scipy.spatial import cKDTree
import sqlite3

# Set up paths
DATA_DIR = Path("../../data")
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"
OUTPUT_DIR = Path("../../outputs")

# Create directories if they don't exist
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

def download_nhd_data():
    """
    Download National Hydrography Dataset for Washington State
    
    Note: This is a placeholder function. In practice, you would:
    1. Access USGS National Map or similar service
    2. Download WA state hydrography data
    3. Or use web services like WFS to query specific regions
    """
    print("Downloading NHD data for Washington State...")
    print("Note: This requires access to USGS National Map or similar service")
    
    # For now, we'll use a placeholder approach
    # In real implementation, you would download from:
    # https://www.usgs.gov/core-science-systems/ngp/national-hydrography
    
    # Example web service approach (would need actual endpoints):
    # wa_bounds = get_wa_state_bounds()
    # nhd_url = f"https://hydro.nationalmap.gov/arcgis/rest/services/nhd/MapServer/WFSServer"
    # waterways = gpd.read_file(nhd_url, bbox=wa_bounds)
    
    print("TODO: Implement actual NHD data download")
    print("For now, proceeding with placeholder data structure...")
    
    return None

def calculate_water_proximity(plots_gdf, waterways_gdf):
    """
    Calculate distance to nearest water feature for each FIA plot
    
    Args:
        plots_gdf: GeoDataFrame of FIA plots
        waterways_gdf: GeoDataFrame of water features
    
    Returns:
        plots_gdf with added water proximity columns
    """
    print("Calculating distance to nearest water features...")
    
    if waterways_gdf is None:
        print("No waterways data available - creating placeholder distance values")
        # Create placeholder distances for demonstration
        np.random.seed(42)  # For reproducible results
        plots_gdf['dist_to_water_m'] = np.random.exponential(scale=1000, size=len(plots_gdf))
        plots_gdf['nearest_water_type'] = np.random.choice(
            ['stream', 'river', 'lake', 'wetland'], 
            size=len(plots_gdf),
            p=[0.6, 0.2, 0.1, 0.1]
        )
        return plots_gdf
    
    # Convert to projected CRS for accurate distance calculations
    plots_proj = plots_gdf.to_crs('EPSG:5070')  # Albers Equal Area
    waterways_proj = waterways_gdf.to_crs('EPSG:5070')
    
    # Create spatial index for efficient distance calculations
    water_points = [geom.centroid for geom in waterways_proj.geometry]
    water_coords = np.array([(point.x, point.y) for point in water_points])
    
    plot_coords = np.array([(geom.x, geom.y) for geom in plots_proj.geometry])
    
    # Use KDTree for fast nearest neighbor search
    tree = cKDTree(water_coords)
    distances, indices = tree.query(plot_coords)
    
    # Add results to plots GeoDataFrame
    plots_gdf['dist_to_water_m'] = distances
    plots_gdf['nearest_water_type'] = waterways_proj.iloc[indices]['feature_type'].values
    
    return plots_gdf

def classify_moisture_zones(plots_gdf):
    """
    Classify plots into moisture zones based on water proximity
    
    Args:
        plots_gdf: GeoDataFrame with distance to water calculations
    
    Returns:
        plots_gdf with added moisture zone classifications
    """
    print("Classifying moisture zones...")
    
    def get_moisture_zone(distance):
        """Classify moisture zones based on distance to water"""
        if distance <= 100:
            return 'Riparian'
        elif distance <= 500:
            return 'Near-water'
        elif distance <= 1500:
            return 'Moderate'
        else:
            return 'Dry'
    
    plots_gdf['moisture_zone'] = plots_gdf['dist_to_water_m'].apply(get_moisture_zone)
    
    # Calculate additional moisture-related metrics
    plots_gdf['water_influence'] = np.exp(-plots_gdf['dist_to_water_m'] / 1000)  # Exponential decay
    plots_gdf['riparian_binary'] = (plots_gdf['dist_to_water_m'] <= 100).astype(int)
    
    return plots_gdf

def analyze_mushroom_water_associations():
    """
    Analyze mushroom observations in relation to water proximity
    """
    print("Analyzing mushroom-water associations...")
    
    # Load mushroom data
    mushroom_file = OUTPUT_DIR / "gbif_mushroom_data_filtered.csv"
    if not mushroom_file.exists():
        print(f"Mushroom data not found at {mushroom_file}")
        return None
    
    mushrooms = pd.read_csv(mushroom_file)
    
    # Convert to GeoDataFrame
    mushroom_gdf = gpd.GeoDataFrame(
        mushrooms,
        geometry=[Point(xy) for xy in zip(mushrooms['longitude'], mushrooms['latitude'])],
        crs='EPSG:4326'
    )
    
    # For now, create placeholder water proximity analysis
    # In real implementation, this would use actual waterways data
    np.random.seed(42)
    mushroom_gdf['dist_to_water_m'] = np.random.exponential(scale=800, size=len(mushroom_gdf))
    mushroom_gdf['moisture_zone'] = mushroom_gdf['dist_to_water_m'].apply(
        lambda x: 'Riparian' if x <= 100 else 'Near-water' if x <= 500 else 'Moderate' if x <= 1500 else 'Dry'
    )
    
    # Analyze species-specific water associations
    species_water_stats = mushroom_gdf.groupby('species_name').agg({
        'dist_to_water_m': ['mean', 'median', 'std'],
        'moisture_zone': lambda x: x.mode().iloc[0] if len(x) > 0 else 'Unknown'
    }).round(1)
    
    print("\nSpecies-Water Proximity Analysis:")
    print(species_water_stats.head(10))
    
    # Seasonal water association analysis
    if 'month' in mushroom_gdf.columns:
        seasonal_water = mushroom_gdf.groupby('month')['dist_to_water_m'].mean()
        print(f"\nSeasonal Water Proximity (mean distance):")
        for month, dist in seasonal_water.items():
            print(f"  Month {month}: {dist:.0f}m")
    
    return mushroom_gdf

def main():
    """Main execution function"""
    print("=== Waterways Data Extraction for Mushroom Habitat Modeling ===\n")
    
    # Step 1: Download/access waterways data
    waterways_gdf = download_nhd_data()
    
    # Step 2: Load FIA plot data
    forest_data_file = OUTPUT_DIR / "comprehensive_forest_data_wa.geojson"
    if not forest_data_file.exists():
        print(f"Error: Forest data not found at {forest_data_file}")
        print("Please run extract_comprehensive_forest_data.py first")
        return
    
    print("Loading FIA forest data...")
    plots_gdf = gpd.read_file(forest_data_file)
    print(f"Loaded {len(plots_gdf)} forest plots")
    
    # Step 3: Calculate water proximity for FIA plots
    plots_with_water = calculate_water_proximity(plots_gdf, waterways_gdf)
    
    # Step 4: Classify moisture zones
    plots_final = classify_moisture_zones(plots_with_water)
    
    # Step 5: Analyze mushroom-water associations
    mushroom_water_analysis = analyze_mushroom_water_associations()
    
    # Step 6: Save results
    output_file = OUTPUT_DIR / "forest_data_with_water_proximity.geojson"
    plots_final.to_file(output_file, driver='GeoJSON')
    print(f"\nResults saved to: {output_file}")
    
    # Step 7: Generate summary statistics
    print("\n=== Water Proximity Summary ===")
    print(f"Distance to water statistics (meters):")
    print(f"  Mean: {plots_final['dist_to_water_m'].mean():.0f}")
    print(f"  Median: {plots_final['dist_to_water_m'].median():.0f}")
    print(f"  Min: {plots_final['dist_to_water_m'].min():.0f}")
    print(f"  Max: {plots_final['dist_to_water_m'].max():.0f}")
    
    print(f"\nMoisture zone distribution:")
    zone_counts = plots_final['moisture_zone'].value_counts()
    for zone, count in zone_counts.items():
        pct = 100 * count / len(plots_final)
        print(f"  {zone}: {count:,} plots ({pct:.1f}%)")
    
    print(f"\nForest type vs. water proximity:")
    forest_water_stats = plots_final.groupby('FOREST_TYPE_PREDICTED')['dist_to_water_m'].agg(['mean', 'median'])
    print(forest_water_stats.round(0))
    
    print("\n=== Next Steps ===")
    print("1. Acquire actual NHD data for Washington State")
    print("2. Implement real waterways data processing")
    print("3. Validate water proximity calculations with field knowledge")
    print("4. Analyze mushroom species-specific water dependencies")
    print("5. Integrate with forest type models for habitat suitability")

if __name__ == "__main__":
    main() 