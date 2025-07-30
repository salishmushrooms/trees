#!/usr/bin/env python3
"""
Extract Coastal Forest Species Data

Extracts plot data for Sitka spruce and lodgepole pine in coastal areas
to support coastal forest bioregion mapping.
"""

import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import json
from pathlib import Path

# Configuration
SITKA_SPRUCE_CODE = 98
LODGEPOLE_PINE_CODE = 108
MAX_COASTAL_ELEVATION_FT = 1000
COASTAL_LONGITUDE_THRESHOLD = -123.0  # West of this line

# Database paths
DB_PATHS = {
    'WA': 'data/raw/trees_SQLite_FIADB_WA.db',
    'OR': 'data/raw/trees_SQLite_FIADB_OR.db',
}

def extract_coastal_species_from_csv():
    """Extract from existing comprehensive forest data files"""
    print("Extracting coastal species from existing data files...")
    
    # First check if we have comprehensive forest data
    wa_data_path = Path("outputs/geojson/comprehensive_forest_data_wa.geojson")
    
    coastal_plots = []
    
    if wa_data_path.exists():
        print(f"Loading data from {wa_data_path}")
        gdf = gpd.read_file(wa_data_path)
        
        # Filter for coastal areas
        coastal_gdf = gdf[
            (gdf['LON'] < COASTAL_LONGITUDE_THRESHOLD) &
            (gdf['ELEV'] <= MAX_COASTAL_ELEVATION_FT)
        ]
        
        # Check for Sitka spruce columns
        sitka_cols = [col for col in coastal_gdf.columns if 'SITKA' in col.upper()]
        if sitka_cols:
            sitka_plots = coastal_gdf[coastal_gdf[sitka_cols[0]] > 0]
            print(f"Found {len(sitka_plots)} plots with Sitka spruce")
            coastal_plots.append(sitka_plots)
        
        # For lodgepole pine, we need to be more careful as it includes inland varieties
        # Filter by both species presence AND coastal location
        pine_cols = [col for col in coastal_gdf.columns if 'LODGEPOLE' in col.upper() or ('PINE' in col.upper() and '108' in str(coastal_gdf.columns))]
        
        print(f"Available columns: {list(coastal_gdf.columns[:20])}...")  # Show first 20 columns
        
    # Also check plot carbon data which might have different columns
    plot_carbon_path = Path("outputs/plot_carbon_percentiles_latest_surveys.geojson")
    if plot_carbon_path.exists():
        print(f"\nChecking plot carbon data at {plot_carbon_path}")
        carbon_gdf = gpd.read_file(plot_carbon_path)
        
        # Filter for coastal
        coastal_carbon = carbon_gdf[
            (carbon_gdf.geometry.x < COASTAL_LONGITUDE_THRESHOLD) &
            (carbon_gdf['ELEVATION_FT'] <= MAX_COASTAL_ELEVATION_FT)
        ]
        
        print(f"Found {len(coastal_carbon)} coastal low-elevation plots")
        
        # Check for Sitka spruce presence
        if 'SITKA_SPRUCE_AG_CARBON' in coastal_carbon.columns:
            sitka_plots = coastal_carbon[coastal_carbon['SITKA_SPRUCE_AG_CARBON'] > 0]
            print(f"Found {len(sitka_plots)} coastal plots with Sitka spruce")
        
        # Check for lodgepole pine presence (shore pine at coast)
        if 'LODGEPOLE_PINE_AG_CARBON' in coastal_carbon.columns:
            pine_plots = coastal_carbon[coastal_carbon['LODGEPOLE_PINE_AG_CARBON'] > 0]
            print(f"Found {len(pine_plots)} coastal plots with lodgepole/shore pine")
        
        # Create output with coastal indicator
        coastal_carbon['is_coastal'] = True
        coastal_carbon['coastal_distance_deg'] = COASTAL_LONGITUDE_THRESHOLD - coastal_carbon.geometry.x
        
        return coastal_carbon
    
    return None

def create_coastal_species_summary(coastal_data):
    """Create summary statistics for coastal species"""
    if coastal_data is None or len(coastal_data) == 0:
        return None
    
    # Handle different elevation column names
    elev_col = 'ELEVATION_FT' if 'ELEVATION_FT' in coastal_data.columns else 'ELEV'
    
    summary = {
        'total_coastal_plots': len(coastal_data),
        'elevation_range': {
            'min': coastal_data[elev_col].min(),
            'max': coastal_data[elev_col].max(),
            'mean': coastal_data[elev_col].mean()
        },
        'latitude_range': {
            'min': coastal_data.geometry.y.min(),
            'max': coastal_data.geometry.y.max()
        },
        'longitude_range': {
            'min': coastal_data.geometry.x.min(),
            'max': coastal_data.geometry.x.max()
        }
    }
    
    # Add forest type distribution if available
    if 'FORTYPCD' in coastal_data.columns:
        forest_types = coastal_data['FORTYPCD'].value_counts().to_dict()
        summary['forest_types'] = forest_types
    
    return summary

def main():
    """Main extraction function"""
    print("=== Extracting Coastal Forest Species Data ===\n")
    
    # Create output directory
    output_dir = Path("outputs/coastal_analysis")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Extract coastal data
    coastal_data = extract_coastal_species_from_csv()
    
    if coastal_data is not None and len(coastal_data) > 0:
        # Save coastal plots
        output_file = output_dir / "coastal_forest_plots.geojson"
        coastal_data.to_file(output_file, driver='GeoJSON')
        print(f"\nSaved {len(coastal_data)} coastal plots to {output_file}")
        
        # Create and save summary
        summary = create_coastal_species_summary(coastal_data)
        if summary:
            summary_file = output_dir / "coastal_forest_summary.json"
            with open(summary_file, 'w') as f:
                json.dump(summary, f, indent=2)
            print(f"Saved summary to {summary_file}")
            
            # Print summary
            print("\n=== Coastal Forest Summary ===")
            print(f"Total coastal plots: {summary['total_coastal_plots']}")
            print(f"Elevation range: {summary['elevation_range']['min']:.0f} - {summary['elevation_range']['max']:.0f} ft")
            print(f"Latitude range: {summary['latitude_range']['min']:.2f} - {summary['latitude_range']['max']:.2f}")
            print(f"Longitude range: {summary['longitude_range']['min']:.2f} - {summary['longitude_range']['max']:.2f}")
    else:
        print("No coastal data found. You may need to run extract_comprehensive_forest_data.py first")
        print("Or adjust the coastal criteria (elevation/longitude thresholds)")

if __name__ == "__main__":
    main()