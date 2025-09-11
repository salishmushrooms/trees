#!/usr/bin/env python3

import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import matplotlib.pyplot as plt
import numpy as np

def analyze_chanterelle_ecoregions():
    """
    Analyze chanterelle observations by ecoregions, filtering out obscured locations.
    Creates summary charts grouped by species and ecoregion codes (L3 and L4).
    """
    
    # Read chanterelle data
    print("Loading chanterelle observations...")
    chanterelles = pd.read_csv('/Users/JJC/trees/data/mushrooms/chanterelles.csv')
    
    print(f"Total observations: {len(chanterelles)}")
    print(f"Columns: {list(chanterelles.columns)}")
    
    # Check geoprivacy values
    print(f"\nGeoprivacy values:")
    print(chanterelles['geoprivacy'].value_counts(dropna=False))
    
    # Filter out obscured observations (where geoprivacy is not null/empty)
    # Keep only observations with accurate coordinates (no geoprivacy restrictions)
    accurate_obs = chanterelles[chanterelles['geoprivacy'].isna() | (chanterelles['geoprivacy'] == '')]
    
    print(f"\nObservations with accurate coordinates: {len(accurate_obs)}")
    print(f"Filtered out {len(chanterelles) - len(accurate_obs)} obscured observations")
    
    # Create GeoDataFrame from accurate observations
    print("\nCreating spatial points...")
    geometry = [Point(xy) for xy in zip(accurate_obs['longitude'], accurate_obs['latitude'])]
    gdf_chanterelles = gpd.GeoDataFrame(accurate_obs, geometry=geometry, crs='EPSG:4326')
    
    # Load Washington ecoregions shapefile
    print("Loading Washington ecoregions...")
    ecoregions = gpd.read_file('/Users/JJC/trees/data/ecoregions/wa_eco_l4/wa_eco_l4.shp')
    
    print(f"Ecoregions columns: {list(ecoregions.columns)}")
    print(f"Number of ecoregions: {len(ecoregions)}")
    
    # Ensure same CRS for spatial join
    ecoregions = ecoregions.to_crs('EPSG:4326')
    
    # Spatial join to assign ecoregion codes to each observation
    print("Performing spatial join...")
    joined = gpd.sjoin(gdf_chanterelles, ecoregions, how='left', predicate='within')
    
    # Check how many observations got matched
    matched = joined.dropna(subset=['US_L4CODE'])
    print(f"Observations matched to ecoregions: {len(matched)} / {len(gdf_chanterelles)}")
    
    # Create summary by species and L4 ecoregions
    print("\nCreating L4 summary...")
    l4_summary = matched.groupby(['scientific_name', 'US_L4CODE', 'US_L4NAME']).size().reset_index(name='count')
    l4_summary = l4_summary.sort_values(['scientific_name', 'count'], ascending=[True, False])
    
    # Create summary by species and L3 ecoregions  
    print("Creating L3 summary...")
    l3_summary = matched.groupby(['scientific_name', 'US_L3CODE', 'US_L3NAME']).size().reset_index(name='count')
    l3_summary = l3_summary.sort_values(['scientific_name', 'count'], ascending=[True, False])
    
    # Display summaries
    print("\n" + "="*80)
    print("SUMMARY BY SPECIES AND L4 ECOREGIONS")
    print("="*80)
    for species in l4_summary['scientific_name'].unique():
        species_data = l4_summary[l4_summary['scientific_name'] == species]
        print(f"\n{species}:")
        for _, row in species_data.iterrows():
            print(f"  L4 {row['US_L4CODE']}: {row['US_L4NAME']} - {row['count']} observations")
    
    print("\n" + "="*80)
    print("SUMMARY BY SPECIES AND L3 ECOREGIONS")
    print("="*80)
    for species in l3_summary['scientific_name'].unique():
        species_data = l3_summary[l3_summary['scientific_name'] == species]
        print(f"\n{species}:")
        for _, row in species_data.iterrows():
            print(f"  L3 {row['US_L3CODE']}: {row['US_L3NAME']} - {row['count']} observations")
    
    # Create simple bar charts for visualization
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    # L4 ecoregions chart
    for i, species in enumerate(l4_summary['scientific_name'].unique()):
        species_data = l4_summary[l4_summary['scientific_name'] == species]
        ax1.bar([f"{code}" for code in species_data['US_L4CODE']], 
                species_data['count'], 
                label=species.replace('Cantharellus ', ''), 
                alpha=0.7)
    
    ax1.set_title('Chanterelle Observations by L4 Ecoregion Code')
    ax1.set_xlabel('US L4 Ecoregion Code')
    ax1.set_ylabel('Number of Observations')
    ax1.legend()
    ax1.grid(axis='y', alpha=0.3)
    
    # L3 ecoregions chart
    for i, species in enumerate(l3_summary['scientific_name'].unique()):
        species_data = l3_summary[l3_summary['scientific_name'] == species]
        ax2.bar([f"{code}" for code in species_data['US_L3CODE']], 
                species_data['count'], 
                label=species.replace('Cantharellus ', ''), 
                alpha=0.7)
    
    ax2.set_title('Chanterelle Observations by L3 Ecoregion Code')
    ax2.set_xlabel('US L3 Ecoregion Code')
    ax2.set_ylabel('Number of Observations')
    ax2.legend()
    ax2.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('/Users/JJC/trees/outputs/chanterelle_ecoregion_summary.png', dpi=300, bbox_inches='tight')
    print(f"\nVisualization saved to: /Users/JJC/trees/outputs/chanterelle_ecoregion_summary.png")
    
    # Save summary tables
    l4_summary.to_csv('/Users/JJC/trees/outputs/chanterelle_l4_ecoregion_summary.csv', index=False)
    l3_summary.to_csv('/Users/JJC/trees/outputs/chanterelle_l3_ecoregion_summary.csv', index=False)
    
    print(f"L4 summary saved to: /Users/JJC/trees/outputs/chanterelle_l4_ecoregion_summary.csv")
    print(f"L3 summary saved to: /Users/JJC/trees/outputs/chanterelle_l3_ecoregion_summary.csv")
    
    return l4_summary, l3_summary, matched

if __name__ == "__main__":
    l4_summary, l3_summary, matched_data = analyze_chanterelle_ecoregions()