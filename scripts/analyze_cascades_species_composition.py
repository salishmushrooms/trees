#!/usr/bin/env python3
"""
Analyze species composition for Eastern vs High Cascades bioregion classification

This script analyzes FIA plot data to identify optimal species composition
thresholds for distinguishing Eastern Cascades (dry forests) from High 
Cascades (mesic/subalpine forests).
"""

import geopandas as gpd
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import json

# File paths
PLOT_DATA = "outputs/plot_carbon_percentiles_latest_surveys.geojson"
EASTERN_BOUNDARY = "outputs/bioregions/eastern_cascades_broad.shp"
HIGH_BOUNDARY = "outputs/bioregions/high_cascades_broad.shp"

def load_and_analyze_plots():
    """Load plot data and analyze species composition patterns"""
    
    print("=== Analyzing Cascades Species Composition ===")
    
    # Load plot data
    print(f"Loading plot data from {PLOT_DATA}")
    plots = gpd.read_file(PLOT_DATA)
    print(f"Loaded {len(plots)} total plots")
    
    # Load boundaries
    eastern_boundary = gpd.read_file(EASTERN_BOUNDARY).to_crs(plots.crs)
    high_boundary = gpd.read_file(HIGH_BOUNDARY).to_crs(plots.crs)
    
    # Find plots within each boundary
    eastern_plots = plots[plots.within(eastern_boundary.unary_union)]
    high_plots = plots[plots.within(high_boundary.unary_union)]
    
    print(f"\nPlots within boundaries:")
    print(f"  Eastern Cascades boundary: {len(eastern_plots)} plots")
    print(f"  High Cascades boundary: {len(high_plots)} plots")
    
    # Define species groups
    eastern_species = [
        'PONDEROSA_PINE_AG_CARBON',
        'LODGEPOLE_PINE_AG_CARBON', 
        'WESTERN_LARCH_AG_CARBON',
        'WESTERN_JUNIPER_AG_CARBON'
    ]
    
    high_cascades_species = [
        'PACIFIC_SILVER_FIR_AG_CARBON',
        'SUBALPINE_FIR_AG_CARBON',
        'MOUNTAIN_HEMLOCK_AG_CARBON',
        'NOBLE_FIR_AG_CARBON',
        'SHASTA_RED_FIR_AG_CARBON',
        'ENGELMANN_SPRUCE_AG_CARBON'
    ]
    
    intermediate_species = [
        'GRAND_FIR_AG_CARBON',
        'DOUGLAS_FIR_AG_CARBON',
        'WESTERN_HEMLOCK_AG_CARBON',
        'WHITE_FIR_AG_CARBON'
    ]
    
    # Calculate species group carbon for all plots
    plots['eastern_carbon'] = plots[eastern_species].sum(axis=1)
    plots['high_cascades_carbon'] = plots[high_cascades_species].sum(axis=1)
    plots['intermediate_carbon'] = plots[intermediate_species].sum(axis=1)
    
    # Calculate ratios
    plots['eastern_ratio'] = plots['eastern_carbon'] / (plots['TOTAL_CARBON_AG'] + 0.1)
    plots['high_cascades_ratio'] = plots['high_cascades_carbon'] / (plots['TOTAL_CARBON_AG'] + 0.1)
    
    # Analyze each region
    print("\n=== Eastern Cascades Plots Analysis ===")
    analyze_region(eastern_plots, "Eastern Cascades", eastern_species, high_cascades_species)
    
    print("\n=== High Cascades Plots Analysis ===")
    analyze_region(high_plots, "High Cascades", eastern_species, high_cascades_species)
    
    # Find optimal thresholds
    print("\n=== Optimal Classification Thresholds ===")
    
    # Eastern Cascades criteria
    eastern_threshold = plots['eastern_carbon'].quantile(0.7)
    print(f"\nEastern Cascades suggested criteria:")
    print(f"  - Eastern species carbon > {eastern_threshold:.0f} kg/ha")
    print(f"  - OR Ponderosa pine carbon > {plots['PONDEROSA_PINE_AG_CARBON'].quantile(0.6):.0f} kg/ha")
    print(f"  - OR Western larch present (> 0 kg/ha)")
    print(f"  - AND Western hemlock < {plots['WESTERN_HEMLOCK_AG_CARBON'].quantile(0.3):.0f} kg/ha")
    
    # High Cascades criteria  
    high_threshold = plots['high_cascades_carbon'].quantile(0.7)
    print(f"\nHigh Cascades suggested criteria:")
    print(f"  - High Cascades species carbon > {high_threshold:.0f} kg/ha")
    print(f"  - OR Pacific silver fir > {plots['PACIFIC_SILVER_FIR_AG_CARBON'].quantile(0.6):.0f} kg/ha")
    print(f"  - OR Mountain hemlock > {plots['MOUNTAIN_HEMLOCK_AG_CARBON'].quantile(0.6):.0f} kg/ha")
    print(f"  - AND Ponderosa pine < {plots['PONDEROSA_PINE_AG_CARBON'].quantile(0.2):.0f} kg/ha")
    
    # Test classification
    test_classification(plots, eastern_plots, high_plots)
    
    return plots

def analyze_region(region_plots, region_name, eastern_species, high_cascades_species):
    """Analyze species composition for a specific region"""
    
    if len(region_plots) == 0:
        print(f"No plots in {region_name}")
        return
    
    print(f"\nSpecies presence in {region_name} (% of plots with species):")
    
    # Eastern species
    print("\nEastern/Dry species:")
    for species in eastern_species:
        present = (region_plots[species] > 0).sum()
        pct = present / len(region_plots) * 100
        avg_carbon = region_plots[region_plots[species] > 0][species].mean() if present > 0 else 0
        print(f"  {species.replace('_AG_CARBON', '')}: {pct:.1f}% of plots, avg {avg_carbon:.0f} kg/ha when present")
    
    # High Cascades species
    print("\nHigh Cascades/Subalpine species:")
    for species in high_cascades_species:
        present = (region_plots[species] > 0).sum()
        pct = present / len(region_plots) * 100
        avg_carbon = region_plots[region_plots[species] > 0][species].mean() if present > 0 else 0
        print(f"  {species.replace('_AG_CARBON', '')}: {pct:.1f}% of plots, avg {avg_carbon:.0f} kg/ha when present")
    
    # Summary statistics
    print(f"\nSummary for {region_name}:")
    print(f"  Mean elevation: {region_plots['ELEVATION_FT'].mean():.0f} ft")
    print(f"  Mean total carbon: {region_plots['TOTAL_CARBON_AG'].mean():.0f} kg/ha")
    print(f"  Mean eastern species carbon: {region_plots['eastern_carbon'].mean():.0f} kg/ha")
    print(f"  Mean high cascades species carbon: {region_plots['high_cascades_carbon'].mean():.0f} kg/ha")

def test_classification(all_plots, eastern_truth, high_truth):
    """Test species-based classification against boundary-defined plots"""
    
    print("\n=== Testing Species-Based Classification ===")
    
    # Define classification criteria
    eastern_criteria = (
        # Primary indicators
        ((all_plots['PONDEROSA_PINE_AG_CARBON'] > 20) | 
         (all_plots['LODGEPOLE_PINE_AG_CARBON'] > 30) |
         (all_plots['WESTERN_LARCH_AG_CARBON'] > 0)) &
        # Exclusion criteria
        (all_plots['WESTERN_HEMLOCK_AG_CARBON'] < 50) &
        (all_plots['PACIFIC_SILVER_FIR_AG_CARBON'] < 30)
    )
    
    high_criteria = (
        # Primary indicators
        ((all_plots['PACIFIC_SILVER_FIR_AG_CARBON'] > 20) |
         (all_plots['SUBALPINE_FIR_AG_CARBON'] > 20) |
         (all_plots['MOUNTAIN_HEMLOCK_AG_CARBON'] > 20) |
         (all_plots['NOBLE_FIR_AG_CARBON'] > 10)) &
        # Exclusion criteria
        (all_plots['PONDEROSA_PINE_AG_CARBON'] < 20) &
        (all_plots['WESTERN_JUNIPER_AG_CARBON'] == 0)
    )
    
    # Apply classification
    eastern_classified = all_plots[eastern_criteria]
    high_classified = all_plots[high_criteria]
    
    print(f"\nClassification results:")
    print(f"  Eastern Cascades: {len(eastern_classified)} plots classified")
    print(f"  High Cascades: {len(high_classified)} plots classified")
    
    # Check overlap
    overlap = eastern_classified.index.intersection(high_classified.index)
    print(f"  Overlap (classified as both): {len(overlap)} plots")
    
    # Save classified plots for visualization
    eastern_classified.to_file('outputs/eastern_cascades_species_classified_plots.geojson')
    high_classified.to_file('outputs/high_cascades_species_classified_plots.geojson')
    
    print("\nSaved classified plots for visualization in QGIS")

if __name__ == "__main__":
    plots = load_and_analyze_plots()
    
    # Create summary report
    summary = {
        'analysis_date': pd.Timestamp.now().isoformat(),
        'eastern_cascades_indicators': {
            'primary_species': ['ponderosa_pine', 'lodgepole_pine', 'western_larch'],
            'exclusion_species': ['western_hemlock', 'pacific_silver_fir'],
            'typical_carbon_range': '50-200 kg/ha eastern species'
        },
        'high_cascades_indicators': {
            'primary_species': ['pacific_silver_fir', 'subalpine_fir', 'mountain_hemlock', 'noble_fir'],
            'exclusion_species': ['ponderosa_pine', 'western_juniper'],
            'typical_carbon_range': '50-300 kg/ha subalpine species'
        }
    }
    
    with open('outputs/cascades_species_analysis_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    
    print("\n✅ Analysis complete! Check outputs for classified plots.")