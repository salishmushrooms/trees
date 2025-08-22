#!/usr/bin/env python3
"""
Create High Cascades Forest Bioregion using Species Composition

This script creates a High Cascades forest bioregion by:
1. Loading FIA plot data with species carbon information
2. Identifying plots dominated by subalpine/true fir species
3. Creating buffers around these plots to form continuous bioregion
4. Applying tree cover constraints within the bioregion
5. Generating output suitable for Mapbox integration

Species-based approach for accurate ecological classification.
"""

import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.geometry import Point, Polygon, MultiPolygon
from shapely.ops import unary_union
import rasterio
from rasterio.mask import mask
from rasterio.features import shapes
from rasterio.warp import calculate_default_transform, reproject, Resampling
import warnings
from pathlib import Path
import json

warnings.filterwarnings('ignore')

# Configuration - High Cascades Forest Parameters
MIN_TREE_COVER_PCT = 20       # Forest definition threshold
BUFFER_DISTANCE_KM = 5         # Buffer around indicator plots (km)
MIN_PLOTS_FOR_CLUSTER = 3      # Minimum plots to form a cluster
OUTPUT_RESOLUTION = 120        # Match elevation data resolution

# File paths
PLOT_DATA = "outputs/plot_carbon_percentiles_latest_surveys.geojson"
BASE_BOUNDARY_SHP = "/Users/JJC/trees/outputs/bioregions/high_cascades_broad.shp"
TREE_COVER_RASTER = "outputs/mapbox_masks/pnw_tree_cover_30m_full.tif"
OUTPUT_DIR = Path("outputs/bioregions")

def load_and_classify_plots():
    """Load plot data and classify High Cascades plots based on species composition"""
    
    print("=== Loading and Classifying High Cascades Plots ===")
    
    # Load plot data
    print(f"Loading plot data from {PLOT_DATA}")
    plots = gpd.read_file(PLOT_DATA)
    
    if plots.crs != 'EPSG:4326':
        plots = plots.to_crs('EPSG:4326')
    
    print(f"Loaded {len(plots)} total plots")
    
    # Load boundary for spatial filtering
    boundary = gpd.read_file(BASE_BOUNDARY_SHP).to_crs('EPSG:4326')
    boundary_geom = boundary.unary_union
    
    # Filter plots within general boundary area
    plots_in_region = plots[plots.within(boundary_geom.buffer(0.1))]  # Small buffer for edge plots
    print(f"Plots within/near boundary: {len(plots_in_region)}")
    
    # Define High Cascades species criteria
    high_cascades_criteria = (
        # Primary subalpine/true fir indicators - at least one must be significant
        (
            (plots_in_region['PACIFIC_SILVER_FIR_AG_CARBON'] > 30) |    # Significant Pacific silver fir
            (plots_in_region['SUBALPINE_FIR_AG_CARBON'] > 30) |         # Significant subalpine fir
            (plots_in_region['MOUNTAIN_HEMLOCK_AG_CARBON'] > 30) |      # Significant mountain hemlock
            (plots_in_region['NOBLE_FIR_AG_CARBON'] > 20) |             # Significant noble fir
            (plots_in_region['SHASTA_RED_FIR_AG_CARBON'] > 20) |        # Significant Shasta red fir
            (plots_in_region['ENGELMANN_SPRUCE_AG_CARBON'] > 20)        # Significant Engelmann spruce
        ) &
        # Exclusion criteria - dry species should be limited
        (
            (plots_in_region['PONDEROSA_PINE_AG_CARBON'] < 30) &        # Limited ponderosa pine
            (plots_in_region['WESTERN_JUNIPER_AG_CARBON'] == 0) &       # No western juniper
            (plots_in_region['WESTERN_LARCH_AG_CARBON'] < 50)           # Limited western larch
        ) &
        # Ensure reasonable total carbon (exclude clearcuts/young stands)
        (plots_in_region['TOTAL_CARBON_AG'] > 40)
    )
    
    # Alternative criteria for mixed high-elevation stands
    mixed_high_criteria = (
        # Subalpine species ratio approach
        (
            (plots_in_region['PACIFIC_SILVER_FIR_AG_CARBON'] + 
             plots_in_region['SUBALPINE_FIR_AG_CARBON'] + 
             plots_in_region['MOUNTAIN_HEMLOCK_AG_CARBON'] +
             plots_in_region['NOBLE_FIR_AG_CARBON'] +
             plots_in_region['ENGELMANN_SPRUCE_AG_CARBON']) / 
            (plots_in_region['TOTAL_CARBON_AG'] + 0.1) > 0.4  # 40% of carbon from subalpine species
        ) &
        (plots_in_region['TOTAL_CARBON_AG'] > 40)
    )
    
    # Additional criteria for high-elevation Douglas-fir/grand fir stands
    high_elevation_mixed = (
        # High elevation indicated by elevation if available
        ((plots_in_region['ELEVATION_FT'] > 4000) if 'ELEVATION_FT' in plots_in_region.columns else True) &
        # Mixed with true firs
        ((plots_in_region['GRAND_FIR_AG_CARBON'] > 50) | 
         (plots_in_region['WHITE_FIR_AG_CARBON'] > 50)) &
        # Some subalpine species present
        ((plots_in_region['PACIFIC_SILVER_FIR_AG_CARBON'] > 0) |
         (plots_in_region['SUBALPINE_FIR_AG_CARBON'] > 0)) &
        # Limited dry species
        (plots_in_region['PONDEROSA_PINE_AG_CARBON'] < 20)
    )
    
    # Combine criteria
    high_cascades_plots = plots_in_region[high_cascades_criteria | mixed_high_criteria | high_elevation_mixed]
    
    print(f"\nClassified {len(high_cascades_plots)} High Cascades indicator plots")
    
    # Summarize species composition
    if len(high_cascades_plots) > 0:
        print("\nSpecies composition of classified plots:")
        
        for species in ['PACIFIC_SILVER_FIR', 'SUBALPINE_FIR', 'MOUNTAIN_HEMLOCK', 
                       'NOBLE_FIR', 'ENGELMANN_SPRUCE', 'GRAND_FIR', 'DOUGLAS_FIR']:
            col_name = f'{species}_AG_CARBON'
            if col_name in high_cascades_plots.columns:
                present = (high_cascades_plots[col_name] > 0).sum()
                if present > 0:
                    avg_carbon = high_cascades_plots[high_cascades_plots[col_name] > 0][col_name].mean()
                    print(f"  {species}: {present}/{len(high_cascades_plots)} plots ({present/len(high_cascades_plots)*100:.1f}%), "
                          f"avg {avg_carbon:.0f} kg/ha when present")
    
    return high_cascades_plots, boundary

def create_bioregion_from_plots(indicator_plots, boundary_gdf):
    """Create bioregion polygons from indicator plots using buffers"""
    
    print("\n=== Creating Bioregion from Indicator Plots ===")
    
    if len(indicator_plots) == 0:
        print("No indicator plots found - cannot create bioregion")
        return None
    
    # Project to Albers for accurate buffering
    plots_albers = indicator_plots.to_crs('EPSG:5070')
    boundary_albers = boundary_gdf.to_crs('EPSG:5070')
    
    # Create buffers around each plot
    buffer_meters = BUFFER_DISTANCE_KM * 1000
    print(f"Creating {BUFFER_DISTANCE_KM}km buffers around {len(plots_albers)} indicator plots")
    
    plot_buffers = plots_albers.geometry.buffer(buffer_meters)
    
    # Merge overlapping buffers
    merged_buffers = unary_union(plot_buffers)
    
    # Clip to boundary
    clipped_region = merged_buffers.intersection(boundary_albers.unary_union)
    
    # Convert back to GeoDataFrame
    if clipped_region.is_empty:
        print("Warning: No valid region after clipping to boundary")
        return None
    
    # Handle both Polygon and MultiPolygon
    if clipped_region.geom_type == 'Polygon':
        polygons = [clipped_region]
    else:
        polygons = list(clipped_region.geoms)
    
    bioregion_gdf = gpd.GeoDataFrame(
        [{'region': 'high_cascades', 'method': 'species_based'} for _ in polygons],
        geometry=polygons,
        crs='EPSG:5070'
    )
    
    # Filter by minimum area
    bioregion_gdf['area_sqkm'] = bioregion_gdf.geometry.area / 1000000
    large_polygons = bioregion_gdf[bioregion_gdf['area_sqkm'] > 1.0]  # Min 1 km²
    
    print(f"Created {len(large_polygons)} bioregion polygons")
    
    # Convert back to WGS84
    bioregion_wgs84 = large_polygons.to_crs('EPSG:4326')
    
    return bioregion_wgs84

def apply_tree_cover_constraint(bioregion_gdf):
    """Apply tree cover constraint to refine bioregion boundaries"""
    
    print("\n=== Applying Tree Cover Constraint ===")
    
    if bioregion_gdf is None or len(bioregion_gdf) == 0:
        return None
    
    # Get combined geometry
    bioregion_geom = bioregion_gdf.unary_union
    
    # Load tree cover raster
    with rasterio.open(TREE_COVER_RASTER) as tcc_src:
        # Mask to bioregion
        try:
            tcc_masked, tcc_transform = mask(tcc_src, [bioregion_geom], crop=True, nodata=tcc_src.nodata)
            tcc_data = tcc_masked[0]
            
            # Apply tree cover threshold
            valid_tcc = ~np.isnan(tcc_data)
            forest_mask = (tcc_data >= MIN_TREE_COVER_PCT) & valid_tcc
            
            print(f"Tree cover pixels ≥ {MIN_TREE_COVER_PCT}%: {forest_mask.sum()}")
            print(f"Approximate forest area: {forest_mask.sum() * 30 * 30 / 1000000:.0f} km²")
            
            # Convert mask to polygons
            forest_polygons = []
            for geom, value in shapes(forest_mask.astype(np.uint8), 
                                    mask=forest_mask, 
                                    transform=tcc_transform):
                if value == 1:
                    from shapely.geometry import shape
                    poly = shape(geom)
                    if poly.is_valid and poly.area > 0:
                        forest_polygons.append(poly)
            
            print(f"Extracted {len(forest_polygons)} forest polygons")
            
            if not forest_polygons:
                print("Warning: No forest polygons extracted")
                return bioregion_gdf
            
            # Create GeoDataFrame
            forest_gdf = gpd.GeoDataFrame(
                geometry=forest_polygons,
                crs=tcc_src.crs
            ).to_crs('EPSG:4326')
            
            # Clean up small fragments
            forest_gdf['area_sqkm'] = forest_gdf.geometry.area * 111 * 111
            large_forests = forest_gdf[forest_gdf['area_sqkm'] > 0.1]  # Min 0.1 km²
            
            # Merge polygons
            if len(large_forests) > 0:
                merged_forest = unary_union(large_forests.geometry)
                
                final_gdf = gpd.GeoDataFrame(
                    [{
                        'region_name': 'High Cascades Forest',
                        'region_code': 'high_cascades_species',
                        'description': f'High Cascades forest based on species composition with ≥{MIN_TREE_COVER_PCT}% tree cover',
                        'min_tree_cover_pct': MIN_TREE_COVER_PCT,
                        'classification_method': 'species_composition',
                        'indicator_plots': len(indicator_plots),
                        'created_date': pd.Timestamp.now().isoformat()
                    }],
                    geometry=[merged_forest],
                    crs='EPSG:4326'
                )
                
                return final_gdf
            
        except Exception as e:
            print(f"Error applying tree cover constraint: {e}")
            return bioregion_gdf
    
    return None

def main():
    """Main processing function"""
    print("=== Creating High Cascades Forest Bioregion (Species-Based) ===")
    
    # Load and classify plots
    indicator_plots, boundary = load_and_classify_plots()
    
    if len(indicator_plots) == 0:
        print("No High Cascades indicator plots found")
        return
    
    # Save indicator plots for review
    indicator_file = OUTPUT_DIR / "high_cascades_indicator_plots.geojson"
    indicator_plots.to_file(indicator_file, driver='GeoJSON')
    print(f"\nSaved indicator plots to {indicator_file}")
    
    # Create bioregion from plots
    bioregion_raw = create_bioregion_from_plots(indicator_plots, boundary)
    
    if bioregion_raw is None:
        print("Failed to create bioregion from plots")
        return
    
    # Apply tree cover constraint
    bioregion_final = apply_tree_cover_constraint(bioregion_raw)
    
    if bioregion_final is None:
        print("Failed to apply tree cover constraint")
        return
    
    # Calculate final area
    bioregion_final['area_sqkm'] = bioregion_final.geometry.area * 111 * 111
    total_area = bioregion_final['area_sqkm'].sum()
    
    print(f"\n=== Final Bioregion Statistics ===")
    print(f"Total area: {total_area:.0f} km²")
    print(f"Based on {len(indicator_plots)} indicator plots")
    
    # Save final bioregion
    output_file = OUTPUT_DIR / "high_cascades_species_based.geojson"
    bioregion_final.to_file(output_file, driver='GeoJSON')
    
    # Create summary
    summary = {
        'bioregion_name': 'High Cascades Forest (Species-Based)',
        'creation_date': pd.Timestamp.now().isoformat(),
        'method': 'species_composition_classification',
        'parameters': {
            'min_tree_cover_pct': MIN_TREE_COVER_PCT,
            'buffer_distance_km': BUFFER_DISTANCE_KM,
            'min_plots_for_cluster': MIN_PLOTS_FOR_CLUSTER
        },
        'species_criteria': {
            'primary_indicators': [
                'pacific_silver_fir > 30 kg/ha',
                'subalpine_fir > 30 kg/ha',
                'mountain_hemlock > 30 kg/ha',
                'noble_fir > 20 kg/ha',
                'engelmann_spruce > 20 kg/ha'
            ],
            'exclusion_criteria': [
                'ponderosa_pine < 30 kg/ha',
                'western_juniper = 0 kg/ha',
                'western_larch < 50 kg/ha'
            ],
            'mixed_criteria': '40% carbon from subalpine species'
        },
        'results': {
            'indicator_plots': len(indicator_plots),
            'total_area_sqkm': float(total_area),
            'polygon_count': len(bioregion_final)
        }
    }
    
    summary_file = OUTPUT_DIR / "high_cascades_species_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"\n✅ High Cascades species-based bioregion created successfully!")
    print(f"   Output: {output_file}")
    print(f"   Indicator plots: {indicator_file}")
    print(f"   Summary: {summary_file}")

if __name__ == "__main__":
    main()