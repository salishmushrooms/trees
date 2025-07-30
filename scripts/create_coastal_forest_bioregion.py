#!/usr/bin/env python3
"""
Create Coastal Forest Bioregion Vector

This script creates a coastal forest bioregion polygon that:
1. Extends from Port Angeles on the north Olympic Peninsula
2. Wraps around the Pacific Coast 
3. Extends south to approximately 41.94882°, -124.20744°
4. Uses Sitka spruce and shore pine (Pinus contorta) as indicator species
5. Applies elevation constraints (typically <1000 ft)
6. Uses tree cover data for accurate land/sea boundaries
"""

import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.geometry import Point, Polygon, MultiPolygon
from shapely.ops import unary_union
import rasterio
from rasterio.mask import mask
import warnings
from pathlib import Path
from tqdm import tqdm
import json

warnings.filterwarnings('ignore')

# Configuration
NORTH_BOUNDARY = 48.5  # North of Port Angeles
SOUTH_BOUNDARY = 41.94882  # Southern extent
WEST_BOUNDARY = -125.0  # Pacific Ocean boundary
EAST_BOUNDARY = -122.5  # Inland boundary (will be refined by elevation/species)
MAX_ELEVATION_FT = 1000  # Maximum elevation for coastal forests
COASTAL_BUFFER_KM = 50  # Maximum distance from coast in km

# Species codes
SITKA_SPRUCE = 98
LODGEPOLE_PINE = 108  # Includes shore pine subspecies

# File paths
PLOT_DATA = "outputs/plot_carbon_percentiles_latest_surveys.geojson"
COASTAL_PLOTS = "outputs/coastal_analysis/coastal_forest_plots.geojson"  # Pre-extracted data
TCC_MASK_5PCT = "cache/species_masks/pnw_tcc_mask_5pct_30m_0.5sqkm.geojson"
ELEVATION_MASK = "cache/species_masks/pnw_elevation_mask_0_1000ft.geojson"  # May need to create
OUTPUT_DIR = Path("outputs/bioregions")

def load_coastal_species_plots():
    """Load plots containing Sitka spruce or shore pine at low elevations"""
    print("Loading coastal species plots...")
    
    # Check if pre-extracted coastal data exists
    coastal_path = Path(COASTAL_PLOTS)
    if coastal_path.exists():
        print(f"Loading pre-extracted data from {coastal_path}")
        coastal_gdf = gpd.read_file(coastal_path)
        
        # Filter for species of interest
        if 'SITKA_SPRUCE_AG_CARBON' in coastal_gdf.columns:
            sitka_mask = coastal_gdf['SITKA_SPRUCE_AG_CARBON'] > 0
        else:
            sitka_mask = False
            
        if 'LODGEPOLE_PINE_AG_CARBON' in coastal_gdf.columns:
            pine_mask = coastal_gdf['LODGEPOLE_PINE_AG_CARBON'] > 0
        else:
            pine_mask = False
        
        # Get plots with either species
        species_plots = coastal_gdf[sitka_mask | pine_mask]
        print(f"Found {len(species_plots)} plots with Sitka spruce or shore pine")
        
        return species_plots
    else:
        # Fallback to loading from main plot data
        gdf = gpd.read_file(PLOT_DATA)
        
        # Filter for coastal region bounds
        coastal_gdf = gdf[
            (gdf.geometry.y >= SOUTH_BOUNDARY) &
            (gdf.geometry.y <= NORTH_BOUNDARY) &
            (gdf.geometry.x <= EAST_BOUNDARY) &
            (gdf.geometry.x >= WEST_BOUNDARY)
        ]
        
        # Filter for low elevation (handle different column names)
        elev_col = 'ELEVATION_FT' if 'ELEVATION_FT' in coastal_gdf.columns else 'ELEV'
        coastal_gdf = coastal_gdf[coastal_gdf[elev_col] <= MAX_ELEVATION_FT]
        
        print(f"Found {len(coastal_gdf)} plots in coastal elevation range")
        
        return coastal_gdf

def create_coastal_boundary_polygon():
    """Create initial coastal boundary polygon"""
    # Define key coastal points from north to south
    coastal_points = [
        # North Olympic Peninsula (Port Angeles area)
        (-123.4, 48.12),  # Port Angeles
        (-124.0, 48.15),  # Strait of Juan de Fuca
        (-124.6, 48.38),  # Cape Flattery area
        
        # Pacific Coast - Olympic Peninsula
        (-124.7, 47.9),   # Northern WA coast
        (-124.4, 47.5),   # Central WA coast
        (-124.1, 47.0),   # Southern WA coast
        
        # Oregon Coast
        (-124.0, 46.2),   # Columbia River mouth
        (-123.9, 45.5),   # Northern OR coast
        (-124.1, 44.5),   # Central OR coast
        (-124.3, 43.5),   # Southern OR coast
        (-124.4, 42.5),   # Far southern OR coast
        
        # California North Coast
        (-124.20744, 41.94882),  # Southern boundary
        
        # Inland boundaries (will be refined)
        (-123.5, 41.94882),  # Inland southern point
        (-123.0, 43.0),      # Inland OR
        (-122.8, 45.0),      # Inland OR/WA
        (-122.5, 47.0),      # Inland WA
        (-122.8, 48.0),      # Inland north WA
        (-123.0, 48.12),     # Back to Port Angeles inland
    ]
    
    return Polygon(coastal_points)

def apply_tree_cover_mask(geometry, tcc_threshold=5):
    """Apply tree cover mask to ensure we only include forested areas"""
    print(f"Applying tree cover mask (>{tcc_threshold}%)...")
    
    # Check if regional TCC mask exists
    tcc_mask_path = Path(TCC_MASK_5PCT)
    if tcc_mask_path.exists():
        try:
            print(f"Loading cached TCC mask from {tcc_mask_path}")
            print(f"File size: {tcc_mask_path.stat().st_size / 1024 / 1024:.1f} MB")
            
            # For large files, read in chunks or skip
            if tcc_mask_path.stat().st_size > 100 * 1024 * 1024:  # 100MB
                print("Warning: TCC mask file is very large, skipping for now")
                print("Consider using a simplified mask or higher minimum area threshold")
                return geometry
                
            tcc_mask_gdf = gpd.read_file(tcc_mask_path)
            
            # Intersect with our coastal boundary
            if not tcc_mask_gdf.empty:
                print(f"Loaded {len(tcc_mask_gdf)} tree cover polygons")
                tcc_geometry = unary_union(tcc_mask_gdf.geometry)
                masked_geometry = geometry.intersection(tcc_geometry)
                return masked_geometry
        except Exception as e:
            print(f"Error loading TCC mask: {e}")
            print("Proceeding without tree cover constraint")
    else:
        print(f"Warning: TCC mask not found at {tcc_mask_path}")
        print("Proceeding without tree cover constraint")
    
    return geometry

def apply_elevation_constraint(geometry):
    """Apply elevation mask to limit to low elevations"""
    print(f"Applying elevation constraint (<{MAX_ELEVATION_FT} ft)...")
    
    # Check for existing elevation mask
    elev_mask_path = Path(ELEVATION_MASK)
    if elev_mask_path.exists():
        print(f"Loading cached elevation mask from {elev_mask_path}")
        elev_mask_gdf = gpd.read_file(elev_mask_path)
        
        if not elev_mask_gdf.empty:
            elev_geometry = unary_union(elev_mask_gdf.geometry)
            masked_geometry = geometry.intersection(elev_geometry)
            return masked_geometry
    else:
        print(f"Warning: Elevation mask not found at {elev_mask_path}")
        print("Consider creating it with: create_regional_elevation_mask.py")
    
    return geometry

def refine_with_species_data(boundary_polygon, species_plots):
    """Refine boundary using actual species occurrence data"""
    print("Refining boundary with species occurrence data...")
    
    if len(species_plots) == 0:
        print("Warning: No species plots found for refinement")
        return boundary_polygon
    
    # Create buffers around species occurrences
    buffer_distance = 0.1  # degrees (~10 km)
    species_buffers = []
    
    for idx, plot in species_plots.iterrows():
        buffer = plot.geometry.buffer(buffer_distance)
        species_buffers.append(buffer)
    
    # Union all buffers
    species_region = unary_union(species_buffers)
    
    # Intersect with original boundary
    refined_boundary = boundary_polygon.intersection(species_region)
    
    return refined_boundary

def simplify_for_web(geometry, tolerance=0.001):
    """Simplify geometry for web mapping"""
    if hasattr(geometry, 'simplify'):
        return geometry.simplify(tolerance)
    return geometry

def main():
    """Main processing function"""
    print("=== Creating Coastal Forest Bioregion ===")
    
    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Load species data
    coastal_plots = load_coastal_species_plots()
    
    # Step 2: Create initial coastal boundary
    print("\nCreating initial coastal boundary polygon...")
    coastal_boundary = create_coastal_boundary_polygon()
    
    # Step 3: Apply tree cover mask
    coastal_forest = apply_tree_cover_mask(coastal_boundary)
    
    # Step 4: Apply elevation constraint
    coastal_forest = apply_elevation_constraint(coastal_forest)
    
    # Step 5: Refine with species data (if available)
    # Note: This step requires species-specific columns in the plot data
    # coastal_forest = refine_with_species_data(coastal_forest, coastal_plots)
    
    # Step 6: Simplify for web use
    print("\nSimplifying geometry for web mapping...")
    coastal_forest_simple = simplify_for_web(coastal_forest, tolerance=0.001)
    
    # Create GeoDataFrame
    coastal_bioregion = gpd.GeoDataFrame(
        [{
            'region_name': 'Coastal Forests',
            'region_code': 'coastal',
            'description': 'Low-elevation Pacific coastal forests characterized by Sitka spruce and shore pine',
            'elevation_max_ft': MAX_ELEVATION_FT,
            'north_bound': NORTH_BOUNDARY,
            'south_bound': SOUTH_BOUNDARY,
            'indicator_species': 'Sitka spruce (Picea sitchensis), Shore pine (Pinus contorta var. contorta)',
            'created_date': pd.Timestamp.now().isoformat(),
            'notes': 'FIA data does not differentiate P. contorta subspecies'
        }],
        geometry=[coastal_forest_simple],
        crs='EPSG:4326'
    )
    
    # Save outputs
    output_file = OUTPUT_DIR / "coastal_forest_bioregion.geojson"
    print(f"\nSaving to {output_file}")
    coastal_bioregion.to_file(output_file, driver='GeoJSON')
    
    # Also save detailed version
    detailed_output = OUTPUT_DIR / "coastal_forest_bioregion_detailed.geojson"
    coastal_bioregion_detailed = coastal_bioregion.copy()
    coastal_bioregion_detailed.geometry = [coastal_forest]  # Non-simplified
    coastal_bioregion_detailed.to_file(detailed_output, driver='GeoJSON')
    
    # Print summary statistics
    print("\n=== Summary ===")
    print(f"Total area: {coastal_bioregion.geometry.area.sum():.4f} square degrees")
    print(f"Approximate area: {coastal_bioregion.geometry.area.sum() * 111 * 111:.0f} km²")
    print(f"Number of coastal plots found: {len(coastal_plots)}")
    
    # Save plot points for reference
    if len(coastal_plots) > 0:
        plot_output = OUTPUT_DIR / "coastal_forest_plots.geojson"
        coastal_plots.to_file(plot_output, driver='GeoJSON')
        print(f"Saved plot data to {plot_output}")
    
    print("\n✅ Coastal forest bioregion created successfully!")
    print(f"   Main output: {output_file}")
    print(f"   Detailed output: {detailed_output}")

if __name__ == "__main__":
    main()