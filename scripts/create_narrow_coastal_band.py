#!/usr/bin/env python3
"""
Create Narrow Coastal Forest Band

Creates a narrow coastal forest bioregion that closely follows the coastline,
similar to the reference bioregional map. Uses a gradient approach where the
band is widest where there are the most coastal species occurrences.
"""

import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.geometry import Point, Polygon, MultiPolygon, LineString
from shapely.ops import unary_union
from pathlib import Path
import warnings

warnings.filterwarnings('ignore')

# Configuration
MAX_COASTAL_WIDTH_KM = 30   # Maximum width of coastal band
MIN_COASTAL_WIDTH_KM = 5    # Minimum width of coastal band
MAX_ELEVATION_FT = 1000     # Maximum elevation for coastal forests

# Geographic bounds
NORTH_BOUND = 48.4
SOUTH_BOUND = 41.95

# High-resolution coastal reference points (longitude, latitude)
DETAILED_COASTAL_POINTS = [
    # Strait of Juan de Fuca - Port Angeles to Pacific
    (-123.0, 48.12),  # Port Angeles
    (-123.2, 48.13),  # East of PA
    (-123.4, 48.14),  # Mid-strait
    (-123.6, 48.16),  # West mid-strait
    (-123.8, 48.20),  # Approaching Neah Bay
    (-124.0, 48.25),  # Neah Bay area
    (-124.3, 48.35),  # Cape Flattery approach
    (-124.6, 48.38),  # Cape Flattery
    
    # Olympic Peninsula Pacific Coast (north to south)
    (-124.65, 48.2),  # South of Cape Flattery
    (-124.7, 47.9),   # Northern WA coast
    (-124.65, 47.7),  # Coastal curve
    (-124.6, 47.5),   # La Push area
    (-124.5, 47.3),   # South of La Push
    (-124.4, 47.1),   # Coastal bend
    (-124.2, 47.0),   # Westport area
    (-124.1, 46.8),   # South of Westport
    
    # Washington/Oregon Coast transition
    (-124.0, 46.3),   # Columbia River mouth area
    (-123.95, 46.1),  # North of Seaside
    (-123.92, 45.9),  # Seaside area
    (-123.9, 45.8),   # South of Seaside
    (-123.9, 45.6),   # Lincoln City approach
    (-123.95, 45.4),  # Lincoln City area
    
    # Central Oregon Coast
    (-124.0, 45.0),   # Central OR coast
    (-124.05, 44.8),  # Newport approach
    (-124.1, 44.6),   # Newport area
    (-124.1, 44.4),   # South of Newport
    (-124.15, 44.2),  # Florence area
    (-124.2, 44.0),   # South of Florence
    (-124.25, 43.8),  # Coos Bay approach
    (-124.3, 43.6),   # Coos Bay area
    
    # Southern Oregon Coast  
    (-124.35, 43.4),  # South of Coos Bay
    (-124.4, 43.2),   # Bandon area
    (-124.42, 43.0),  # South of Bandon
    (-124.4, 42.8),   # Gold Beach area
    (-124.35, 42.6),  # South of Gold Beach
    (-124.3, 42.4),   # Brookings approach
    (-124.25, 42.2),  # Brookings area
    (-124.2, 41.95),  # Southern extent
]

OUTPUT_DIR = Path("outputs/bioregions")

def create_detailed_coastline():
    """Create detailed coastline from reference points"""
    print("Creating detailed coastline reference...")
    
    # Sort points by latitude (north to south)
    points = sorted(DETAILED_COASTAL_POINTS, key=lambda p: p[1], reverse=True)
    coastline = LineString(points)
    
    print(f"Created coastline with {len(points)} reference points")
    return coastline, points

def calculate_species_density_weights():
    """Calculate species density to weight coastal band width"""
    print("Loading coastal species data for density weighting...")
    
    try:
        coastal_plots_path = Path("outputs/coastal_analysis/coastal_forest_plots.geojson")
        if coastal_plots_path.exists():
            plots_gdf = gpd.read_file(coastal_plots_path)
            
            # Create grid cells and count species density
            lat_bins = np.linspace(SOUTH_BOUND, NORTH_BOUND, 15)  # ~0.5 degree bins
            density_weights = {}
            
            for i in range(len(lat_bins) - 1):
                lat_min, lat_max = lat_bins[i], lat_bins[i + 1]
                lat_center = (lat_min + lat_max) / 2
                
                # Count plots in this latitude band
                in_band = plots_gdf[
                    (plots_gdf.geometry.y >= lat_min) & 
                    (plots_gdf.geometry.y < lat_max)
                ]
                
                # Count coastal indicator species
                sitka_count = len(in_band[in_band.get('SITKA_SPRUCE_AG_CARBON', 0) > 0]) if 'SITKA_SPRUCE_AG_CARBON' in in_band.columns else 0
                pine_count = len(in_band[in_band.get('LODGEPOLE_PINE_AG_CARBON', 0) > 0]) if 'LODGEPOLE_PINE_AG_CARBON' in in_band.columns else 0
                
                total_coastal_species = sitka_count + pine_count
                density_weights[lat_center] = total_coastal_species
                
                print(f"  Lat {lat_center:.1f}°: {total_coastal_species} coastal species plots")
            
            return density_weights
            
    except Exception as e:
        print(f"Could not load species data: {e}")
    
    # Default uniform weighting
    return {}

def create_variable_width_coastal_band(coastline_points, density_weights):
    """Create coastal band with variable width based on species density"""
    print("Creating variable-width coastal band...")
    
    band_polygons = []
    
    # Convert km to degrees (approximate)
    km_to_deg = 1.0 / 111.0
    
    for i in range(len(coastline_points) - 1):
        lon1, lat1 = coastline_points[i]
        lon2, lat2 = coastline_points[i + 1]
        
        # Find nearest density weight
        lat_center = (lat1 + lat2) / 2
        nearest_lat = min(density_weights.keys(), key=lambda x: abs(x - lat_center)) if density_weights else lat_center
        species_count = density_weights.get(nearest_lat, 5)  # Default weight
        
        # Calculate width based on species density
        # More species = wider band (up to MAX_COASTAL_WIDTH_KM)
        if species_count > 0:
            width_factor = min(species_count / 20.0, 1.0)  # Normalize to 0-1
            width_km = MIN_COASTAL_WIDTH_KM + (MAX_COASTAL_WIDTH_KM - MIN_COASTAL_WIDTH_KM) * width_factor
        else:
            width_km = MIN_COASTAL_WIDTH_KM
        
        width_deg = width_km * km_to_deg
        
        # Create coastal segment polygon
        # Extend inland (eastward) from coast points
        segment_poly = Polygon([
            (lon1, lat1),                    # Coast point 1
            (lon2, lat2),                    # Coast point 2
            (lon2 + width_deg, lat2),        # Inland point 2  
            (lon1 + width_deg, lat1),        # Inland point 1
            (lon1, lat1)                     # Close polygon
        ])
        
        band_polygons.append(segment_poly)
        
        if i % 10 == 0:  # Progress update
            print(f"  Segment {i+1}/{len(coastline_points)-1}: {width_km:.1f}km width ({species_count} species)")
    
    # Union all segments
    coastal_band = unary_union(band_polygons)
    print(f"Created coastal band from {len(band_polygons)} segments")
    
    return coastal_band

def create_smooth_coastal_band(coastline_points):
    """Create smooth coastal band with consistent narrow width"""
    print("Creating smooth narrow coastal band...")
    
    # Simple approach: create narrow buffer around all coastal points
    km_to_deg = 1.0 / 111.0
    width_deg = 25 * km_to_deg  # 25km width
    
    # Create polygons for north-south segments
    band_polygons = []
    
    # Group points by major coastal segments
    strait_points = [p for p in coastline_points if p[1] >= 48.0]  # Strait of Juan de Fuca
    wa_coast_points = [p for p in coastline_points if 46.0 <= p[1] < 48.0]  # WA Pacific coast
    or_coast_points = [p for p in coastline_points if p[1] < 46.0]  # OR coast
    
    segments = [
        ("Strait of Juan de Fuca", strait_points),
        ("Washington Coast", wa_coast_points), 
        ("Oregon Coast", or_coast_points)
    ]
    
    for segment_name, points in segments:
        if len(points) < 2:
            continue
            
        print(f"  Creating {segment_name} segment ({len(points)} points)")
        
        # Find bounding box for this segment
        lons = [p[0] for p in points]
        lats = [p[1] for p in points]
        
        min_lon, max_lon = min(lons), max(lons)
        min_lat, max_lat = min(lats), max(lats)
        
        # Create coastal band polygon for this segment
        # Extend eastward from the westernmost point
        segment_poly = Polygon([
            (min_lon, min_lat),                    # SW corner
            (min_lon, max_lat),                    # NW corner
            (min_lon + width_deg, max_lat),        # NE corner (inland)
            (min_lon + width_deg, min_lat),        # SE corner (inland)
            (min_lon, min_lat)                     # Close polygon
        ])
        
        band_polygons.append(segment_poly)
    
    # Union all segments
    coastal_band = unary_union(band_polygons)
    return coastal_band

def main():
    """Create narrow coastal forest band"""
    print("=== Creating Narrow Coastal Forest Band ===\n")
    
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Create detailed coastline
    coastline, coastline_points = create_detailed_coastline()
    
    # Step 2: Calculate species density weights
    density_weights = calculate_species_density_weights()
    
    # Step 3: Create narrow coastal band
    if density_weights:
        print("\nUsing variable-width approach based on species density...")
        coastal_band = create_variable_width_coastal_band(coastline_points, density_weights)
    else:
        print("\nUsing smooth narrow band approach...")
        coastal_band = create_smooth_coastal_band(coastline_points)
    
    # Step 4: Validate with species data
    print("\nValidating coastal band...")
    try:
        coastal_plots_path = Path("outputs/coastal_analysis/coastal_forest_plots.geojson")
        if coastal_plots_path.exists():
            plots_gdf = gpd.read_file(coastal_plots_path)
            within_band = plots_gdf[plots_gdf.within(coastal_band)]
            
            sitka_plots = len(within_band[within_band.get('SITKA_SPRUCE_AG_CARBON', 0) > 0]) if 'SITKA_SPRUCE_AG_CARBON' in within_band.columns else 0
            pine_plots = len(within_band[within_band.get('LODGEPOLE_PINE_AG_CARBON', 0) > 0]) if 'LODGEPOLE_PINE_AG_CARBON' in within_band.columns else 0
            
            validation_plots = len(within_band)
            print(f"  Captures {validation_plots}/{len(plots_gdf)} coastal plots ({validation_plots/len(plots_gdf)*100:.1f}%)")
            print(f"  Contains {sitka_plots} Sitka spruce plots, {pine_plots} shore pine plots")
    except Exception as e:
        print(f"  Validation failed: {e}")
        validation_plots = 0
    
    # Step 5: Create final GeoDataFrame
    coastal_bioregion = gpd.GeoDataFrame(
        [{
            'region_name': 'Narrow Coastal Forest Band',
            'region_code': 'coastal_narrow',
            'description': 'Narrow Pacific coastal forest band following coastline closely',
            'max_width_km': MAX_COASTAL_WIDTH_KM,
            'min_width_km': MIN_COASTAL_WIDTH_KM,
            'max_elevation_ft': MAX_ELEVATION_FT,
            'validation_plots': validation_plots,
            'method': 'coastline_following_variable_width',
            'created_date': pd.Timestamp.now().isoformat()
        }],
        geometry=[coastal_band],
        crs='EPSG:4326'
    )
    
    # Step 6: Save output
    output_file = OUTPUT_DIR / "coastal_forest_narrow_band.geojson"
    print(f"\nSaving to {output_file}")
    coastal_bioregion.to_file(output_file, driver='GeoJSON')
    
    # Step 7: Summary
    area_km2 = coastal_bioregion.to_crs('EPSG:5070').geometry.area.sum() / 1e6
    
    print(f"\n=== Summary ===")
    print(f"Narrow coastal band area: {area_km2:.0f} km²")
    print(f"Width range: {MIN_COASTAL_WIDTH_KM}-{MAX_COASTAL_WIDTH_KM} km")
    print(f"Extent: {SOUTH_BOUND}°N to {NORTH_BOUND}°N")
    print(f"Species validation: {validation_plots} plots captured")
    
    print(f"\n✅ Narrow coastal forest band created!")
    print(f"   Output: {output_file}")

if __name__ == "__main__":
    main()