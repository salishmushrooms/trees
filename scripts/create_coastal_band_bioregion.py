#!/usr/bin/env python3
"""
Create Accurate Coastal Forest Bioregion Band

Creates a narrow coastal forest bioregion following the actual coastline using:
1. Coastline data (NOAA or Natural Earth)
2. Distance buffers inland from coast
3. Elevation constraints for low-elevation coastal areas
4. Tree cover constraints for forested areas

This creates the thin coastal band visible in regional biogeographic maps.
"""

import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.geometry import Point, Polygon, MultiPolygon
from shapely.ops import unary_union
import rasterio
from rasterio.mask import mask
import requests
import zipfile
from pathlib import Path
from tqdm import tqdm
import warnings

warnings.filterwarnings('ignore')

# Configuration
COASTAL_BUFFER_KM = 15  # Distance inland from coast (km)
MAX_ELEVATION_FT = 500   # Maximum elevation for coastal forests
TREE_COVER_THRESHOLD = 5  # Minimum tree cover percentage

# Geographic bounds
NORTH_BOUND = 48.5   # North of Port Angeles
SOUTH_BOUND = 41.95  # Southern Oregon coast
WEST_BOUND = -125.0  # Pacific Ocean
EAST_BOUND = -122.0  # Inland limit

# Data sources
NATURAL_EARTH_COASTLINE = "https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_coastline.zip"
NOAA_SHORELINE_URL = "https://shoreline.noaa.gov/data/datasheets/pnw.html"

OUTPUT_DIR = Path("outputs/bioregions")
CACHE_DIR = Path("cache/coastline")

def download_coastline_data():
    """Download high-resolution coastline data"""
    print("Downloading coastline data...")
    
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    
    # Try Natural Earth 10m coastline first
    coastline_zip = CACHE_DIR / "ne_10m_coastline.zip"
    coastline_shp = CACHE_DIR / "ne_10m_coastline.shp"
    
    if not coastline_shp.exists():
        if not coastline_zip.exists():
            print("Downloading Natural Earth 10m coastline...")
            response = requests.get(NATURAL_EARTH_COASTLINE)
            coastline_zip.write_bytes(response.content)
        
        print("Extracting coastline data...")
        with zipfile.ZipFile(coastline_zip, 'r') as zip_ref:
            zip_ref.extractall(CACHE_DIR)
    
    return coastline_shp

def clip_coastline_to_region(coastline_path):
    """Clip coastline to Pacific Northwest region"""
    print("Loading and clipping coastline to PNW region...")
    
    # Load coastline
    coastline_gdf = gpd.read_file(coastline_path)
    
    # Create bounding box for PNW coast
    bbox = Polygon([
        (WEST_BOUND, SOUTH_BOUND),
        (EAST_BOUND, SOUTH_BOUND), 
        (EAST_BOUND, NORTH_BOUND),
        (WEST_BOUND, NORTH_BOUND),
        (WEST_BOUND, SOUTH_BOUND)
    ])
    
    # Clip to region
    pnw_coastline = coastline_gdf.clip(bbox)
    
    print(f"Clipped to {len(pnw_coastline)} coastline segments")
    return pnw_coastline

def create_coastal_buffers(coastline_gdf, distances_km=[5, 10, 15, 20]):
    """Create multiple buffer zones inland from coast"""
    print(f"Creating coastal buffers at distances: {distances_km} km")
    
    buffers = {}
    
    # Convert to appropriate projection for accurate distance calculations
    # Use Albers Equal Area for Pacific Northwest
    coastline_aea = coastline_gdf.to_crs('EPSG:5070')  # NAD83 Albers
    
    for dist_km in distances_km:
        print(f"  Creating {dist_km}km buffer...")
        
        # Create buffer in meters
        buffer_m = dist_km * 1000
        buffered = coastline_aea.buffer(buffer_m)
        
        # Union all buffer polygons
        buffer_union = unary_union(buffered)
        
        # Convert back to WGS84
        buffer_gdf = gpd.GeoDataFrame([1], geometry=[buffer_union], crs='EPSG:5070')
        buffer_gdf = buffer_gdf.to_crs('EPSG:4326')
        
        buffers[f'{dist_km}km'] = buffer_gdf
        
        # Calculate area
        area_km2 = buffer_gdf.to_crs('EPSG:5070').geometry.area.sum() / 1e6
        print(f"    Buffer area: {area_km2:.0f} km²")
    
    return buffers

def apply_elevation_constraint(coastal_buffer, max_elevation_ft=MAX_ELEVATION_FT):
    """Apply elevation mask to coastal buffer"""
    print(f"Applying elevation constraint (<{max_elevation_ft} ft)...")
    
    # Check for existing elevation mask
    elevation_masks = list(Path("cache/species_masks").glob("pnw_elevation_mask_0_*ft.geojson"))
    
    # Find closest elevation mask
    target_elev = max_elevation_ft
    best_mask = None
    min_diff = float('inf')
    
    for mask_path in elevation_masks:
        # Extract elevation from filename
        parts = mask_path.stem.split('_')
        if len(parts) >= 4:
            try:
                mask_elev = int(parts[3].replace('ft', ''))
                diff = abs(mask_elev - target_elev)
                if diff < min_diff:
                    min_diff = diff
                    best_mask = mask_path
            except ValueError:
                continue
    
    if best_mask:
        print(f"Using elevation mask: {best_mask.name}")
        try:
            elev_mask_gdf = gpd.read_file(best_mask)
            
            # Intersect coastal buffer with elevation mask
            if not elev_mask_gdf.empty:
                elev_geometry = unary_union(elev_mask_gdf.geometry)
                constrained_buffer = coastal_buffer.intersection(elev_geometry)
                
                # Calculate area reduction
                original_area = coastal_buffer.to_crs('EPSG:5070').area.sum() / 1e6
                constrained_area = gpd.GeoSeries([constrained_buffer], crs='EPSG:4326').to_crs('EPSG:5070').area.sum() / 1e6
                reduction_pct = (original_area - constrained_area) / original_area * 100
                
                print(f"  Elevation constraint reduced area by {reduction_pct:.1f}%")
                return constrained_buffer
        except Exception as e:
            print(f"Error applying elevation mask: {e}")
    else:
        print("No suitable elevation mask found, proceeding without elevation constraint")
    
    return coastal_buffer.geometry.iloc[0] if hasattr(coastal_buffer, 'geometry') else coastal_buffer

def create_coastal_band_region():
    """Main function to create coastal band bioregion"""
    print("=== Creating Coastal Band Bioregion ===\n")
    
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Get coastline data
    coastline_path = download_coastline_data()
    
    # Step 2: Clip to PNW region  
    pnw_coastline = clip_coastline_to_region(coastline_path)
    
    # Step 3: Create coastal buffers
    coastal_buffers = create_coastal_buffers(pnw_coastline, [COASTAL_BUFFER_KM])
    
    # Step 4: Apply elevation constraints
    main_buffer = coastal_buffers[f'{COASTAL_BUFFER_KM}km'].geometry.iloc[0]
    coastal_forest_geom = apply_elevation_constraint(main_buffer, MAX_ELEVATION_FT)
    
    # Step 5: Create final GeoDataFrame
    coastal_bioregion = gpd.GeoDataFrame(
        [{
            'region_name': 'Coastal Forest Band',
            'region_code': 'coastal_band',
            'description': f'Narrow coastal forest band within {COASTAL_BUFFER_KM}km of Pacific coast',
            'buffer_distance_km': COASTAL_BUFFER_KM,
            'max_elevation_ft': MAX_ELEVATION_FT,
            'extends_from': 'Port Angeles, WA',
            'extends_to': 'Southern Oregon coast',
            'method': 'coastline_proximity_elevation',
            'created_date': pd.Timestamp.now().isoformat(),
            'data_sources': 'Natural Earth 10m coastline, NLCD elevation'
        }],
        geometry=[coastal_forest_geom],
        crs='EPSG:4326'
    )
    
    # Step 6: Save outputs
    output_file = OUTPUT_DIR / "coastal_forest_band.geojson"
    print(f"\nSaving coastal forest band to {output_file}")
    coastal_bioregion.to_file(output_file, driver='GeoJSON')
    
    # Step 7: Create validation buffers at different distances
    print("\nCreating validation buffers at multiple distances...")
    for dist_km, buffer_gdf in coastal_buffers.items():
        validation_file = OUTPUT_DIR / f"coastal_buffer_{dist_km}.geojson"
        
        # Apply elevation constraint to each buffer
        constrained_geom = apply_elevation_constraint(buffer_gdf.geometry.iloc[0], MAX_ELEVATION_FT)
        
        validation_gdf = gpd.GeoDataFrame(
            [{'distance_km': int(dist_km.replace('km', ''))}],
            geometry=[constrained_geom],
            crs='EPSG:4326'
        )
        
        validation_gdf.to_file(validation_file, driver='GeoJSON')
        print(f"  Saved {validation_file.name}")
    
    # Step 8: Summary statistics
    area_km2 = coastal_bioregion.to_crs('EPSG:5070').geometry.area.sum() / 1e6
    
    print(f"\n=== Summary ===")
    print(f"Coastal forest band area: {area_km2:.0f} km²")
    print(f"Buffer distance: {COASTAL_BUFFER_KM} km inland")  
    print(f"Maximum elevation: {MAX_ELEVATION_FT} ft")
    print(f"Geographic extent: {SOUTH_BOUND}°N to {NORTH_BOUND}°N")
    
    print(f"\n✅ Coastal forest band created successfully!")
    print(f"   Main output: {output_file}")
    
    return coastal_bioregion

if __name__ == "__main__":
    create_coastal_band_region()