#!/usr/bin/env python3
"""
Process PRISM precipitation normals for Pacific Northwest bioregion analysis.

This script extracts annual precipitation data from PRISM climate normals,
clips it to the Pacific Northwest extent, and creates both raster and vector
outputs for bioregion boundary refinement.

Key outputs:
- Annual precipitation raster aligned with existing PNW data
- Precipitation isolines for bioregion analysis
- Summary statistics for the region

Usage:
    python scripts/process_prism_precipitation_pnw.py
"""

import os
import sys
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.enums import Resampling
from rasterio.warp import reproject, calculate_default_transform
from rasterio.mask import mask
from rasterio.features import shapes
from shapely.geometry import shape
import matplotlib.pyplot as plt
from shapely.geometry import box
import json
from datetime import datetime

def get_pnw_extent_from_existing_raster():
    """Get the extent and CRS from existing PNW elevation data."""
    elev_path = "outputs/mapbox_masks/pnw_elevation_120m_mapbox.tif"
    
    if not os.path.exists(elev_path):
        raise FileNotFoundError(f"PNW elevation raster not found: {elev_path}")
    
    with rasterio.open(elev_path) as src:
        bounds = src.bounds
        crs = src.crs
        transform = src.transform
        shape = src.shape
        
    print(f"PNW extent from elevation data:")
    print(f"  Bounds: {bounds}")
    print(f"  CRS: {crs}")
    print(f"  Shape: {shape}")
    print(f"  Resolution: {abs(transform.a):.1f}m x {abs(transform.e):.1f}m")
    
    return bounds, crs, transform, shape

def load_prism_annual_precipitation():
    """Load the PRISM annual precipitation normals."""
    prism_path = "data/raw/PRISM_ppt_30yr_normal_800mM4_all_bil/PRISM_ppt_30yr_normal_800mM4_annual_bil.bil"
    
    if not os.path.exists(prism_path):
        raise FileNotFoundError(f"PRISM annual precipitation file not found: {prism_path}")
    
    with rasterio.open(prism_path) as src:
        precip_data = src.read(1)
        precip_transform = src.transform
        precip_crs = src.crs
        precip_nodata = src.nodata
        
        print(f"PRISM precipitation data:")
        print(f"  CRS: {precip_crs}")
        print(f"  Shape: {src.shape}")
        print(f"  Data range: {np.nanmin(precip_data):.1f} to {np.nanmax(precip_data):.1f} mm")
        print(f"  NoData value: {precip_nodata}")
        
        # Convert mm to inches for easier interpretation
        precip_inches = precip_data / 25.4
        print(f"  Data range (inches): {np.nanmin(precip_inches):.1f} to {np.nanmax(precip_inches):.1f}")
        
    return precip_data, precip_transform, precip_crs, precip_nodata

def create_pnw_boundary_from_extent(bounds, crs):
    """Create a boundary polygon from the PNW extent."""
    minx, miny, maxx, maxy = bounds
    boundary = box(minx, miny, maxx, maxy)
    
    boundary_gdf = gpd.GeoDataFrame([1], geometry=[boundary], crs=crs)
    
    # Buffer slightly to ensure we capture edge areas
    boundary_gdf_buffered = boundary_gdf.to_crs('EPSG:5070')  # Albers for buffering
    boundary_gdf_buffered.geometry = boundary_gdf_buffered.geometry.buffer(5000)  # 5km buffer
    boundary_gdf_buffered = boundary_gdf_buffered.to_crs(crs)  # Back to original CRS
    
    return boundary_gdf_buffered

def reproject_and_clip_precipitation(precip_data, precip_transform, precip_crs, precip_nodata, 
                                   target_bounds, target_crs, target_transform, target_shape):
    """Reproject and clip precipitation data to match PNW extent."""
    
    # Create output array
    output_precip = np.full(target_shape, np.nan, dtype=np.float32)
    
    # Reproject precipitation data to target CRS and extent
    reproject(
        source=precip_data,
        destination=output_precip,
        src_transform=precip_transform,
        src_crs=precip_crs,
        src_nodata=precip_nodata,
        dst_transform=target_transform,
        dst_crs=target_crs,
        dst_nodata=np.nan,
        resampling=Resampling.bilinear
    )
    
    # Mask invalid values
    output_precip = np.where(output_precip <= 0, np.nan, output_precip)
    
    return output_precip

def create_precipitation_isolines(precip_data, transform, crs, output_path):
    """Create precipitation isoline vectors for bioregion analysis."""
    
    # Convert to inches for more intuitive thresholds
    precip_inches = precip_data / 25.4
    
    # Define precipitation zones relevant to PNW bioregions
    precip_zones = {
        'very_dry': (0, 15),      # < 15 inches - Eastern OR/WA arid
        'dry': (15, 25),          # 15-25 inches - Eastern Cascades
        'moderate': (25, 50),     # 25-50 inches - Western valleys
        'wet': (50, 80),          # 50-80 inches - Western Cascades  
        'very_wet': (80, 150),    # 80-150 inches - Coastal forests
        'rainforest': (150, 300)  # >150 inches - Olympic rainforest
    }
    
    # Create categorical precipitation map
    precip_zones_array = np.full_like(precip_inches, 0, dtype=np.uint8)
    
    zone_codes = {}
    for i, (zone_name, (min_val, max_val)) in enumerate(precip_zones.items(), 1):
        mask = (precip_inches >= min_val) & (precip_inches < max_val)
        precip_zones_array[mask] = i
        zone_codes[i] = zone_name
        
        pixel_count = np.sum(mask)
        pixel_area_m2 = abs(transform.a) * abs(transform.e)
        area_km2 = pixel_count * pixel_area_m2 / 1e6
        print(f"  {zone_name}: {min_val}-{max_val} inches, {pixel_count:,} pixels, {area_km2:.1f} km²")
    
    # Convert zones to vector polygons
    zone_shapes = []
    for geom, value in shapes(precip_zones_array, mask=precip_zones_array > 0, transform=transform):
        if value in zone_codes:
            zone_shapes.append({
                'geometry': shape(geom),  # Convert dict to Shapely geometry
                'zone_code': value,
                'zone_name': zone_codes[value],
                'precip_range': f"{precip_zones[zone_codes[value]][0]}-{precip_zones[zone_codes[value]][1]} inches"
            })
    
    # Create GeoDataFrame
    if zone_shapes:
        # Separate geometry from other attributes
        geometries = [shape['geometry'] for shape in zone_shapes]
        attributes = [{k: v for k, v in shape.items() if k != 'geometry'} for shape in zone_shapes]
        
        precip_zones_gdf = gpd.GeoDataFrame(attributes, geometry=geometries, crs=crs)
    else:
        # Create empty GeoDataFrame with proper schema
        precip_zones_gdf = gpd.GeoDataFrame(
            columns=['zone_code', 'zone_name', 'precip_range'],
            geometry=[],
            crs=crs
        )
    
    # Calculate areas
    precip_zones_gdf_albers = precip_zones_gdf.to_crs('EPSG:5070')
    precip_zones_gdf['area_km2'] = precip_zones_gdf_albers.geometry.area / 1e6
    
    # Save zones
    precip_zones_gdf.to_file(output_path, driver='GeoJSON')
    print(f"\nPrecipitation zones saved to: {output_path}")
    
    return precip_zones_gdf

def save_precipitation_raster(precip_data, transform, crs, output_path):
    """Save the processed precipitation raster."""
    
    with rasterio.open(
        output_path,
        'w',
        driver='GTiff',
        height=precip_data.shape[0],
        width=precip_data.shape[1],
        count=1,
        dtype=precip_data.dtype,
        crs=crs,
        transform=transform,
        nodata=np.nan,
        compress='lzw'
    ) as dst:
        dst.write(precip_data, 1)
        
        # Add metadata
        dst.update_tags(
            source='PRISM Climate Group',
            description='Annual precipitation normals (1981-2010)',
            units='millimeters',
            processing_date=datetime.now().isoformat()
        )
    
    print(f"Precipitation raster saved to: {output_path}")

def create_summary_statistics(precip_data, precip_zones_gdf, output_path):
    """Create summary statistics for the precipitation data."""
    
    # Convert to inches for summary
    precip_inches = precip_data / 25.4
    valid_precip = precip_inches[~np.isnan(precip_inches)]
    
    summary = {
        'processing_info': {
            'source': 'PRISM Climate Group - 30-year normals (1981-2010)',
            'processing_date': datetime.now().isoformat(),
            'extent': 'Pacific Northwest (WA, OR, ID)',
            'units': 'millimeters (mm)',
            'conversion_note': 'Multiply by 0.0394 to convert to inches'
        },
        'statistics': {
            'total_pixels': int(np.prod(precip_data.shape)),
            'valid_pixels': int(len(valid_precip)),
            'coverage_percent': float(len(valid_precip) / np.prod(precip_data.shape) * 100),
            'min_precip_mm': float(np.nanmin(precip_data)),
            'max_precip_mm': float(np.nanmax(precip_data)),
            'mean_precip_mm': float(np.nanmean(precip_data)),
            'median_precip_mm': float(np.nanmedian(precip_data)),
            'min_precip_inches': float(np.nanmin(precip_inches)),
            'max_precip_inches': float(np.nanmax(precip_inches)),
            'mean_precip_inches': float(np.nanmean(precip_inches)),
            'median_precip_inches': float(np.nanmedian(precip_inches))
        },
        'bioregion_relevance': {
            'coastal_rainforest': '>150 inches (Olympic Peninsula)',
            'wet_coastal': '80-150 inches (Coast Range, Western Cascades)',
            'western_valleys': '25-50 inches (Willamette, Puget Sound)',
            'eastern_cascades': '15-25 inches (Pine/fir transition)',
            'rain_shadow': '<15 inches (Eastern WA/OR)',
            'olympic_rain_shadow': 'Notable gradient from 150+ to <20 inches'
        },
        'precipitation_zones': []
    }
    
    # Add zone statistics
    for _, zone in precip_zones_gdf.iterrows():
        summary['precipitation_zones'].append({
            'zone_name': zone['zone_name'],
            'precipitation_range': zone['precip_range'],
            'area_km2': float(zone['area_km2'])
        })
    
    # Save summary
    with open(output_path, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"Summary statistics saved to: {output_path}")
    return summary

def main():
    """Main processing function."""
    print("Processing PRISM precipitation data for Pacific Northwest bioregions...")
    print("=" * 70)
    
    # Create output directory
    os.makedirs("outputs/climate", exist_ok=True)
    
    # Get PNW extent from existing data
    print("\n1. Loading PNW extent from existing elevation data...")
    pnw_bounds, pnw_crs, pnw_transform, pnw_shape = get_pnw_extent_from_existing_raster()
    
    # Load PRISM precipitation data
    print("\n2. Loading PRISM annual precipitation normals...")
    precip_data, precip_transform, precip_crs, precip_nodata = load_prism_annual_precipitation()
    
    # Reproject and clip precipitation data
    print("\n3. Reprojecting and clipping precipitation to PNW extent...")
    pnw_precip = reproject_and_clip_precipitation(
        precip_data, precip_transform, precip_crs, precip_nodata,
        pnw_bounds, pnw_crs, pnw_transform, pnw_shape
    )
    
    print(f"   Clipped precipitation data:")
    print(f"   Range: {np.nanmin(pnw_precip):.1f} to {np.nanmax(pnw_precip):.1f} mm")
    print(f"   Range: {np.nanmin(pnw_precip)/25.4:.1f} to {np.nanmax(pnw_precip)/25.4:.1f} inches")
    print(f"   Valid pixels: {np.sum(~np.isnan(pnw_precip)):,} / {np.prod(pnw_shape):,}")
    
    # Save precipitation raster
    print("\n4. Saving precipitation raster...")
    raster_output = "outputs/climate/pnw_precipitation_annual_normals.tif"
    save_precipitation_raster(pnw_precip, pnw_transform, pnw_crs, raster_output)
    
    # Create precipitation zones
    print("\n5. Creating precipitation zones for bioregion analysis...")
    zones_output = "outputs/climate/pnw_precipitation_zones.geojson"
    precip_zones_gdf = create_precipitation_isolines(pnw_precip, pnw_transform, pnw_crs, zones_output)
    
    # Create summary statistics
    print("\n6. Creating summary statistics...")
    summary_output = "outputs/climate/pnw_precipitation_summary.json"
    summary = create_summary_statistics(pnw_precip, precip_zones_gdf, summary_output)
    
    print("\n" + "=" * 70)
    print("PRISM precipitation processing completed successfully!")
    print("\nOutput files created:")
    print(f"  - Precipitation raster: {raster_output}")
    print(f"  - Precipitation zones: {zones_output}")
    print(f"  - Summary statistics: {summary_output}")
    
    print(f"\nPrecipitation summary:")
    print(f"  - Range: {summary['statistics']['min_precip_inches']:.1f} to {summary['statistics']['max_precip_inches']:.1f} inches")
    print(f"  - Mean: {summary['statistics']['mean_precip_inches']:.1f} inches")
    print(f"  - Perfect for distinguishing PNW bioregions!")

if __name__ == "__main__":
    main()