"""
Extract black cottonwood distribution from LANDFIRE EVT 2023.

Fixed version: uses WGS84-defined target bounds converted to source Albers CRS
for windowed read, preventing the diagonal clipping artifact caused by using
rectangular Albers bounds directly.

Output: outputs/geojson/black_cottonwood_landfire_evt.geojson
"""

import numpy as np
import geopandas as gpd
import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.features import shapes
from rasterio.crs import CRS
from pyproj import Transformer
from shapely.geometry import shape
from shapely.ops import unary_union
import json
from pathlib import Path
import time

# --- Config ---
SOURCE_TIF = Path('data/raw/LF2023_EVT_CONUS/Tif/LF2023_EVT_CONUS.tif')
OUTPUT_GEOJSON = Path('outputs/geojson/black_cottonwood_landfire_evt.geojson')
OUTPUT_GEOJSON.parent.mkdir(parents=True, exist_ok=True)

# Western US target bounds in WGS84 (generous margin to capture all riparian areas)
TARGET_BOUNDS_WGS84 = (-125.5, 31.0, -101.5, 50.5)  # (min_lon, min_lat, max_lon, max_lat)
TARGET_CRS = CRS.from_epsg(4326)
TARGET_RES = 0.001  # ~111m at equator, ~80m at 45°N (similar to prior ~120m)

# LANDFIRE EVT codes for black cottonwood
# Primary: SAF 222 - Black Cottonwood-Willow (North Pacific riparian)
EVT_PRIMARY = {
    7156: 'North Pacific Lowland Riparian Forest',
    7157: 'North Pacific Lowland Riparian Forest-Hardwood',
    7158: 'North Pacific Montane Riparian Woodland',
    7159: 'North Pacific Montane Riparian Woodland-Hardwood',
    9165: 'North Pacific Shrub Swamp',
}
# Secondary: SAF 235 - Cottonwood-Willow (Rocky Mountain riparian)
EVT_SECONDARY = {
    9002: 'N. Rocky Mtn Foothill Riparian Woodland-Shrubland',
    9003: 'N. Rocky Mtn Lower Montane Riparian Woodland-Shrubland',
    9012: 'N. Rocky Mtn Lower Montane Riparian Woodland',
    9019: 'N. Rocky Mtn Lower Montane Riparian Shrubland',
    9502: 'N. Rocky Mtn Foothill Riparian Woodland-Shrubland (variant)',
    9503: 'N. Rocky Mtn Lower Montane Riparian Woodland-Shrubland (variant)',
    9512: 'N. Rocky Mtn Lower Montane Riparian Woodland (variant)',
    9519: 'N. Rocky Mtn Lower Montane Riparian Shrubland (variant)',
}
ALL_EVT_CODES = {**{k: 'primary' for k in EVT_PRIMARY}, **{k: 'secondary' for k in EVT_SECONDARY}}
MIN_AREA_KM2 = 1.0
SIMPLIFY_TOLERANCE = 0.003  # degrees


def wgs84_bounds_to_source_window(bounds_wgs84, src):
    """Convert WGS84 bounding box to a rasterio Window in source CRS."""
    min_lon, min_lat, max_lon, max_lat = bounds_wgs84
    transformer = Transformer.from_crs('EPSG:4326', src.crs, always_xy=True)

    # Transform all 4 corners to get proper bounds in projected CRS
    corners = [
        transformer.transform(min_lon, min_lat),
        transformer.transform(max_lon, min_lat),
        transformer.transform(min_lon, max_lat),
        transformer.transform(max_lon, max_lat),
    ]
    xs = [c[0] for c in corners]
    ys = [c[1] for c in corners]
    src_bounds = (min(xs), min(ys), max(xs), max(ys))

    return rasterio.windows.from_bounds(*src_bounds, transform=src.transform)


def extract_and_reproject():
    print(f'Opening source: {SOURCE_TIF}')
    t0 = time.time()

    with rasterio.open(SOURCE_TIF) as src:
        print(f'  Source CRS: {src.crs}')
        print(f'  Source shape: {src.shape}')

        # Compute windowed read from WGS84 target bounds
        window = wgs84_bounds_to_source_window(TARGET_BOUNDS_WGS84, src)
        window = window.intersection(rasterio.windows.Window(0, 0, src.width, src.height))
        print(f'  Reading window: {window}')

        src_data = src.read(1, window=window)
        src_transform = src.window_transform(window)
        src_crs = src.crs
        src_nodata = src.nodata

    print(f'  Read {src_data.shape} pixels in {time.time()-t0:.1f}s')

    # Reproject to WGS84
    print('Reprojecting to WGS84...')
    t1 = time.time()
    dst_transform, dst_width, dst_height = calculate_default_transform(
        src_crs, TARGET_CRS,
        src_data.shape[1], src_data.shape[0],
        *rasterio.windows.bounds(
            rasterio.windows.Window(0, 0, src_data.shape[1], src_data.shape[0]),
            src_transform
        ),
        resolution=TARGET_RES
    )
    dst_data = np.zeros((dst_height, dst_width), dtype=src_data.dtype)
    reproject(
        source=src_data,
        destination=dst_data,
        src_transform=src_transform,
        src_crs=src_crs,
        dst_transform=dst_transform,
        dst_crs=TARGET_CRS,
        resampling=Resampling.nearest,
        src_nodata=src_nodata,
        dst_nodata=0,
    )
    print(f'  Reprojected to {dst_data.shape} pixels in {time.time()-t1:.1f}s')

    return dst_data, dst_transform, TARGET_CRS


def apply_morphological_smoothing(mask):
    """Minimal 2x2 smoothing to reduce pixelation without expanding area."""
    from scipy.ndimage import binary_closing
    structure = np.ones((2, 2), dtype=bool)
    return binary_closing(mask, structure=structure, iterations=1)


def polygonize_evt(dst_data, dst_transform, dst_crs):
    print('Extracting cottonwood pixels...')
    features = []

    for evt_code, cottonwood_type in ALL_EVT_CODES.items():
        mask = (dst_data == evt_code).astype(np.uint8)
        pixel_count = mask.sum()
        if pixel_count == 0:
            continue
        print(f'  EVT {evt_code} ({cottonwood_type}): {pixel_count:,} pixels')

        # Minimal morphological smoothing
        mask_bool = apply_morphological_smoothing(mask.astype(bool))
        mask = mask_bool.astype(np.uint8)

        # Polygonize
        for geom, val in shapes(mask, mask=(mask == 1), transform=dst_transform):
            if val == 1:
                features.append({
                    'geometry': geom,
                    'cottonwood_type': cottonwood_type,
                    'evt_code': evt_code,
                })

    print(f'  Total raw polygons: {len(features):,}')
    return features


def build_geodataframe(features, crs):
    print('Building GeoDataFrame...')
    gdf = gpd.GeoDataFrame(
        [{'cottonwood_type': f['cottonwood_type'], 'evt_code': f['evt_code'],
          'geometry': shape(f['geometry'])} for f in features],
        crs=crs
    )

    # Reproject to equal-area for area calculation
    gdf_ea = gdf.to_crs('EPSG:5070')
    gdf['area_km2'] = gdf_ea.geometry.area / 1e6

    # Filter by minimum area
    before = len(gdf)
    gdf = gdf[gdf['area_km2'] >= MIN_AREA_KM2].copy()
    print(f'  After area filter (>={MIN_AREA_KM2} km²): {len(gdf):,} / {before:,}')

    # Add metadata
    gdf['species'] = 'Black Cottonwood'
    gdf['species_code'] = 747
    gdf['source'] = 'LANDFIRE EVT 2023'
    gdf['area_km2'] = gdf['area_km2'].round(2)

    # Simplify geometry
    gdf['geometry'] = gdf['geometry'].simplify(SIMPLIFY_TOLERANCE, preserve_topology=True)

    # Reorder columns
    gdf = gdf[['species', 'species_code', 'cottonwood_type', 'evt_code', 'area_km2', 'source', 'geometry']]

    return gdf


def main():
    t_total = time.time()

    dst_data, dst_transform, dst_crs = extract_and_reproject()
    features = polygonize_evt(dst_data, dst_transform, dst_crs)
    gdf = build_geodataframe(features, dst_crs)

    print(f'\nSummary:')
    print(f'  Total features: {len(gdf):,}')
    print(f'  Bounds: {gdf.total_bounds}')
    print(f'  By type:')
    for t, g in gdf.groupby('cottonwood_type'):
        print(f'    {t}: {len(g):,} polygons, {g["area_km2"].sum():.0f} km²')

    print(f'\nSaving to {OUTPUT_GEOJSON}...')
    gdf.to_file(OUTPUT_GEOJSON, driver='GeoJSON')

    size_mb = OUTPUT_GEOJSON.stat().st_size / 1e6
    print(f'  Saved {size_mb:.1f} MB in {time.time()-t_total:.1f}s total')


if __name__ == '__main__':
    main()
