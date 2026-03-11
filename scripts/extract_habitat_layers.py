"""
Bulk LANDFIRE EVT habitat layer extraction for Forayz mushroom foraging map.

Reads the LANDFIRE 2023 EVT raster ONCE for western US, reprojects ONCE,
then extracts all habitat layers in a single pass.

Outputs (per layer):
  outputs/geojson/{layer_id}.geojson
  outputs/{layer_id}.pmtiles

Usage:
  python scripts/extract_habitat_layers.py                 # all layers
  python scripts/extract_habitat_layers.py --layer garry-oak
  python scripts/extract_habitat_layers.py --skip-existing  # skip if GeoJSON exists
  python scripts/extract_habitat_layers.py --no-upload      # skip R2 upload

Alder gap note:
  LANDFIRE EVT has no red alder (Alnus rubra) type for the lower 48.
  Alder only appears in Alaska EVT classes. Consider FIA plot data (SPCD 351)
  as an alternative for a point-based alder layer.
"""

import argparse
import subprocess
import sys
import time
from pathlib import Path

import geopandas as gpd
import numpy as np
import rasterio
from pyproj import Transformer
from rasterio.crs import CRS
from rasterio.features import shapes
from rasterio.warp import Resampling, calculate_default_transform, reproject
from scipy.ndimage import binary_closing
from shapely.geometry import shape

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SOURCE_TIF = Path("data/raw/LF2023_EVT_CONUS/Tif/LF2023_EVT_CONUS.tif")
GEOJSON_DIR = Path("outputs/geojson")
PMTILES_DIR = Path("outputs")
GEOJSON_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Raster extraction config
# ---------------------------------------------------------------------------
# Western US bounds — generous margin, WGS84
TARGET_BOUNDS_WGS84 = (-125.5, 31.0, -101.5, 50.5)
TARGET_CRS = CRS.from_epsg(4326)
TARGET_RES = 0.001  # ~111m lat, ~80m lon at 45°N

# ---------------------------------------------------------------------------
# Layer definitions
# ---------------------------------------------------------------------------
# Each layer: display_name, description, evt_codes {code: label}, min_area_km2,
#             r2_key, tippecanoe_layer (internal MVT layer name)
LAYERS = {
    # --- Already extracted but included for completeness / re-run capability ---
    "black-cottonwood": {
        "display_name": "Riparian Cottonwood",
        "description": "Cottonwood-willow riparian corridors along streams and rivers",
        "evt_codes": {
            7156: "North Pacific Lowland Riparian Forest",
            7157: "North Pacific Lowland Riparian Shrubland",
            7158: "North Pacific Montane Riparian Woodland",
            7159: "North Pacific Montane Riparian Shrubland",
            9165: "North Pacific Shrub Swamp",
            9002: "Columbia Basin Foothill Riparian Woodland",
            9003: "Great Basin Foothill and Lower Montane Riparian Woodland",
            9012: "Northern Rocky Mountain Lower Montane Riparian Woodland",
            9019: "Rocky Mountain Lower Montane-Foothill Riparian Woodland",
            9502: "Columbia Basin Foothill Riparian Shrubland",
            9503: "Great Basin Foothill and Lower Montane Riparian Shrubland",
            9512: "Northern Rocky Mountain Lower Montane Riparian Shrubland",
            9519: "Rocky Mountain Lower Montane-Foothill Riparian Shrubland",
        },
        "min_area_km2": 1.0,
        "r2_key": "trees/black_cottonwood.pmtiles",
        "tippecanoe_layer": "cottonwood",
    },

    # --- New layers ---
    "coastal-rainforest": {
        "display_name": "Coastal Rainforest",
        "description": "Sitka spruce, western red cedar, and wet hemlock forest along the Pacific coast",
        "evt_codes": {
            7036: "North Pacific Seasonal Sitka Spruce Forest",
            7039: "North Pacific Maritime Mesic-Wet Douglas-fir-Western Hemlock Forest",
            7178: "North Pacific Hypermaritime Western Red-cedar-Western Hemlock Forest",
        },
        "min_area_km2": 1.0,
        "r2_key": "trees/coastal_rainforest.pmtiles",
        "tippecanoe_layer": "coastal_rainforest",
    },

    "west-cascades-conifer": {
        "display_name": "West Cascades Conifer Forest",
        "description": "Maritime Douglas-fir, western hemlock, and silver fir — the primary chanterelle belt",
        "evt_codes": {
            7035: "North Pacific Dry Douglas-fir-(Madrone) Forest and Woodland",
            7037: "North Pacific Maritime Dry-Mesic Douglas-fir-Western Hemlock Forest",
            7042: "North Pacific Mesic Western Hemlock-Silver Fir Forest",
            7174: "North Pacific Dry-Mesic Silver Fir-Western Hemlock-Douglas-fir Forest",
        },
        "min_area_km2": 1.0,
        "r2_key": "trees/west_cascades_conifer.pmtiles",
        "tippecanoe_layer": "west_cascades_conifer",
    },

    "high-cascades-fir": {
        "display_name": "High Cascades & Sierra Fir Zone",
        "description": "Mountain hemlock, silver fir, red fir, and white fir at upper elevations",
        "evt_codes": {
            7038: "North Pacific Maritime Mesic Subalpine Parkland",
            7041: "North Pacific Mountain Hemlock Forest",
            7032: "Mediterranean California Red Fir Forest",
            7172: "Sierran-Intermontane Desert Western White Pine-White Fir Woodland",
        },
        "min_area_km2": 1.0,
        "r2_key": "trees/high_cascades_fir.pmtiles",
        "tippecanoe_layer": "high_cascades_fir",
    },

    "rocky-mtn-douglas-fir": {
        "display_name": "Rocky Mountain Douglas-fir Forest",
        "description": "Interior montane Douglas-fir forests of the Rocky Mountains",
        "evt_codes": {
            7166: "Middle Rocky Mountain Montane Douglas-fir Forest and Woodland",
        },
        "min_area_km2": 1.0,
        "r2_key": "trees/rocky_mtn_douglas_fir.pmtiles",
        "tippecanoe_layer": "rocky_mtn_douglas_fir",
    },

    "rocky-mtn-spruce-fir": {
        "display_name": "Rocky Mountain Spruce-Fir",
        "description": "Subalpine Engelmann spruce and subalpine fir — prime porcini/king bolete habitat",
        "evt_codes": {
            7055: "Rocky Mountain Subalpine Dry-Mesic Spruce-Fir Forest and Woodland",
            7056: "Rocky Mountain Subalpine Mesic-Wet Spruce-Fir Forest and Woodland",
        },
        "min_area_km2": 2.0,
        "r2_key": "trees/rocky_mtn_spruce_fir.pmtiles",
        "tippecanoe_layer": "rocky_mtn_spruce_fir",
    },

    "ponderosa-pine": {
        "display_name": "Ponderosa Pine Woodland",
        "description": "Open ponderosa and Jeffrey pine forests — east Cascades, Rockies, and Sierra Nevada",
        "evt_codes": {
            7031: "California Montane Jeffrey Pine-(Ponderosa Pine) Woodland",
            7053: "Northern Rocky Mountain Ponderosa Pine Woodland and Savanna",
            7054: "Southern Rocky Mountain Ponderosa Pine Woodland",
            7060: "East Cascades Ponderosa Pine Forest and Woodland",
            7179: "Northwestern Great Plains-Black Hills Ponderosa Pine Woodland and Savanna",
        },
        "min_area_km2": 2.0,
        "r2_key": "trees/ponderosa_pine.pmtiles",
        "tippecanoe_layer": "ponderosa_pine",
    },

    "lodgepole-pine": {
        "display_name": "Lodgepole Pine Forest",
        "description": "Lodgepole pine forests of the Rockies and Sierra — primary matsutake habitat",
        "evt_codes": {
            7050: "Rocky Mountain Lodgepole Pine Forest",
            7058: "Sierra Nevada Subalpine Lodgepole Pine Forest and Woodland",
            7167: "Rocky Mountain Poor-Site Lodgepole Pine Forest",
        },
        "min_area_km2": 5.0,  # Higher threshold — massive footprint
        "r2_key": "trees/lodgepole_pine.pmtiles",
        "tippecanoe_layer": "lodgepole_pine",
    },

    "garry-oak": {
        "display_name": "Garry Oak Woodland",
        "description": "Oregon white oak (Garry oak) woodlands — rare PNW habitat, spring morels, chanterelles",
        "evt_codes": {
            7008: "North Pacific Oak Woodland",
            7262: "East Cascades Oak Forest and Woodland",
            7263: "East Cascades Oak-Ponderosa Pine Forest and Woodland",
        },
        "min_area_km2": 0.5,  # Lower threshold — these are fragmented and rare
        "r2_key": "trees/garry_oak.pmtiles",
        "tippecanoe_layer": "garry_oak",
    },

    "california-oak": {
        "display_name": "California Oak Woodland",
        "description": "Valley oaks, blue oaks, black oaks, and live oaks of California",
        "evt_codes": {
            7029: "Mediterranean California Mixed Oak Woodland",
            7112: "California Central Valley Mixed Oak Savanna",
            7113: "California Coastal Live Oak Woodland and Savanna",
            7118: "Southern California Oak Woodland and Savanna",
            7260: "Mediterranean California Lower Montane Black Oak Forest and Woodland",
            7261: "Mediterranean California Lower Montane Black Oak-Conifer Forest and Woodland",
            7264: "California Lower Montane Blue Oak Woodland and Savanna",
            7265: "California Lower Montane Blue Oak-Foothill Pine Woodland and Savanna",
        },
        "min_area_km2": 0.5,
        "r2_key": "trees/california_oak.pmtiles",
        "tippecanoe_layer": "california_oak",
    },
}


# ---------------------------------------------------------------------------
# Raster helpers
# ---------------------------------------------------------------------------

def wgs84_bounds_to_window(bounds_wgs84, src):
    """Convert WGS84 bounding box to a rasterio Window in source CRS.
    Converts all 4 corners to handle non-rectangular projections correctly."""
    min_lon, min_lat, max_lon, max_lat = bounds_wgs84
    transformer = Transformer.from_crs("EPSG:4326", src.crs, always_xy=True)
    corners = [
        transformer.transform(min_lon, min_lat),
        transformer.transform(max_lon, min_lat),
        transformer.transform(min_lon, max_lat),
        transformer.transform(max_lon, max_lat),
    ]
    xs = [c[0] for c in corners]
    ys = [c[1] for c in corners]
    src_bounds = (min(xs), min(ys), max(xs), max(ys))
    win = rasterio.windows.from_bounds(*src_bounds, transform=src.transform)
    return win.intersection(rasterio.windows.Window(0, 0, src.width, src.height))


def load_western_us_raster():
    """Read LANDFIRE EVT raster for western US and reproject to WGS84.
    Returns (data_array, transform, crs)."""
    print(f"Opening source: {SOURCE_TIF}")
    t0 = time.time()

    with rasterio.open(SOURCE_TIF) as src:
        print(f"  Source CRS: {src.crs}, shape: {src.shape}")
        window = wgs84_bounds_to_window(TARGET_BOUNDS_WGS84, src)
        print(f"  Reading window: col={int(window.col_off)}, row={int(window.row_off)}, "
              f"w={int(window.width)}, h={int(window.height)}")
        src_data = src.read(1, window=window)
        src_transform = src.window_transform(window)
        src_crs = src.crs
        src_nodata = src.nodata

    print(f"  Read {src_data.shape} pixels in {time.time() - t0:.1f}s")

    print("Reprojecting to WGS84...")
    t1 = time.time()
    dst_transform, dst_width, dst_height = calculate_default_transform(
        src_crs, TARGET_CRS,
        src_data.shape[1], src_data.shape[0],
        *rasterio.windows.bounds(
            rasterio.windows.Window(0, 0, src_data.shape[1], src_data.shape[0]),
            src_transform,
        ),
        resolution=TARGET_RES,
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
    print(f"  Reprojected to {dst_data.shape} pixels in {time.time() - t1:.1f}s")
    return dst_data, dst_transform, TARGET_CRS


# ---------------------------------------------------------------------------
# Per-layer extraction
# ---------------------------------------------------------------------------

def extract_layer(layer_id, layer_cfg, raster_data, raster_transform, raster_crs):
    """Extract one habitat layer from the pre-loaded raster."""
    evt_codes = layer_cfg["evt_codes"]
    min_area = layer_cfg["min_area_km2"]

    print(f"\n{'='*60}")
    print(f"Layer: {layer_cfg['display_name']} ({layer_id})")
    print(f"  EVT codes: {list(evt_codes.keys())}")
    t0 = time.time()

    features = []
    for evt_code, evt_label in evt_codes.items():
        mask = (raster_data == evt_code).astype(np.uint8)
        px = int(mask.sum())
        if px == 0:
            print(f"  {evt_code}: 0 pixels — skipping")
            continue
        print(f"  {evt_code} ({evt_label[:50]}): {px:,} pixels")

        # Minimal morphological smoothing — 2×2, 1 iteration
        smooth = binary_closing(mask.astype(bool), structure=np.ones((2, 2)), iterations=1).astype(np.uint8)

        for geom, val in shapes(smooth, mask=(smooth == 1), transform=raster_transform):
            if val == 1:
                features.append({"geometry": geom, "evt_code": evt_code, "evt_label": evt_label})

    if not features:
        print(f"  No features found — skipping layer")
        return None

    print(f"  Raw polygons: {len(features):,}")

    # Build GeoDataFrame
    gdf = gpd.GeoDataFrame(
        [{"evt_code": f["evt_code"], "evt_label": f["evt_label"],
          "geometry": shape(f["geometry"])} for f in features],
        crs=raster_crs,
    )

    # Area filter (equal-area projection)
    gdf_ea = gdf.to_crs("EPSG:5070")
    gdf["area_km2"] = (gdf_ea.geometry.area / 1e6).round(2)
    before = len(gdf)
    gdf = gdf[gdf["area_km2"] >= min_area].copy()
    print(f"  After area filter (>={min_area} km²): {len(gdf):,} / {before:,} polygons")

    # Add metadata
    gdf["habitat"] = layer_cfg["display_name"]
    gdf["source"] = "LANDFIRE EVT 2023"
    gdf["geometry"] = gdf["geometry"].simplify(0.003, preserve_topology=True)
    gdf = gdf[["habitat", "evt_code", "evt_label", "area_km2", "source", "geometry"]]

    # Save GeoJSON
    out_path = GEOJSON_DIR / f"{layer_id}.geojson"
    gdf.to_file(out_path, driver="GeoJSON")
    size_mb = out_path.stat().st_size / 1e6
    print(f"  Saved GeoJSON: {out_path} ({size_mb:.1f} MB, {len(gdf):,} features, "
          f"bounds={[round(x, 2) for x in gdf.total_bounds]}) in {time.time()-t0:.1f}s")
    return out_path


def build_pmtiles(layer_id, layer_cfg, geojson_path):
    """Run tippecanoe to produce PMTiles from GeoJSON."""
    # Use underscore filename to match R2 key convention (e.g. black_cottonwood.pmtiles)
    out_path = PMTILES_DIR / f"{layer_id.replace('-', '_')}.pmtiles"
    cmd = [
        "tippecanoe",
        "-o", str(out_path),
        "-z", "12", "-Z", "4",
        "--drop-densest-as-needed",
        "--simplification=10",
        "-l", layer_cfg["tippecanoe_layer"],
        "--force",
        str(geojson_path),
    ]
    print(f"  Running tippecanoe...")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  ERROR: tippecanoe failed:\n{result.stderr[-500:]}")
        return None
    size_mb = out_path.stat().st_size / 1e6
    print(f"  PMTiles: {out_path} ({size_mb:.1f} MB)")
    return out_path


def upload_to_r2(layer_id, layer_cfg, pmtiles_path):
    """Upload PMTiles to R2 via the morel-maps R2Uploader."""
    morel_pipeline = Path(__file__).parent.parent.parent / "morel-maps/scripts/pipeline_v2"
    sys.path.insert(0, str(morel_pipeline))

    import os
    import yaml

    with open(morel_pipeline / "config.yaml") as f:
        config = yaml.safe_load(f)

    creds_file = Path.home() / ".config/morel-maps/r2_credentials"
    with open(creds_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith("export "):
                line = line[7:]
            if "=" in line and not line.startswith("#"):
                k, v = line.split("=", 1)
                os.environ[k] = v.strip('"').strip("'")

    from output.r2_uploader import R2Uploader
    uploader = R2Uploader(config)

    r2_key = layer_cfg["r2_key"]
    layer_name = r2_key.replace(".pmtiles", "")  # strip extension for upload_pmtiles
    result = uploader.upload_pmtiles(pmtiles_path=str(pmtiles_path), layer_name=layer_name, force=True)
    if result.success:
        print(f"  Uploaded: {result.url} ({result.size_mb:.1f} MB)")
    else:
        print(f"  Upload failed: {result.error}")
    return result.success


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Extract LANDFIRE habitat layers for Forayz")
    parser.add_argument("--layer", help="Extract only this layer ID")
    parser.add_argument("--skip-existing", action="store_true",
                        help="Skip layers where GeoJSON already exists")
    parser.add_argument("--no-upload", action="store_true", help="Skip R2 upload")
    parser.add_argument("--geojson-only", action="store_true", help="Skip tippecanoe and upload")
    args = parser.parse_args()

    # Select layers to process
    if args.layer:
        if args.layer not in LAYERS:
            print(f"Unknown layer: {args.layer}. Available: {list(LAYERS.keys())}")
            sys.exit(1)
        selected = {args.layer: LAYERS[args.layer]}
    else:
        selected = LAYERS

    # Filter skip-existing
    if args.skip_existing:
        to_run = {}
        for lid, lcfg in selected.items():
            geojson = GEOJSON_DIR / f"{lid}.geojson"
            if geojson.exists():
                print(f"Skipping {lid} — GeoJSON exists ({geojson.stat().st_size/1e6:.1f} MB)")
            else:
                to_run[lid] = lcfg
        selected = to_run

    if not selected:
        print("Nothing to do.")
        return

    print(f"\nLayers to extract: {list(selected.keys())}")
    print(f"Source: {SOURCE_TIF}\n")

    t_total = time.time()

    # Load raster ONCE for all layers
    raster_data, raster_transform, raster_crs = load_western_us_raster()

    # Extract each layer
    geojson_paths = {}
    for layer_id, layer_cfg in selected.items():
        result = extract_layer(layer_id, layer_cfg, raster_data, raster_transform, raster_crs)
        if result:
            geojson_paths[layer_id] = result

    if args.geojson_only:
        print(f"\nDone (GeoJSON only). Total: {time.time()-t_total:.1f}s")
        return

    # Build PMTiles
    pmtiles_paths = {}
    print(f"\n{'='*60}")
    print("Building PMTiles...")
    for layer_id, geojson_path in geojson_paths.items():
        layer_cfg = selected[layer_id]
        print(f"\n{layer_cfg['display_name']}")
        result = build_pmtiles(layer_id, layer_cfg, geojson_path)
        if result:
            pmtiles_paths[layer_id] = result

    if args.no_upload:
        print(f"\nDone (no upload). Total: {time.time()-t_total:.1f}s")
        return

    # Upload to R2
    print(f"\n{'='*60}")
    print("Uploading to R2...")
    for layer_id, pmtiles_path in pmtiles_paths.items():
        layer_cfg = selected[layer_id]
        print(f"\n{layer_cfg['display_name']}")
        upload_to_r2(layer_id, layer_cfg, pmtiles_path)

    print(f"\n{'='*60}")
    print(f"All done. Total time: {time.time()-t_total:.1f}s")
    print(f"\nNext step: update Cloudflare Worker with new layer routes, then deploy.")
    print("  cd ~/cloudflare-workers/forayz-tiles && npx wrangler deploy")


if __name__ == "__main__":
    main()
