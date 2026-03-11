# Habitat Layers — Forayz Mushroom Foraging Map

*Created: March 2026*

These are ecosystem/habitat layers for the Forayz mushroom foraging app, covering the Western US.
The layers are intentionally habitat-focused (not mushroom-specific) — foragers decide what
tree associations they care about.

## Live Layers

All served via: `https://tiles-api.salishmushrooms.com/tilejson/{layer-id}.json`

| Layer ID | Display Name | PMTiles | Size | Features | Notes |
|----------|-------------|---------|------|----------|-------|
| `black-cottonwood` | Riparian Cottonwood | `trees/black_cottonwood.pmtiles` | 1.4MB | 1,402 | Streams, river corridors |
| `coastal-rainforest` | Coastal Rainforest | `trees/coastal_rainforest.pmtiles` | 4.1MB | 2,124 | Sitka spruce, red cedar, wet hemlock |
| `west-cascades-conifer` | West Cascades Conifer Forest | `trees/west_cascades_conifer.pmtiles` | 8.1MB | 3,798 | Primary chanterelle belt |
| `high-cascades-fir` | High Cascades & Sierra Fir Zone | `trees/high_cascades_fir.pmtiles` | 2.4MB | 1,476 | Mountain hemlock, red fir, white fir |
| `rocky-mtn-douglas-fir` | Rocky Mountain Douglas-fir Forest | `trees/rocky_mtn_douglas_fir.pmtiles` | 5.1MB | 2,342 | Interior montane Rockies |
| `rocky-mtn-spruce-fir` | Rocky Mountain Spruce-Fir | `trees/rocky_mtn_spruce_fir.pmtiles` | 9.2MB | 4,230 | Subalpine porcini/king bolete belt |
| `ponderosa-pine` | Ponderosa Pine Woodland | `trees/ponderosa_pine.pmtiles` | 11.8MB | 3,049 | Matsutake, morels post-burn |
| `lodgepole-pine` | Lodgepole Pine Forest | `trees/lodgepole_pine.pmtiles` | 5.3MB | 944 | Primary matsutake habitat |
| `garry-oak` | Garry Oak Woodland | `trees/garry_oak.pmtiles` | 0.2MB | 319 | Rare PNW habitat, fragmented |
| `california-oak` | California Oak Woodland | `trees/california_oak.pmtiles` | 1.8MB | 2,568 | Valley, blue, black, live oaks |

**Total R2 storage**: ~49 MB across 10 layers

## Data Source

**LANDFIRE 2023 Existing Vegetation Type (EVT)**
- File: `data/raw/LF2023_EVT_CONUS/Tif/LF2023_EVT_CONUS.tif` (2.8 GB, EPSG:5070, 30m)
- Lookup: `data/raw/LF2024_EVT.csv` (1,069 vegetation types)
- Coverage: CONUS (we extract Western US only)

## EVT Codes by Layer

### black-cottonwood
| Code | Name |
|------|------|
| 7156 | North Pacific Lowland Riparian Forest |
| 7157 | North Pacific Lowland Riparian Shrubland |
| 7158 | North Pacific Montane Riparian Woodland |
| 7159 | North Pacific Montane Riparian Shrubland |
| 9165 | North Pacific Shrub Swamp |
| 9002 | Columbia Basin Foothill Riparian Woodland |
| 9003 | Great Basin Foothill and Lower Montane Riparian Woodland |
| 9012 | Northern Rocky Mountain Lower Montane Riparian Woodland |
| 9019 | Rocky Mountain Lower Montane-Foothill Riparian Woodland |
| 9502 | Columbia Basin Foothill Riparian Shrubland |
| 9503 | Great Basin Foothill and Lower Montane Riparian Shrubland |
| 9512 | Northern Rocky Mountain Lower Montane Riparian Shrubland |
| 9519 | Rocky Mountain Lower Montane-Foothill Riparian Shrubland |

### coastal-rainforest
| Code | Name |
|------|------|
| 7036 | North Pacific Seasonal Sitka Spruce Forest |
| 7039 | North Pacific Maritime Mesic-Wet Douglas-fir-Western Hemlock Forest |
| 7178 | North Pacific Hypermaritime Western Red-cedar-Western Hemlock Forest |

### west-cascades-conifer
| Code | Name |
|------|------|
| 7035 | North Pacific Dry Douglas-fir-(Madrone) Forest and Woodland |
| 7037 | North Pacific Maritime Dry-Mesic Douglas-fir-Western Hemlock Forest |
| 7042 | North Pacific Mesic Western Hemlock-Silver Fir Forest |
| 7174 | North Pacific Dry-Mesic Silver Fir-Western Hemlock-Douglas-fir Forest |

### high-cascades-fir
| Code | Name |
|------|------|
| 7038 | North Pacific Maritime Mesic Subalpine Parkland |
| 7041 | North Pacific Mountain Hemlock Forest |
| 7032 | Mediterranean California Red Fir Forest |
| 7172 | Sierran-Intermontane Desert Western White Pine-White Fir Woodland |

### rocky-mtn-douglas-fir
| Code | Name |
|------|------|
| 7166 | Middle Rocky Mountain Montane Douglas-fir Forest and Woodland |

### rocky-mtn-spruce-fir
| Code | Name |
|------|------|
| 7055 | Rocky Mountain Subalpine Dry-Mesic Spruce-Fir Forest and Woodland |
| 7056 | Rocky Mountain Subalpine Mesic-Wet Spruce-Fir Forest and Woodland |

### ponderosa-pine
| Code | Name |
|------|------|
| 7031 | California Montane Jeffrey Pine-(Ponderosa Pine) Woodland |
| 7053 | Northern Rocky Mountain Ponderosa Pine Woodland and Savanna |
| 7054 | Southern Rocky Mountain Ponderosa Pine Woodland |
| 7060 | East Cascades Ponderosa Pine Forest and Woodland |
| 7179 | Northwestern Great Plains-Black Hills Ponderosa Pine Woodland and Savanna |

### lodgepole-pine
| Code | Name |
|------|------|
| 7050 | Rocky Mountain Lodgepole Pine Forest |
| 7058 | Sierra Nevada Subalpine Lodgepole Pine Forest and Woodland |
| 7167 | Rocky Mountain Poor-Site Lodgepole Pine Forest |

### garry-oak
| Code | Name |
|------|------|
| 7008 | North Pacific Oak Woodland |
| 7262 | East Cascades Oak Forest and Woodland |
| 7263 | East Cascades Oak-Ponderosa Pine Forest and Woodland |

### california-oak
| Code | Name |
|------|------|
| 7029 | Mediterranean California Mixed Oak Woodland |
| 7112 | California Central Valley Mixed Oak Savanna |
| 7113 | California Coastal Live Oak Woodland and Savanna |
| 7118 | Southern California Oak Woodland and Savanna |
| 7260 | Mediterranean California Lower Montane Black Oak Forest and Woodland |
| 7261 | Mediterranean California Lower Montane Black Oak-Conifer Forest and Woodland |
| 7264 | California Lower Montane Blue Oak Woodland and Savanna |
| 7265 | California Lower Montane Blue Oak-Foothill Pine Woodland and Savanna |

## Processing Pipeline

```
LANDFIRE EVT TIF (EPSG:5070, 30m, 2.8 GB)
  → Windowed read: western US (WGS84 bounds → Albers window, ~90s)
  → Reproject to WGS84 at 0.001° (~111m) resolution (~42s)
  → Per layer: mask EVT codes → morphological smooth (2×2, 1 iter)
             → polygonize → filter by min_area_km2 → simplify (0.003°)
             → save GeoJSON
  → Tippecanoe: z4-z12, --drop-densest-as-needed, --simplification=10
  → Upload to R2: trees/{layer_name}.pmtiles
```

**Key fix**: WGS84 bounds are converted to Albers CRS before computing the
raster window. Without this, the output has a diagonal rectangular clipping
artifact (parallelogram in WGS84 space).

**Total runtime**: ~27 minutes for all 10 layers (single raster read/reproject)

**Extraction script**: `scripts/extract_habitat_layers.py`
- Re-run any single layer: `python scripts/extract_habitat_layers.py --layer garry-oak`
- Re-run all, skip existing: `python scripts/extract_habitat_layers.py --skip-existing`
- GeoJSON only (no tippecanoe/upload): add `--geojson-only`

## Cloudflare Worker

Layer routes defined in `~/cloudflare-workers/forayz-tiles/src/index.js`.
Deploy after any worker changes: `cd ~/cloudflare-workers/forayz-tiles && npx wrangler deploy`

## Known Gaps

### Red Alder (Alnus rubra)
LANDFIRE EVT has **no red alder type for the lower 48**. Alder only appears
in Alaska EVT classes. Red alder is a component species in riparian forests
but never dominates enough for its own EVT class.

Potential alternatives:
- **FIA plot data**: Filter `comprehensive_forest_data_pnw.geojson` for SPCD 351
  (red alder) — would be points, not polygons
- **USFS Forest Inventory**: Species-level rasters if they become available

### Bigleaf Maple (Acer macrophyllum)
Same issue — not a dominant species in LANDFIRE EVT classification.
Only 3 maple EVT types exist, none relevant to Pacific Northwest.
FIA plot data (SPCD 312) would be the alternative.

### Aspen (Populus tremuloides)
Available in LANDFIRE (codes 7009, 7011, 7061) but not yet extracted.
Rocky Mountain aspen groves are significant morel habitat after disturbance.
Could add as a future layer.
