# Black Cottonwood Visual Layer — Summary & Future Work

*Created: March 2026*

## Accomplishments

### 1. Bug Fix: SPCD 746 → 747
- `extract_comprehensive_forest_data.py` was using SPCD 746 (quaking aspen) labeled as "Black Cottonwood"
- Correct code is **SPCD 747** (Populus balsamifera ssp. trichocarpa)
- `create_species_distribution_maps_cached.py` already had the correct code (747)
- Fixed in extraction script with explanatory comment

### 2. Mapbox → MapLibre Migration (Comments)
- Updated references in `CLAUDE.md`, `create_species_distribution_maps_cached.py`
- Current stack: **MapLibre + PMTiles** (served via Cloudflare R2)
- Legacy Mapbox/MBTiles references may remain in older scripts

### 3. New Cottonwood Visual Layer (LANDFIRE EVT)
Created a high-quality cottonwood distribution layer using LANDFIRE 2023 Existing Vegetation Type data.

**Output files:**
- `outputs/geojson/black_cottonwood_landfire_evt.geojson` (3.7 MB) — 2,534 polygons, 8,749 km²
- `outputs/black_cottonwood_landfire.pmtiles` (2.4 MB) — z4-z12, layer name "cottonwood"

**Properties:** `species`, `species_code` (747), `cottonwood_type` (primary/secondary), `area_km2`, `source`

**Breakdown:**
- Primary (SAF 222 — Black Cottonwood-Willow): 1,362 polygons, 4,664 km²
- Secondary (SAF 235 — Cottonwood-Willow): 1,172 polygons, 4,085 km²

## Data Sources Evaluated

### LANDFIRE EVT 2023 (Selected)
- **What**: 30m CONUS raster classifying actual vegetation types
- **Why it works**: Shows where cottonwood *actually grows*, not just where it could grow
- **Source file**: `data/raw/LF2023_EVT_CONUS/Tif/LF2023_EVT_CONUS.tif` (2.8 GB, EPSG:5070)
- **Lookup table**: `data/raw/LF2024_EVT.csv` (1,070 vegetation types)
- **FIA validation**: 2% false positive rate (high specificity). Only 7% of FIA cottonwood plots fall within polygons, which is expected — cottonwood grows in narrow riparian corridors while FIA plots are randomly distributed and coordinates are fuzzed ±0.5-1 mile.

### USGS Woody Riparian Species Habitat Suitability (Rejected)
- **What**: 100m raster, ensemble score 0-10
- **File**: `data/raw/usgs_riparian/X10th_combined.tif` (63 MB, ESRI:102008)
- **Why rejected**: Too broad for visual layer. Score==10 covered 86K km² in PNW. Even with elevation filter (<4000 ft), 131K km². Shows climate/terrain suitability, not actual presence.
- **Validation**: 85% of FIA cottonwood plots scored >=9 (good as a suitability model, just not for visualization)

## LANDFIRE EVT Codes Used

### Primary — SAF 222: Black Cottonwood-Willow
| Code | Name | PNW Area |
|------|------|----------|
| 7156 | North Pacific Lowland Riparian Forest | 2,853 km² |
| 7158 | North Pacific Montane Riparian Woodland | 824 km² |
| 7157 | North Pacific Lowland Riparian Forest-Hardwood | smaller |
| 7159 | North Pacific Montane Riparian Woodland-Hardwood | smaller |
| 9165 | North Pacific Shrub Swamp | smaller |

### Secondary — SAF 235: Cottonwood-Willow
| Code | Name | PNW Area |
|------|------|----------|
| 9012 | N. Rocky Mtn Lower Montane Riparian Woodland | 2,234 km² |
| 9512 | N. Rocky Mtn Lower Montane Riparian Woodland (variant) | smaller |
| 9019, 9519 | N. Rocky Mtn Lower Montane Riparian Shrubland | smaller |
| 9002, 9502, 9003, 9503 | Additional Rocky Mtn riparian types | smaller |

## Processing Pipeline

```
LANDFIRE TIF (EPSG:5070, 30m)
  → Reproject to WGS84 (~120m)
  → Minimal morphological smoothing (2×2 kernel, 1 iteration)
  → Polygonize
  → Filter ≥1 km²
  → Simplify (0.003° tolerance)
  → GeoJSON
  → Tippecanoe → PMTiles (z4-z12)
```

## Large File Cleanup

These files can be deleted to recover ~3.7 GB:

| File | Size | Notes |
|------|------|-------|
| `data/raw/LF2023_EVT_CONUS.zip` | 3.6 GB | TIF already extracted |
| `data/raw/usgs_riparian/X10th_combined.zip` | 59 MB | TIF already extracted |
| `data/raw/usgs_riparian/cottonwood_score10_pnw_250m.tif` | 569 KB | Intermediate, superseded |
| `outputs/geojson/black_cottonwood_usgs_riparian.geojson` | 1.1 MB | Superseded by LANDFIRE version |

**Keep** (useful for future work):
- `data/raw/LF2023_EVT_CONUS/Tif/LF2023_EVT_CONUS.tif` (2.8 GB) — full CONUS EVT raster
- `data/raw/LF2024_EVT.csv` — lookup table for all 1,070 EVT codes
- `data/raw/usgs_riparian/X10th_combined.tif` (63 MB) — valid suitability model, may be useful as a filter

## Available LANDFIRE Vegetation Types for Future Work

The LANDFIRE EVT dataset covers 62 PNW-relevant vegetation types. Key species with multiple EVT types:

| Species | EVT Types | Notes |
|---------|-----------|-------|
| Oak (Quercus) | 154 types | Includes Oregon white oak woodlands/savannas |
| Hemlock | 16 types | Western hemlock, mountain hemlock forests |
| Cedar | 15 types | Western red cedar forests |
| Maple | 14 types | Bigleaf maple mixed forests |
| Douglas-fir | 8 types | Most abundant PNW conifer |
| Aspen | 7 types | Montane/boreal deciduous forests |
| Alder | 5 types | Red alder riparian/early succession |

This same LANDFIRE extraction pipeline could produce visual layers for any of these species using the same approach: identify relevant EVT codes → extract from CONUS TIF → polygonize → simplify → PMTiles.

## Lessons Learned

1. **Suitability ≠ Presence**: Climate/terrain suitability models (USGS) are far too broad for visual layers. LANDFIRE EVT maps actual vegetation, which is what we need.
2. **FIA undersamples riparian**: Random plot placement misses narrow corridors. Only 7% overlap with known cottonwood areas is expected, not a failure.
3. **Morphological smoothing matters**: Too aggressive (iterations=2, 3×3) expanded area 3.7× and tripled file size. Minimal smoothing (iterations=1, 2×2) preserved accuracy.
4. **LANDFIRE is a gold mine**: 1,070 vegetation types at 30m. One 2.8 GB download enables layers for dozens of species with the same pipeline.
5. **Tippecanoe PMTiles workflow**: `tippecanoe -o out.pmtiles -z 12 -Z 4 --drop-densest-as-needed --simplification=10 -l layername` produces compact, web-ready tiles.

## Next Steps

- [ ] Upload PMTiles to R2 (`tiles.salishmushrooms.com`)
- [ ] Clean up large files listed above
- [ ] Consider LANDFIRE layers for other tree species (oak, alder, hemlock, etc.)
- [ ] Integrate with morel-maps app if cottonwood/riparian context is relevant to foraging
