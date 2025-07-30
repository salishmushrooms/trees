# Tree Cover Data Usage Guide

## Overview

This guide documents how tree canopy cover (TCC) data is organized and used throughout the Pacific Northwest Forest Habitat Mapping Project.

## Data Source

- **Dataset**: NLCD Tree Canopy Cover 2021 v2021-4
- **Resolution**: 30m native
- **Coverage**: CONUS-wide (includes all PNW states)
- **Values**: 0-100% tree canopy cover per pixel
- **File**: `data/raw/nlcd_tcc_conus_2021_v2021-4.tif`

## Caching System

### Regional Masks (Recommended)

Regional TCC masks cover the entire Pacific Northwest and can be reused across all species:

```
cache/species_masks/
├── pnw_tcc_mask_5pct_30m_2.0sqkm.geojson     # 5% threshold (recommended)
├── pnw_tcc_mask_10pct.geojson                # 10% threshold
├── pnw_tcc_mask_20pct.geojson                # 20% threshold (conservative)
└── pnw_tcc_mask_5pct_30m_0.5sqkm.geojson     # 5% with smaller min area
```

**Naming Convention**: `pnw_tcc_mask_{threshold}pct_{resolution}m_{min_area}sqkm.geojson`

### Creating Regional Masks

```bash
# Standard 5% mask (recommended for most species)
python scripts/create_regional_tcc_mask.py --threshold 5 --resolution 30 --min-area 2.0

# List existing masks
python scripts/create_regional_tcc_mask.py --list-existing

# Fill gaps in partial mask
python scripts/create_regional_tcc_mask.py --threshold 5 --fill-gaps
```

## Usage in Species Distribution Mapping

### Default Behavior (Regional Mask)

```bash
# Uses cached regional mask automatically
python scripts/visualization/create_species_distribution_maps_cached.py \
    --species douglas-fir \
    --tcc-threshold 5
```

### Force Species-Specific Mask

```bash
# Creates new mask for this species only (slower)
python scripts/visualization/create_species_distribution_maps_cached.py \
    --species douglas-fir \
    --use-species-tcc-mask
```

## Tree Cover Mask Processing

### Key Functions

1. **`create_tcc_mask()`** in `create_species_distribution_maps_cached.py`:
   - Creates masks of areas ABOVE the tree cover threshold
   - Uses regional mask if available
   - Falls back to species-specific processing

2. **`apply_tcc_constraints()`**:
   - Intersects species distribution with tree cover mask
   - Removes areas below threshold
   - Logs area changes for tracking

### Processing Steps

1. **Load TCC Raster**: Opens NLCD raster at native 30m resolution
2. **Apply Threshold**: Creates binary mask (1 = above threshold, 0 = below)
3. **Vectorize**: Converts raster mask to vector polygons
4. **Filter by Size**: Removes polygons smaller than minimum area
5. **Simplify**: Optional morphological smoothing for cleaner boundaries
6. **Cache**: Saves as GeoJSON for reuse

## Performance Optimization

### Caching Benefits

- **Initial Creation**: 30-60 minutes for regional mask
- **Subsequent Use**: 2-5 minutes per species (90% faster)
- **File Sizes**: 50-200MB depending on threshold and minimum area

### Memory Management

- Processes in chunks to avoid memory issues
- Automatic cleanup of temporary files
- Progress tracking for long operations

## Mapbox Integration

Tree cover masks are also used for web mapping constraints:

```
outputs/mapbox_masks/
├── pnw_tree_cover_30m_full.tif      # Full resolution
├── pnw_tree_cover_60m_half.tif      # Half resolution
├── pnw_tree_cover_90m_third.tif     # Third resolution
└── pnw_tree_cover_120m_quarter.tif  # Quarter resolution
```

## Common Use Cases

### 1. Coastal Forest Mapping
```python
# Low threshold for capturing shore pine on dunes
tcc_threshold = 5  # 5% captures sparse coastal forests
```

### 2. Dense Forest Species
```python
# Higher threshold for old-growth associated species
tcc_threshold = 20  # 20% for dense canopy species
```

### 3. Edge/Transitional Species
```python
# Medium threshold for oak woodlands, etc.
tcc_threshold = 10  # 10% balances open/closed canopy
```

## Troubleshooting

### "Regional TCC mask not found"
Create the mask first:
```bash
python scripts/create_regional_tcc_mask.py --threshold 5
```

### Memory Issues
- Increase `--min-area` to reduce polygon count
- Use coarser `--resolution` (60m instead of 30m)
- Process with `--no-morphological` to skip smoothing

### Incomplete Coverage
```bash
# Fill gaps in existing mask
python scripts/create_regional_tcc_mask.py --threshold 5 --fill-gaps
```

## Best Practices

1. **Use Regional Masks**: Create once, use for all species
2. **Standard Thresholds**: 5% (sparse), 10% (mixed), 20% (dense)
3. **Consistent Parameters**: Same resolution and min area across analyses
4. **Document Choices**: Note threshold used in species metadata
5. **Visual Validation**: Check results in QGIS before finalizing

## Related Documentation

- `archive/REGIONAL_TCC_MASK_README.md`: Detailed mask creation guide
- `scripts/create_regional_tcc_mask.py`: Mask creation script
- `scripts/visualization/create_species_distribution_maps_cached.py`: Main usage
- `CLAUDE.md`: Project-wide caching strategy