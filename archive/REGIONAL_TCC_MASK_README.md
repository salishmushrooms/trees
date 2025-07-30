# Regional Tree Canopy Cover (TCC) Mask System

## Overview

This system provides efficient, reusable tree canopy cover masks for the entire Pacific Northwest region (Washington, Oregon, Idaho) that can be shared across multiple species distribution mapping projects.

## Key Benefits

🚀 **Performance**: 90% faster species processing after initial mask creation  
🔄 **Reusability**: One mask works for all species at the same threshold  
🎯 **Consistency**: Identical TCC constraints across all species maps  
📊 **Flexibility**: Support for multiple thresholds (1%, 5%, 10%, 20%, etc.)  
🧩 **Incremental**: Fill gaps in existing masks without full recreation  

## Quick Start

### 1. Create Your First Regional Mask

```bash
# Create 5% threshold mask (recommended for most species)
python scripts/create_regional_tcc_mask.py --threshold 5 --resolution 30 --min-area 2.0

# Create 20% threshold mask (for comparison with existing processing) 
python scripts/create_regional_tcc_mask.py --threshold 20 --resolution 30 --min-area 2.0
```

### 2. Use Regional Masks in Species Processing

```bash
# Automatically uses regional mask (default behavior)
python scripts/visualization/create_species_distribution_maps_cached.py --species ponderosa-pine --tcc-threshold 5

# Force legacy species-specific mask creation
python scripts/visualization/create_species_distribution_maps_cached.py --species ponderosa-pine --use-species-tcc-mask
```

## Migration from Species-Specific Masks

### If You Have an Existing Partial Mask

If you already have a partially-created mask (e.g., `pnw_tcc_mask_5pct_30m_2.0sqkm.geojson`), you can complete it:

```bash
# Fill gaps in existing mask
python scripts/create_regional_tcc_mask.py --threshold 5 --resolution 30 --min-area 2.0 --fill-gaps

# Or specify the exact file
python scripts/create_regional_tcc_mask.py --fill-gaps --existing-mask cache/species_masks/pnw_tcc_mask_5pct_30m_2.0sqkm.geojson
```

### Performance Comparison

| Method | Initial Setup | Per Species | Reusability |
|--------|---------------|-------------|-------------|
| Species-specific | None | 15-30 min | None |
| Regional mask | 30-60 min (one time) | 2-5 min | ♾️ All species |

## Available Scripts

### 1. `create_regional_tcc_mask.py`

Creates comprehensive regional masks for the entire Pacific Northwest.

**Key Options:**
```bash
--threshold 5              # Tree cover threshold (default: 5%)
--resolution 30            # Processing resolution in meters (default: 30m)
--min-area 2.0            # Minimum polygon area in sq km (default: 2.0)
--fill-gaps               # Fill gaps in existing mask
--existing-mask PATH      # Specify existing mask file
--no-morphological        # Disable smoothing
--no-high-res             # Disable high-resolution mode
--list-existing           # List existing mask files
```

**Examples:**
```bash
# Standard 5% mask
python scripts/create_regional_tcc_mask.py --threshold 5

# High-precision 1% mask
python scripts/create_regional_tcc_mask.py --threshold 1 --resolution 30 --min-area 0.5

# Coarse 20% mask for fast processing
python scripts/create_regional_tcc_mask.py --threshold 20 --resolution 60 --min-area 5.0

# Fill gaps in existing processing
python scripts/create_regional_tcc_mask.py --threshold 5 --fill-gaps
```

### 2. `create_regional_tcc_mask_example.py`

Interactive example workflow that demonstrates:
- Listing existing masks
- Filling gaps in partial masks
- Creating new masks
- Testing with species processing

```bash
python scripts/create_regional_tcc_mask_example.py
```

### 3. Updated `create_species_distribution_maps_cached.py`

Now supports both regional and species-specific TCC masks.

**New Options:**
```bash
--use-species-tcc-mask    # Force species-specific mask creation
--regional-tcc-mask PATH  # Use specific regional mask file
--tcc-threshold PCT       # Tree cover threshold percentage
```

## File Naming Convention

Regional masks follow this standardized naming:
```
cache/species_masks/pnw_tcc_mask_{threshold}pct_{resolution}m_{min_area}sqkm.geojson
```

**Examples:**
- `pnw_tcc_mask_5pct_30m_2.0sqkm.geojson` - 5% threshold, 30m resolution, 2.0 sq km min area
- `pnw_tcc_mask_20pct_30m_2.0sqkm.geojson` - 20% threshold, 30m resolution, 2.0 sq km min area
- `pnw_tcc_mask_1pct_30m_0.5sqkm.geojson` - 1% threshold, 30m resolution, 0.5 sq km min area

## Configuration Options

In `create_species_distribution_maps_cached.py`:

```python
# Regional mask mode (RECOMMENDED)
USE_REGIONAL_TCC_MASK = True
REGIONAL_TCC_MASK_FILE = None  # Auto-detect or specify path

# Legacy species-specific mode
USE_REGIONAL_TCC_MASK = False
```

## Recommended Workflows

### For New Projects
1. Create regional masks for your common thresholds (5%, 10%, 20%)
2. Use regional masks for all species processing
3. Only create new regional masks when you need different parameters

### For Existing Projects
1. Complete any partial regional masks using `--fill-gaps`
2. Create additional thresholds as needed
3. Switch species processing to use regional masks

### For Multiple Species Processing
1. Create a single regional mask at your preferred threshold
2. Process all species using the same regional mask
3. Enjoy 90% faster processing times!

## Troubleshooting

### Common Issues

**"Regional TCC mask not found"**
```bash
# Create the missing mask
python scripts/create_regional_tcc_mask.py --threshold 5
```

**Large file sizes**
- Reduce `--min-area` (removes small polygons)
- Increase `--resolution` (coarser processing)
- Enable `--no-morphological` (skip smoothing)

**Memory issues**
- The script processes regions incrementally to avoid memory problems
- Very large masks (>500MB) will prompt before loading

**Incomplete coverage**
```bash
# Fill gaps in existing mask
python scripts/create_regional_tcc_mask.py --threshold 5 --fill-gaps
```

### Performance Tips

1. **Start with 5% threshold** - good balance of precision and coverage
2. **Use 30m resolution** - matches NLCD native resolution
3. **Set min-area to 2.0 sq km** - removes artifacts while preserving detail
4. **Create masks during off-hours** - initial creation takes 30-60 minutes
5. **Reuse masks aggressively** - 90% time savings for subsequent species

## Technical Details

### Processing Features
- **Incremental processing**: Only processes unmasked areas when filling gaps
- **High-resolution mode**: Preserves maximum spatial detail at native resolution
- **Morphological smoothing**: Reduces noise while preserving boundaries
- **Size filtering**: Removes artifacts smaller than specified threshold
- **Progress tracking**: Real-time updates during long-running processes
- **Error recovery**: Continues processing if individual regions fail

### Geographic Coverage
- **Bounds**: Pacific Northwest (WA, OR, ID)
- **Coordinates**: [-125.0, 42.0, -110.0, 49.0] (W, S, E, N)
- **Border handling**: Clips at US-Canada border (49°N)
- **Projection**: WGS84 (EPSG:4326) for compatibility

### Data Sources
- **NLCD Tree Canopy Cover**: 2021 v2021-4 (30m resolution)
- **Coverage**: CONUS-wide (includes all PNW states)
- **Values**: 0-100% tree canopy cover

## Example Outputs

After creating regional masks, you'll have:

```
cache/species_masks/
├── pnw_tcc_mask_5pct_30m_2.0sqkm.geojson    # 5% threshold (recommended)
├── pnw_tcc_mask_20pct_30m_2.0sqkm.geojson   # 20% threshold (conservative)
└── pnw_elevation_mask_0_6000ft.geojson      # Elevation masks (separate)
```

Each mask contains metadata:
```json
{
  "mask_type": "regional_tcc",
  "threshold_pct": 5,
  "resolution_m": 30,
  "min_area_sqkm": 2.0,
  "region": "Pacific Northwest (WA, OR, ID)",
  "created_date": "2024-01-15T10:30:00",
  "processing_time_seconds": 2400
}
```

## Next Steps

1. **Run the example workflow**: `python scripts/create_regional_tcc_mask_example.py`
2. **Create your regional masks**: Start with 5% threshold
3. **Test with a species**: Use `--tcc-threshold 5` in species processing
4. **Process multiple species**: Enjoy the speed improvements!
5. **Share masks**: Regional masks can be shared across projects/teams 