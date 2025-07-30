# QGIS Integration Guide

This guide documents the manual habitat editing workflow using QGIS for the Pacific Northwest Forest Habitat Mapping Project.

## Overview

The project uses a semi-automated workflow that combines programmatic species distribution mapping with manual refinement in QGIS. This approach leverages expert knowledge to improve habitat boundaries based on plot-level forest inventory data.

## Workflow Steps

### 1. Initial Data Generation

Generate species-specific EDIT_ME files using the distribution mapping script:

```bash
cd scripts/visualization
python create_species_distribution_maps_cached.py --species douglas-fir --workflow generate
```

This creates several files in `outputs/geojson/`:
- `{species}_merged_habitat_{buffer}_EDIT_ME.geojson` - Main habitat polygons for editing
- `{species}_individual_buffers_{buffer}.geojson` - Individual plot buffers for reference
- `{species}_plot_points.geojson` - Original plot locations
- `outputs/{species}/QGIS_EDITING_INSTRUCTIONS_{species}.md` - Species-specific editing guide

### 2. QGIS Setup

#### Required Layers
1. **Base Layer**: `outputs/plot_carbon_percentiles_latest_surveys.geojson`
   - Contains all FIA plot data with carbon percentiles
   - Use for understanding forest composition and species presence

2. **Species Habitat Layer**: `{species}_merged_habitat_{buffer}_EDIT_ME.geojson`
   - The primary editing layer
   - Contains preliminary habitat boundaries

3. **Reference Layers** (optional):
   - Individual buffer circles
   - Plot points
   - Elevation/terrain basemaps

#### Styling Recommendations

For the carbon percentile layer:
1. Style by carbon percentile value
2. Size points proportionally to percentile (0-100)
3. Color by dominant species or forest type
4. Enable labels showing species codes

### 3. Manual Editing Process

#### Key Editing Principles

1. **Expand Boundaries** where:
   - High-density clusters of the target species exist
   - Suitable elevation ranges are present
   - Continuous forest cover connects populations

2. **Contract Boundaries** where:
   - Species is absent or rare
   - Unsuitable habitat conditions exist
   - Natural barriers (ridges, valleys) create isolation

3. **Create New Polygons** for:
   - Isolated populations
   - Disjunct suitable habitat areas

#### QGIS Editing Tools

- **Node Tool**: Fine-tune polygon boundaries
- **Add Feature**: Create new habitat polygons
- **Split Features**: Divide large polygons at natural boundaries
- **Merge Features**: Combine adjacent suitable areas

### 4. Saving Edited Files

**Critical**: Save edited files with the correct naming convention:
```
{species}_merged_habitat_{buffer}_edited.geojson
```

The `_edited` suffix triggers the staged workflow processing.

### 5. Processing Edited Files

Process individual species:
```bash
python create_species_distribution_maps_cached.py --species black-cottonwood --workflow staged
```

Or bulk process all edited files:
```bash
cd scripts
python bulk_process_staged_species.py
```

## Data Interpretation

### Plot Carbon Percentiles

The `plot_carbon_percentiles_latest_surveys.geojson` file contains:
- **PLOT_ID**: Unique FIA plot identifier
- **LAT/LON**: Plot coordinates
- **Species columns**: Carbon percentile (0-100) for each species
- **INVYR**: Inventory year
- **ELEV**: Elevation in feet

Higher percentile values indicate greater relative abundance of that species at the plot.

### Species-Specific Considerations

Each species has unique habitat requirements documented in the auto-generated QGIS editing instructions:
- Elevation ranges
- Moisture preferences  
- Associated species
- Geographic limitations

## Quality Control

### Before Saving

1. Verify habitat polygons align with known species ranges
2. Check for gaps in continuous habitat
3. Ensure boundaries follow natural features
4. Compare with reference data and field knowledge

### Common Issues

- **Over-generalization**: Habitat polygons too broad
- **Under-representation**: Missing known populations
- **Edge effects**: Artificial boundaries at state/region edges
- **Elevation conflicts**: Habitat extending beyond species limits

## Integration with Mapbox

Edited and processed files are converted to Mapbox tiles for web visualization:

```bash
cd scripts/visualization
./create_mapbox_tilesets.sh
```

The resulting `.mbtiles` file contains all species layers with:
- Optimized geometry for web rendering
- Proper zoom level constraints
- Attribute preservation for popups

## Tips and Best Practices

1. **Work at Multiple Scales**: Zoom in for detail, zoom out for context
2. **Use Transparency**: Set habitat polygons to 50% opacity to see underlying plots
3. **Reference Multiple Sources**: Cross-check with elevation data, forest type maps
4. **Document Decisions**: Note rationale for major boundary changes
5. **Iterative Refinement**: Multiple editing passes improve accuracy

## Troubleshooting

### File Not Processing
- Verify exact filename format with `_edited.geojson` suffix
- Check file is saved in correct directory
- Ensure valid GeoJSON format

### Memory Issues with Large Species
- Edit in sections
- Simplify geometry before saving
- Use QGIS simplification tools if needed

### Projection Problems
- Project uses WGS84 (EPSG:4326) for all GeoJSON files
- QGIS should auto-detect, but verify CRS if issues arise