# Bioregion Creation Guide

## Overview

This document provides a comprehensive guide for creating bioregions based on our successful coastal forest bioregion implementation. The process combines raster data analysis with species-specific ecological constraints to create natural, irregular boundaries.

## Key Process Components

### 1. Data Sources Required

#### Essential Raster Data
- **Tree Cover Data**: `outputs/mapbox_masks/pnw_tree_cover_30m_full.tif` (30m resolution)
- **Elevation Data**: `outputs/mapbox_masks/pnw_elevation_120m_mapbox.tif` (120m resolution)
- **Base Boundary**: User-defined shapefile for general area of interest

#### Species Validation Data
- **Plot Data**: `outputs/plot_carbon_percentiles_latest_surveys.geojson`
- **FIA Database**: For species occurrence validation

### 2. Core Script Architecture

#### Primary Script
**File**: `scripts/create_[bioregion_name]_shapefile_constrained.py`

#### Key Functions
1. **Boundary Loading**: Load user-defined base boundary shapefile
2. **Raster Constraints**: Apply elevation and tree cover filters within boundary
3. **Polygon Extraction**: Convert raster mask to vector polygons
4. **Artifact Removal**: Clean up raster-to-vector conversion artifacts
5. **Species Validation**: Validate against known species occurrence data

### 3. Configuration Parameters

#### Typical Bioregion Parameters
```python
# Bioregion-specific constraints
MAX_ELEVATION_FT = 1000        # Species-dependent elevation limit
MIN_TREE_COVER_PCT = 10        # Forest definition threshold
OUTPUT_RESOLUTION = 120        # Match elevation data resolution

# Artifact removal settings
min_area_sqkm = 0.08          # Minimum polygon area (80,000 m²)
inward_buffer = -0.006        # 600m buffer to avoid edge artifacts
```

## Implementation Workflow

### Phase 1: Base Boundary Creation
1. Create initial boundary shapefile using QGIS or other GIS tools
2. Define broad area of interest based on ecological knowledge
3. Save as `outputs/bioregions/[bioregion_name]_broad.shp`

### Phase 2: Raster Constraint Application
1. **Load base boundary** and convert to processing CRS (Albers Equal Area)
2. **Apply inward buffer** to avoid edge artifacts (typically 600m)
3. **Load elevation data** and apply species-specific elevation constraints
4. **Load tree cover data** and apply forest definition threshold
5. **Combine constraints** using logical AND operations

### Phase 3: Polygon Processing
1. **Convert raster to polygons** using `rasterio.features.shapes()`
2. **Remove edge artifacts** through multiple filtering steps:
   - Single-pixel removal (< 1.5 × pixel area)
   - Rectangular strip removal (high aspect ratio + thin width)
   - Minimum area filtering (typically 0.08 km²)
3. **Merge polygons** into coherent bioregion geometry

### Phase 4: Validation and Output
1. **Species validation** against known occurrence data
2. **Generate summary statistics** and metadata
3. **Save final outputs** in multiple formats

## Common Issues and Solutions

### 1. Edge Artifacts
**Problem**: Thin rectangular strips appear along boundary edges
**Symptoms**: 
- Single-pixel wide rectangles (~72-120m wide)
- Perfect rectangular shapes following raster grid
- Areas as small as 7,000-35,000 m²

**Solutions Applied**:
- Inward buffer on base boundary (600m typical)
- Single-pixel removal filter
- Rectangular artifact detection and removal
- Increased minimum area threshold (0.08 km²)

### 2. Raster Processing Edge Effects
**Problem**: Boundary pixels get assigned values due to resampling
**Symptoms**: Valid elevations/tree cover along boundaries that should be excluded

**Solutions**:
```python
# Apply aggressive inward buffer
boundary_geom = boundary_geom.buffer(-0.006)  # 600m inward

# Multi-stage artifact removal
isolated_mask = gdf_geo['area_sqkm'] < single_pixel_threshold
is_thin_strip = (aspect_ratio > 10) & (min_dimension < 150)
```

### 3. Tree Cover NoData Handling
**Problem**: Water areas may have 0% tree cover instead of NoData
**Symptoms**: Offshore areas included in bioregion

**Solutions**:
```python
# Ensure proper NoData handling during resampling
tcc_data = np.full_like(elev_data, np.nan, dtype=np.float32)
reproject(..., dst_nodata=np.nan)
```

## File Organization

### Script Location
```
scripts/create_[bioregion_name]_shapefile_constrained.py
```

### Input Files
```
outputs/bioregions/[bioregion_name]_broad.shp          # Base boundary
outputs/mapbox_masks/pnw_elevation_120m_mapbox.tif     # Elevation data
outputs/mapbox_masks/pnw_tree_cover_30m_full.tif       # Tree cover data
outputs/plot_carbon_percentiles_latest_surveys.geojson # Validation data
```

### Output Files
```
outputs/bioregions/[bioregion_name]_constrained.geojson    # Final bioregion
outputs/bioregions/[bioregion_name]_summary.json          # Metadata
outputs/bioregions/[species]_validation_plots.geojson     # Validation points
```

## Species-Specific Parameters

### Coastal Forest Example
```python
MAX_ELEVATION_FT = 1000        # Coastal-influenced species limit
MIN_TREE_COVER_PCT = 10        # Coastal forest definition
TARGET_SPECIES = "shore pine"  # Pinus contorta var. contorta
VALIDATION_SPCD = 108          # Lodgepole/shore pine species code
```

### Future Bioregion Templates
Consider these parameters for other bioregions:

#### Alpine/Subalpine Forest
```python
MAX_ELEVATION_FT = 7000        # High elevation limit
MIN_ELEVATION_FT = 3000        # Exclude lowland forests
MIN_TREE_COVER_PCT = 15        # Denser forest requirement
```

#### Riparian Forest
```python
MAX_DISTANCE_TO_WATER_M = 500  # Proximity to water bodies
MIN_TREE_COVER_PCT = 20        # Riparian forest definition
SLOPE_MAX_DEGREES = 15         # Relatively flat areas
```

## Validation Approach

### Species Validation
1. **Load validation plots** for target species
2. **Filter by bioregion constraints** (elevation, etc.)
3. **Calculate capture rate** (% of plots within bioregion)
4. **Target capture rate**: >70% for good validation

### Quality Metrics
- **Area calculation**: Final bioregion area in km²
- **Capture rate**: % of target species plots included
- **Polygon count**: Number of separate forest patches
- **Artifact count**: Number of removed edge artifacts

## Debugging and Troubleshooting

### Debug Output Recommendations
```python
# Elevation range validation
valid_elevations = elev_data_ft[elevation_mask]
print(f"Elevation range: {np.min(valid_elevations):.1f}ft to {np.max(valid_elevations):.1f}ft")

# Tree cover statistics
valid_tcc = tcc_data[tree_cover_mask]
print(f"Tree cover range: {np.min(valid_tcc):.1f}% to {np.max(valid_tcc):.1f}%")

# Artifact removal tracking
print(f"Removed {isolated_count} single-pixel polygons")
print(f"Removed {strip_count} rectangular edge artifacts")
```

### Common Debug Checks
1. **Boundary loading**: Verify CRS and geometry validity
2. **Raster alignment**: Ensure elevation and tree cover align properly
3. **Constraint application**: Check that masks have reasonable pixel counts
4. **Polygon extraction**: Verify reasonable polygon count and sizes

## Best Practices

### 1. Incremental Development
- Start with basic tree cover + elevation constraints only
- Add complexity (water exclusion, shape filters) incrementally
- Test each addition separately to isolate issues

### 2. Parameter Tuning
- Start with permissive parameters and tighten as needed
- Use visual inspection in QGIS to validate results
- Adjust minimum area thresholds based on observed artifacts

### 3. Documentation
- Record all parameter changes and their effects
- Document species-specific ecological rationale
- Maintain consistent file naming conventions

### 4. Validation
- Always validate against known species occurrence data
- Use multiple validation metrics (area, capture rate, visual inspection)
- Compare results with existing ecological boundaries where available

## Multi-Bioregion Management

### Combining Multiple Bioregions

Once you have created multiple individual bioregions, use the automated combination script to merge them into a single comprehensive dataset:

**Script**: `scripts/combine_bioregions.py`

#### Usage
```bash
# Automatically combines all *_constrained.geojson files
python scripts/combine_bioregions.py
```

#### Features
- **Automatic discovery**: Finds all `*_constrained.geojson` files in `outputs/bioregions/`
- **Standardized metadata**: Adds consistent region identifiers and timestamps
- **Area calculations**: Computes total area statistics for each region
- **Validation**: Checks geometry validity and required fields
- **Summary reporting**: Creates detailed JSON summary with statistics

#### Outputs
```
outputs/bioregions/all_bioregions_combined.geojson    # Combined dataset
outputs/bioregions/bioregions_combination_summary.json    # Statistics & metadata
```

#### Standardized Fields Added
- `region_name`: Full region name (e.g., "coastal_forest")
- `region_code`: 8-character code for joining (e.g., "COASTAL")
- `area_km2`: Calculated area in square kilometers
- `bioregion_type`: Always "forest_bioregion"
- `dataset`: "PNW_Forest_Bioregions"
- `created_date`: ISO timestamp of combination

#### Web Mapping Integration
The combined file is optimized for web mapping platforms:
- Use `region_name` field for layer styling and filtering
- Use `region_code` for data joins and analysis
- All geometries in WGS84 (EPSG:4326) for web compatibility

### Future Bioregion Workflow

When creating new bioregions:

1. **Create individual bioregion** using the standard workflow
2. **Save as** `outputs/bioregions/[bioregion_name]_constrained.geojson`
3. **Run combination script** to update the master dataset
4. **Validate** the combined output using the built-in validation

The script automatically handles:
- CRS standardization
- Metadata consistency
- Area calculations
- Summary statistics
- Error handling

## Future Enhancements

### Potential Improvements
1. **Morphological operations** for smoother boundaries
2. **Multi-criteria decision analysis** for complex constraints
3. **Automated parameter optimization** based on validation metrics
4. **Climate data integration** for additional ecological constraints

### Scalability Considerations
- **Memory management** for larger rasters
- **Parallel processing** for multiple bioregions
- **Database optimization** for species validation queries
- **Automated combination** triggered by new bioregion creation

---

*This guide is based on the successful implementation of the coastal forest bioregion and should be updated as new bioregions are created and lessons learned.*