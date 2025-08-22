# Bioregion Creation Guide

## Overview

This document provides a comprehensive guide for creating bioregions based on our successful coastal forest bioregion implementation. The process combines raster data analysis with species-specific ecological constraints to create natural, irregular boundaries.

## Key Process Components

### 1. Data Sources Required

#### Essential Raster Data
- **Tree Cover Data**: `outputs/mapbox_masks/pnw_tree_cover_30m_full.tif` (30m resolution)
- **Elevation Data**: `outputs/mapbox_masks/pnw_elevation_120m_mapbox.tif` (120m resolution)
- **Base Boundary**: User-defined shapefile for general area of interest

#### Climate Data (Optional but Recommended)
- **Precipitation Data**: `outputs/climate/pnw_precipitation_annual_normals.tif` (PRISM 30-year normals, 120m resolution)
- **Precipitation Zones**: `outputs/climate/pnw_precipitation_zones.geojson` (Pre-classified precipitation categories)

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

# Climate constraints (when using precipitation data)
MIN_PRECIPITATION_IN = 80      # Minimum annual precipitation (inches)
MAX_PRECIPITATION_IN = 150     # Maximum annual precipitation (inches)

# Artifact removal settings
min_area_sqkm = 0.08          # Minimum polygon area (80,000 m²)
inward_buffer = -0.006        # 600m buffer to avoid edge artifacts
```

## Implementation Workflow

### Raster-Based Approach (Coastal, Olympic Mountains)

#### Phase 1: Base Boundary Creation
1. Create initial boundary shapefile using QGIS or other GIS tools
2. Define broad area of interest based on ecological knowledge
3. Save as `outputs/bioregions/[bioregion_name]_broad.shp`

#### Phase 2: Raster Constraint Application
1. **Load base boundary** and convert to processing CRS (Albers Equal Area)
2. **Apply inward buffer** to avoid edge artifacts (typically 600m)
3. **Load elevation data** and apply species-specific elevation constraints
4. **Load tree cover data** and apply forest definition threshold
5. **(Optional) Load precipitation data** and apply climate constraints
   - Use for coastal rainforest bioregions (>150 inches)
   - Distinguish wet/dry Cascades boundaries (50-80 inch transition)
   - Identify rain shadow effects (<20 inches)
6. **Combine constraints** using logical AND operations

#### Phase 3: Polygon Processing
1. **Convert raster to polygons** using `rasterio.features.shapes()`
2. **Remove edge artifacts** through multiple filtering steps:
   - Single-pixel removal (< 1.5 × pixel area)
   - Rectangular strip removal (high aspect ratio + thin width)
   - Minimum area filtering (typically 0.08 km²)
3. **Merge polygons** into coherent bioregion geometry

#### Phase 4: Validation and Output
1. **Species validation** against known occurrence data
2. **Generate summary statistics** and metadata
3. **Save final outputs** in multiple formats

### Species-Based Approach (Eastern & High Cascades)

#### Phase 1: Plot Classification
1. **Load FIA plot data** with species carbon information
2. **Define species criteria** for target bioregion
3. **Apply classification logic** to identify indicator plots
4. **Filter spatially** to base boundary area

#### Phase 2: Bioregion Generation
1. **Create buffers** around indicator plots (typically 5km radius)
2. **Merge overlapping buffers** to form continuous regions
3. **Clip to base boundary** to respect ecological boundaries
4. **Filter by minimum area** to remove small fragments

#### Phase 3: Tree Cover Refinement
1. **Apply tree cover constraint** to buffered regions
2. **Convert to polygons** using forest mask
3. **Remove small fragments** below minimum area threshold
4. **Merge final polygons** into coherent bioregion

#### Phase 4: Validation and Output
1. **Save indicator plots** for review and validation
2. **Calculate area statistics** and plot counts
3. **Generate summary metadata** with species criteria
4. **Export final bioregion** and supporting files

## Species-Based Classification Methodology

### Rationale for Eastern & High Cascades

Traditional elevation-based bioregion delineation fails to accurately separate Eastern Cascades dry forests from High Cascades mesic forests because:

- **Elevation overlap**: Ponderosa pine and true firs can occur at similar elevations (2000-4000ft)
- **Local variability**: Microclimates create species mixing independent of elevation
- **Rain shadow effects**: Not strictly elevation-dependent
- **Fire history**: Disturbance patterns affect species composition more than elevation

### Species Composition Approach

The species-based methodology uses FIA plot data to identify forest types based on actual tree species presence and carbon storage:

#### Eastern Cascades Classification Logic
```python
eastern_cascades_criteria = (
    # Primary dry forest indicators (OR condition)
    (ponderosa_pine > 20) | (lodgepole_pine > 30) | 
    (western_larch > 0) | (western_juniper > 10)
) & (
    # Exclusion of mesic species (AND condition)
    (western_hemlock < 50) & (pacific_silver_fir < 40) & 
    (western_redcedar < 50)
) & (
    # Forest definition
    total_carbon > 30
)

# Alternative: Species ratio approach
dry_species_ratio = dry_carbon / total_carbon > 0.4
```

#### High Cascades Classification Logic
```python
high_cascades_criteria = (
    # Primary subalpine indicators (OR condition)
    (pacific_silver_fir > 30) | (subalpine_fir > 30) |
    (mountain_hemlock > 30) | (noble_fir > 20) |
    (engelmann_spruce > 20)
) & (
    # Exclusion of dry species (AND condition)
    (ponderosa_pine < 30) & (western_juniper == 0) &
    (western_larch < 50)
) & (
    # Forest definition
    total_carbon > 40
)

# Alternative: Subalpine species ratio approach
subalpine_ratio = subalpine_carbon / total_carbon > 0.4
```

### Buffer-Based Region Creation

After plot classification, continuous bioregions are created through:

1. **Spatial clustering**: 5km buffers around indicator plots
2. **Overlap merging**: Union of overlapping buffers
3. **Boundary clipping**: Intersection with broad ecological boundaries
4. **Tree cover constraint**: Final refinement using canopy cover data

### Advantages of Species-Based Approach

- **Ecologically accurate**: Based on actual species presence
- **Climate-responsive**: Reflects moisture/temperature gradients
- **Disturbance-aware**: Captures fire and management effects
- **Validation-ready**: Directly tied to measurable forest attributes
- **Transferable**: Can be applied to other regions with FIA data

## Common Issues and Solutions

### 1. Elevation Data Units - **CRITICAL**
**Problem**: Elevation raster may already be in feet instead of meters
**File**: `outputs/mapbox_masks/pnw_elevation_120m_mapbox.tif` is in **feet**, not meters

**Symptoms**:
- Elevation constraints capturing wrong elevation zones  
- 1000-6000ft constraint only finding 300-1800ft areas
- Bioregions dominated by unexpectedly low elevations

**Solution**:
```python
# WRONG - Double converts feet to feet  
elev_data_ft = convert_elevation_to_feet(elev_data)

# CORRECT - Data is already in feet
elev_data_ft = elev_data  # Data is already in feet
```

**⚠️ Always verify elevation data units before applying conversion functions!**

### 2. Edge Artifacts
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

#### Raster-Based Scripts
```
scripts/create_[bioregion_name]_shapefile_constrained.py
```

#### Species-Based Scripts  
```
scripts/create_[bioregion_name]_species_based.py
scripts/analyze_cascades_species_composition.py  # Analysis helper
```

### Input Files

#### Raster-Based Bioregions
```
outputs/bioregions/[bioregion_name]_broad.shp          # Base boundary
outputs/mapbox_masks/pnw_elevation_120m_mapbox.tif     # Elevation data (feet)
outputs/mapbox_masks/pnw_tree_cover_30m_full.tif       # Tree cover data
outputs/plot_carbon_percentiles_latest_surveys.geojson # Validation data
```

#### Species-Based Bioregions
```
outputs/bioregions/[bioregion_name]_broad.shp          # Base boundary
outputs/plot_carbon_percentiles_latest_surveys.geojson # Species composition data
outputs/mapbox_masks/pnw_tree_cover_30m_full.tif       # Tree cover constraint
```

### Output Files

#### Raster-Based Outputs
```
outputs/bioregions/[bioregion_name]_constrained.geojson    # Final bioregion
outputs/bioregions/[bioregion_name]_summary.json          # Metadata
outputs/bioregions/[species]_validation_plots.geojson     # Validation points
```

#### Species-Based Outputs
```
outputs/bioregions/[bioregion_name]_species_based.geojson  # Final bioregion
outputs/bioregions/[bioregion_name]_indicator_plots.geojson # Classified plots
outputs/bioregions/[bioregion_name]_species_summary.json   # Species criteria & stats
```

## Species-Specific Parameters

### Coastal Forest Example
```python
MAX_ELEVATION_FT = 1000        # Coastal-influenced species limit
MIN_TREE_COVER_PCT = 10        # Coastal forest definition
MIN_PRECIPITATION_IN = 50      # Coastal moisture requirement
TARGET_SPECIES = "shore pine"  # Pinus contorta var. contorta
VALIDATION_SPCD = 108          # Lodgepole/shore pine species code
```

### Olympic Mountains
```python
MIN_ELEVATION_FT = 500         # Exclude lowland areas
MAX_ELEVATION_FT = 6000        # Exclude alpine peaks
MIN_TREE_COVER_PCT = 30        # Denser montane forest requirement
MIN_PRECIPITATION_IN = 80      # High precipitation requirement
```

### Olympic Rainforest (Precipitation-Based)
```python
MIN_ELEVATION_FT = 0           # Sea level to montane
MAX_ELEVATION_FT = 3000        # Below subalpine zone
MIN_TREE_COVER_PCT = 50        # Dense rainforest canopy
MIN_PRECIPITATION_IN = 150     # Extreme precipitation threshold
# Creates temperate rainforest bioregion
```

### Eastern Cascades (Elevation-Based - Original Approach)
```python
MAX_ELEVATION_FT = 3500        # Rain shadow effect upper limit
MIN_TREE_COVER_PCT = 20        # Drier forest requirement
MAX_PRECIPITATION_IN = 25      # Dry side precipitation limit
# Note: This approach was enhanced with species-based classification
```

### Eastern Cascades (Species-Based - Recommended)
```python
# Primary dry forest indicators
PONDEROSA_PINE_MIN = 20        # kg/ha carbon threshold
LODGEPOLE_PINE_MIN = 30        # kg/ha carbon threshold
WESTERN_LARCH_MIN = 0          # Any presence indicates Eastern Cascades
WESTERN_JUNIPER_MIN = 10       # kg/ha carbon threshold

# Exclusion criteria (mesic species limits)
WESTERN_HEMLOCK_MAX = 50       # kg/ha carbon threshold
PACIFIC_SILVER_FIR_MAX = 40    # kg/ha carbon threshold
WESTERN_REDCEDAR_MAX = 50      # kg/ha carbon threshold
```

### High Cascades (Species-Based - Recommended)
```python
# Primary subalpine/true fir indicators
PACIFIC_SILVER_FIR_MIN = 30    # kg/ha carbon threshold
SUBALPINE_FIR_MIN = 30         # kg/ha carbon threshold
MOUNTAIN_HEMLOCK_MIN = 30      # kg/ha carbon threshold
NOBLE_FIR_MIN = 20             # kg/ha carbon threshold
ENGELMANN_SPRUCE_MIN = 20      # kg/ha carbon threshold

# Exclusion criteria (dry species limits)
PONDEROSA_PINE_MAX = 30        # kg/ha carbon threshold
WESTERN_JUNIPER_MAX = 0        # No juniper in High Cascades
WESTERN_LARCH_MAX = 50         # kg/ha carbon threshold
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
- **Mapbox tiles creation**: Automatically generates mbtiles using Tippecanoe (max zoom 14)

#### Outputs
```
outputs/bioregions/all_bioregions_combined.geojson         # Combined dataset
outputs/bioregions/bioregions_combination_summary.json     # Statistics & metadata
outputs/mapbox_mbtiles/pnw_forest_bioregions.mbtiles       # Mapbox tiles (max zoom 14)
```

#### Standardized Fields Added
- `region_name`: Full region name (e.g., "coastal_forest")
- `region_code`: 8-character code for joining (e.g., "COASTAL")
- `area_km2`: Calculated area in square kilometers
- `bioregion_type`: Always "forest_bioregion"
- `dataset`: "PNW_Forest_Bioregions"
- `created_date`: ISO timestamp of combination

#### Web Mapping Integration
The script creates both GeoJSON and Mapbox tiles optimized for web mapping:

**GeoJSON Features:**
- Use `region_name` field for layer styling and filtering
- Use `region_code` for data joins and analysis
- All geometries in WGS84 (EPSG:4326) for web compatibility

**Mapbox Tiles Features:**
- Layer name: `bioregions`
- Max zoom level: 14 (suitable for regional analysis)
- Includes properties: `region_name`, `region_code`, `area_km2`, `bioregion_type`
- Ready for direct upload to Mapbox Studio

**Prerequisites:**
- Tippecanoe must be installed: `brew install tippecanoe`

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

## Using Precipitation Data for Enhanced Bioregion Definition

### Why Precipitation Matters
Precipitation is a primary driver of forest composition in the Pacific Northwest, creating distinct ecological zones:

- **Olympic Rainforest**: 150-240+ inches annually - supports massive temperate rainforest
- **Coastal Forests**: 80-150 inches - wet coastal conifer forests
- **Western Cascades**: 50-80 inches - productive Douglas-fir/hemlock forests
- **Western Valleys**: 25-50 inches - oak woodlands and prairie margins
- **Eastern Cascades**: 15-25 inches - dry pine/fir transition forests
- **Rain Shadow**: <15 inches - shrub-steppe and grasslands

### Integrating Precipitation Constraints

#### Loading Precipitation Data
```python
# Load precipitation raster (already aligned to PNW extent)
precip_path = "outputs/climate/pnw_precipitation_annual_normals.tif"
with rasterio.open(precip_path) as src:
    precip_data = src.read(1)
    precip_transform = src.transform
    
# Convert mm to inches for intuitive thresholds
precip_inches = precip_data / 25.4

# Apply precipitation constraints
precip_mask = (precip_inches >= MIN_PRECIPITATION_IN) & 
              (precip_inches <= MAX_PRECIPITATION_IN)
```

#### Example: Olympic Rainforest Bioregion
```python
# Combine elevation, tree cover, and extreme precipitation
rainforest_mask = (
    (elev_data_ft >= 0) & (elev_data_ft <= 3000) &
    (tcc_data >= 50) &
    (precip_inches >= 150)
)
```

#### Example: Eastern Cascades Dry Forest
```python
# Use precipitation to refine rain shadow boundaries
dry_forest_mask = (
    (elev_data_ft >= 1000) & (elev_data_ft <= 5000) &
    (tcc_data >= 20) &
    (precip_inches >= 15) & (precip_inches <= 30)
)
```

### Pre-Processed Precipitation Resources

1. **Precipitation Raster**: `outputs/climate/pnw_precipitation_annual_normals.tif`
   - PRISM 30-year normals (1981-2010)
   - 120m resolution matching elevation data
   - Values in millimeters (divide by 25.4 for inches)

2. **Precipitation Zones**: `outputs/climate/pnw_precipitation_zones.geojson`
   - Pre-classified precipitation categories
   - Useful for validation and visualization

3. **Mapbox Tiles**: `outputs/climate/pnw_precipitation_zones.mbtiles`
   - Ready for web mapping overlay
   - Max zoom level 14

### Processing Script
Use `scripts/process_prism_precipitation_pnw.py` to regenerate precipitation data from raw PRISM files if needed.

## Future Enhancements

### Potential Improvements
1. **Morphological operations** for smoother boundaries
2. **Multi-criteria decision analysis** for complex constraints
3. **Automated parameter optimization** based on validation metrics
4. **Additional climate variables** (temperature, fog frequency, snow depth)

### Scalability Considerations
- **Memory management** for larger rasters
- **Parallel processing** for multiple bioregions
- **Database optimization** for species validation queries
- **Automated combination** triggered by new bioregion creation

---

*This guide is based on the successful implementation of the coastal forest bioregion and should be updated as new bioregions are created and lessons learned.*