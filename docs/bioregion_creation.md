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

#### Key Functions (Latest Version)
1. **Boundary Loading**: Load user-defined base boundary shapefile
2. **Dynamic NLCD Clipping**: Create cached boundary-specific NLCD subset from full CONUS dataset
3. **Raster Constraints**: Apply elevation and tree cover filters with proper CRS handling
4. **Polygon Extraction**: Convert raster mask to vector polygons
5. **Two-Step Clipping**: Remove edge artifacts using precise boundary intersection + inward buffer
6. **Multi-Polygon Output**: Keep individual polygons for performance (no merging)

### 3. Configuration Parameters

#### Typical Bioregion Parameters (Latest Version)
```python
# Bioregion-specific constraints (Urban Agricultural Example)
MIN_ELEVATION_FT = 1           # Minimum elevation
MAX_ELEVATION_FT = 1500        # Maximum elevation for urban/agricultural areas
MIN_TREE_COVER_PCT = 0         # Minimum tree cover (urban/agricultural)
MAX_TREE_COVER_PCT = 50        # Maximum tree cover (not dense forest)
OUTPUT_RESOLUTION = 120        # Match elevation data resolution

# Data source optimization
TREE_COVER_RASTER = "data/raw/nlcd_tcc_conus_2021_v2021-4.tif"  # Full CONUS dataset
# Automatically creates: "outputs/mapbox_masks/nlcd_tcc_urban_agricultural_extent.tif"

# Performance optimization settings
min_area_sqkm = 0.5           # Larger minimum for fewer polygons (was 0.08)
inward_buffer = 0             # No pre-buffer (use two-step clipping instead)
skip_merging = True           # Keep individual polygons for performance

# *** NEW: Geometry Simplification Settings (Latest 2025) ***
simplification_tolerance = 0.001      # ~100m simplification for web performance
coordinate_precision_decimal_places = 3  # 3 decimal places = ~100m precision
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

#### Phase 3: Polygon Processing (Latest Approach)
1. **Convert raster to polygons** using `rasterio.features.shapes()`
2. **Two-step clipping** for edge artifact removal:
   - Step 1: Clip to original boundary (creates partial pixel artifacts)
   - Step 2: Apply 50m inward buffer and re-clip (removes edge artifacts)
3. **Multi-stage filtering**:
   - 2000 sq meter minimum during clipping
   - Single-pixel removal (< 1.5 × pixel area)
   - Rectangular strip removal (high aspect ratio + thin width)
   - Final minimum area filtering (0.5 km² typical)
4. **Geometry simplification and precision reduction** (NEW 2025):
   - Simplify geometry with 100m tolerance (0.001 degrees)
   - Round coordinates to 3 decimal places (~100m precision)
   - Reduces file size by ~90% while maintaining regional accuracy
5. **Skip merging** for performance (keep individual polygons)

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

### 1.5. NLCD CRS Updates (January 2025) - **CRITICAL**
**Problem**: The updated NLCD tree cover dataset (`nlcd_tcc_conus_2021_v2021-4.tif`) uses a different coordinate reference system (likely Albers Equal Area Conic) instead of the previously assumed WGS84.

**Symptoms**:
- CRS-related errors during raster masking operations
- Failed reprojections when combining elevation and tree cover data  
- Boundary geometry incompatibility with tree cover raster

**Scripts Updated** (January 2025):
- ✅ `create_urban_agricultural_shapefile_constrained.py` (reference implementation)
- ✅ `create_eastern_cascades_shapefile_constrained.py`
- ✅ `create_high_cascades_shapefile_constrained.py`  
- ✅ `create_coastal_forest_shapefile_constrained.py`
- ✅ `create_olympic_mountains_shapefile_constrained.py`

**Solution Implemented**: All scripts now include automatic CRS detection and reprojection:

```python
# Enhanced create_clipped_nlcd_if_needed() function:
with rasterio.open(TREE_COVER_RASTER) as src:
    print(f"NLCD CRS: {src.crs}")
    print(f"Boundary CRS: {boundary_gdf.crs}")
    
    # Reproject boundary to match NLCD CRS
    boundary_reproj = boundary_gdf.to_crs(src.crs)
    boundary_geom_reproj = boundary_reproj.geometry.iloc[0]
    boundary_buffered_reproj = boundary_geom_reproj.buffer(1000)  # 1000m buffer in projected units

# Enhanced apply_raster_constraints_within_boundary() function:
with rasterio.open(tcc_path) as tcc_src:
    # Handle CRS mismatch - reproject boundary to match tree cover data if needed
    if tcc_src.crs != elev_crs:
        print(f"Reprojecting boundary from {elev_crs} to {tcc_src.crs} for tree cover processing")
        boundary_gdf_tcc = gpd.GeoDataFrame([1], geometry=[boundary_geom], crs=elev_crs)
        boundary_gdf_tcc = boundary_gdf_tcc.to_crs(tcc_src.crs)
        boundary_geom_tcc = boundary_gdf_tcc.geometry.iloc[0]
    else:
        boundary_geom_tcc = boundary_geom
        
    # Use the correctly reprojected boundary for all tree cover operations
    tcc_masked, _ = mask(tcc_src, [boundary_geom_tcc], crop=True, nodata=tcc_src.nodata)
```

**⚠️ All bioregion creation scripts now handle NLCD CRS automatically - no manual intervention required!**

### 2. Edge Artifacts (SOLVED - Latest Approach)
**Problem**: Thin rectangular strips appear along boundary edges due to raster/vector misalignment
**Symptoms**: 
- ~500m thick border around entire perimeter
- Single-pixel wide rectangles (~72-120m wide)
- Perfect rectangular shapes following raster grid
- Areas as small as 7,000-35,000 m²

**Root Cause**: Raster pixels don't align perfectly with vector boundaries, creating partial pixels that become thin polygons when converted to vectors.

**Latest Solution - Two-Step Clipping (Urban Agricultural Script v2024)**:
```python
def clip_to_original_boundary(gdf, boundary_gdf):
    # Step 1: Initial clipping to original boundary
    # (Creates partial pixel artifacts)
    
    # Step 2: Inward buffer and re-clip to remove edge artifacts
    buffer_degrees = -0.00045  # ~50m inward buffer
    buffered_boundary = original_boundary.buffer(buffer_degrees)
    
    # Re-clip to buffered boundary removes all edge artifacts
    final_clipped_geom = geom.intersection(buffered_boundary)
```

**Advantages of Two-Step Approach**:
- **Eliminates buffer guesswork**: No need to pre-determine buffer size
- **Preserves real data**: Only removes problematic edge artifacts
- **Cleaner results**: No ~500m borders around boundaries
- **Faster processing**: Eliminates need for complex merging operations

### 4. Geometry Simplification Standard (NEW 2025)
**Problem**: High-precision geometries create massive file sizes (50+ MB) and cause Mapbox upload failures
**Solution**: Standardized 100m precision across all bioregion scripts

#### Implementation Code Template
```python
# Simplify geometry for web use (consistent with bioregion combination)
tolerance = 0.001  # ~100m simplification
large_polygons['geometry'] = large_polygons.geometry.simplify(tolerance)

# Round coordinates to 3 decimal places (~100m precision)
from shapely.ops import transform
def round_coordinates(geom):
    """Round coordinates to 3 decimal places (~100m precision)"""
    def round_coords(x, y, z=None):
        rounded_x = round(x, 3)
        rounded_y = round(y, 3)
        return (rounded_x, rounded_y) if z is None else (rounded_x, rounded_y, round(z, 3))
    return transform(round_coords, geom)

large_polygons['geometry'] = large_polygons['geometry'].apply(round_coordinates)

# Update metadata
summary['parameters'].update({
    'simplification_tolerance': 0.001,
    'coordinate_precision_decimal_places': 3
})
```

#### Apply to All Scripts
**Required Updates**: Apply this pattern to all bioregion creation scripts:
- `create_coast_range_shapefile_constrained.py` ✅
- `create_western_cascades_shapefile_constrained.py` ✅  
- `create_urban_agricultural_shapefile_constrained.py` ✅
- `create_eastern_cascades_shapefile_constrained.py` ✅ (Updated Jan 2025)
- `create_high_cascades_shapefile_constrained.py` ✅ (Updated Jan 2025)
- `create_olympic_mountains_shapefile_constrained.py` ✅ (Updated Jan 2025)
- `create_coastal_forest_shapefile_constrained.py` ✅ (Updated Jan 2025)
- Any other bioregion creation scripts

#### Benefits
- **File size reduction**: 90%+ smaller files (53MB → 5MB typical)
- **Faster processing**: Quicker uploads, downloads, and rendering
- **Mapbox compatibility**: Eliminates "vector layers" metadata errors
- **Consistent precision**: Uniform 100m accuracy across all bioregions
- **Web performance**: Optimal for regional-scale mapping applications

#### Precision Justification
100m precision is optimal for bioregion mapping because:
- **Ecological relevance**: Bioregions are broad habitat zones, not precise boundaries
- **Data source alignment**: Matches 120m elevation data resolution
- **Web mapping standards**: Appropriate for zoom levels 0-13
- **Processing efficiency**: Balances accuracy with performance

### 2. NLCD Tree Cover Data Handling (Latest 2025 Approach)
**Problem**: Using regional NLCD subset (WA/OR only) limits coverage and requires manual file management. Additionally, the updated NLCD dataset uses a different CRS than originally expected.

**Latest Solution - Dynamic NLCD Clipping with CRS Handling**:
```python
def create_clipped_nlcd_if_needed(boundary_gdf):
    # Use full CONUS NLCD dataset: data/raw/nlcd_tcc_conus_2021_v2021-4.tif
    # Automatically handle CRS differences and clip to boundary extent
    with rasterio.open(TREE_COVER_RASTER) as src:
        print(f"NLCD CRS: {src.crs}")
        print(f"Boundary CRS: {boundary_gdf.crs}")
        
        # Reproject boundary to match NLCD CRS
        boundary_reproj = boundary_gdf.to_crs(src.crs)
        boundary_geom_reproj = boundary_reproj.geometry.iloc[0]
        boundary_buffered_reproj = boundary_geom_reproj.buffer(1000)  # 1000m buffer in projected units
        
        # Mask and crop to boundary
        clipped_data, clipped_transform = mask(
            src, 
            [boundary_buffered_reproj], 
            crop=True, 
            nodata=src.nodata
        )
```

**Critical CRS Handling in Processing**:
```python
# In apply_raster_constraints_within_boundary():
with rasterio.open(tcc_path) as tcc_src:
    # Handle CRS mismatch - reproject boundary to match tree cover data if needed
    if tcc_src.crs != elev_crs:
        print(f"Reprojecting boundary from {elev_crs} to {tcc_src.crs} for tree cover processing")
        boundary_gdf_tcc = gpd.GeoDataFrame([1], geometry=[boundary_geom], crs=elev_crs)
        boundary_gdf_tcc = boundary_gdf_tcc.to_crs(tcc_src.crs)
        boundary_geom_tcc = boundary_gdf_tcc.geometry.iloc[0]
    else:
        boundary_geom_tcc = boundary_geom
        
    # Use the correctly reprojected boundary for all tree cover operations
    tcc_masked, _ = mask(tcc_src, [boundary_geom_tcc], crop=True, nodata=tcc_src.nodata)
```

**Benefits**:
- **Complete coverage**: Uses full CONUS dataset ensuring no gaps
- **Automatic CRS detection**: Detects coordinate system mismatches and reprojects as needed
- **Robust processing**: Handles NLCD data in Albers Equal Area or other projected coordinate systems
- **Performance optimization**: Creates boundary-specific cached files
- **Efficient processing**: Only clips once, reuses cached version

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

## Validation Approach (Ignored for now)

### Species Validation (Ignored for now)
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

## Latest Performance Improvements (Urban Agricultural Script v2024)

The urban agricultural bioregion script represents the most advanced version of the bioregion creation pipeline, incorporating significant performance and accuracy improvements:

### Key Improvements

#### 1. Two-Step Edge Artifact Removal
**Problem Solved**: Previous approach used pre-buffers that caused data gaps and still left edge artifacts
**Solution**: Two-step clipping process that eliminates the root cause of edge artifacts
```python
# Step 1: Clip to exact boundary (creates partial pixel artifacts)
# Step 2: Apply 50m inward buffer and re-clip (removes all edge artifacts)
buffer_degrees = -0.00045  # ~50m at PNW latitude
```

#### 2. Dynamic NLCD Data Management  
**Problem Solved**: Regional NLCD subsets require manual file management and may have coverage gaps
**Solution**: Automatic clipping from full CONUS dataset with caching
```python
def create_clipped_nlcd_if_needed(boundary_gdf):
    # Uses full CONUS dataset, clips once, caches for reuse
    # Handles CRS reprojection automatically
```

#### 3. Performance Optimization
**Problem Solved**: Polygon merging with 8,000+ polygons caused excessive processing time (hours)
**Solution**: Skip merging, output individual polygons
- **Before**: Hours of processing for `unary_union()` on 8,774 polygons  
- **After**: Seconds to complete, individual polygons work perfectly in GIS/web mapping

#### 4. Enhanced Filtering
**Multi-stage approach** eliminates artifacts more effectively:
```python
# During clipping: 2000 sq meter minimum
# Post-processing: Single-pixel removal, edge artifact removal, final area filter
```

### Script Comparison: Before vs After

| Aspect | Original Approach | Latest Approach (Urban Agricultural) |
|--------|-------------------|-------------------------------------|
| **Edge Artifacts** | Pre-buffer + post-filtering | Two-step clipping |
| **NLCD Data** | Regional subset files | Dynamic CONUS clipping |
| **Performance** | Hours (merging bottleneck) | Seconds (skip merging) |
| **Polygon Count** | 1 merged polygon | Individual polygons |
| **CRS Handling** | Manual coordination | Automatic detection/reprojection |
| **Caching** | No caching | Automatic NLCD caching |
| **Artifact Removal** | Single-pass filtering | Multi-stage progressive filtering |

### Recommended Migration

For new bioregions, use the urban agricultural script pattern:
1. **Copy** `scripts/create_urban_agricultural_shapefile_constrained.py`
2. **Modify** constraint parameters for your target bioregion
3. **Update** file paths and output names
4. **Retain** the two-step clipping and performance optimizations

## Best Practices

### 1. Incremental Development
- Start with basic tree cover + elevation constraints only
- Add complexity (water exclusion, shape filters) incrementally
- Test each addition separately to isolate issues

### 2. Parameter Tuning
- Start with permissive parameters and tighten as needed
- Use visual inspection in QGIS to validate results
- Adjust minimum area thresholds based on observed artifacts

### 3. Geometry Optimization (NEW 2025)
- **Apply simplification early**: Use 100m tolerance (0.001 degrees) in individual scripts
- **Standardize precision**: Round coordinates to 3 decimal places (~100m precision)
- **Document settings**: Include simplification parameters in metadata
- **Test file sizes**: Aim for <10MB for combined regional datasets

### 4. Documentation
- Record all parameter changes and their effects
- Document species-specific ecological rationale
- Maintain consistent file naming conventions

### 5. Validation (Ignored for now)
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
outputs/bioregions/all_bioregions_combined.geojson         # Combined dataset (~100m precision)
outputs/bioregions/bioregions_combination_summary.json     # Statistics & metadata
outputs/mapbox_mbtiles/pnw_forest_bioregions.mbtiles       # Mapbox tiles (max zoom 13)
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
- Max zoom level: 13 (suitable for regional analysis, reduced for performance)
- Includes properties: `region_name`, `region_code`, `area_km2`, `bioregion_type`
- Geometry precision: ~100m (optimized for web performance)
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