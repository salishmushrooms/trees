# Pacific Northwest Forest and Mushroom Habitat Modeling Project

## Project Overview
**Primary Goal**: Create predictive models for mushroom fruiting patterns across the Pacific Northwest using forest inventory data, environmental covariates, and mushroom observation records.

**Evolution**: This project began as Douglas-fir carbon mapping and has evolved into comprehensive mushroom habitat modeling for the Pacific Northwest.

## Project Organization

### Directory Structure
```
trees/
├── docs/                           # All documentation
│   ├── PROJECT_PLAN.md
│   ├── MUSHROOM_HABITAT_MODELING_PLAN.md
│   └── FIA_PLOT_STRUCTURE.md
├── data/
│   ├── raw/                        # Original source data
│   │   ├── trees_SQLite_FIADB_WA.db
│   │   ├── trees_SQLite_FIADB_OR.db
│   │   ├── nlcd_tcc_conus_2021_v2021-4.tif
│   │   └── *.pdf (FIADB guides)
│   └── processed/                  # Cleaned/processed data
├── scripts/
│   ├── data_extraction/           # Database queries & data extraction
│   ├── analysis/                  # Spatial analysis & processing
│   └── raster_creation/           # Raster generation scripts
└── outputs/
    ├── rasters/                   # Generated raster products
    ├── geojson/                   # Spatial vector outputs
    └── csv/                       # Tabular data outputs
```

### File Naming Conventions
- **Raw data**: Original filenames preserved
- **Scripts**: `verb_object_method.py` (e.g., `extract_douglas_fir_carbon.py`)
- **Outputs**: `species_data_wa_resolution_method.format` (e.g., `douglas_fir_carbon_wa_1000m_r2km.tif`)
- **Documentation**: `TOPIC_DESCRIPTION.md` (all caps for major docs)

## Target Species (Updated)

### Chanterelles
- *Cantharellus formosus* (Pacific Golden Chanterelle)
- *Cantharellus cascadensis* (Cascade Chanterelle)  
- *Cantharellus subalpinus* (Subalpine Chanterelle)

### Boletes
- *Boletus edulis* (King Bolete)
- *Boletus rex-veris* (Spring King)
- *Boletus barrowsii* (White King)
- *Boletus regineus* (Queen Bolete)

### Morels (Expanded List)
- *Morchella importuna* (Landscape Morel - urban associations)
- *Morchella snyderi* (Snyder's Morel)
- *Morchella tridentina* 
- *Morchella americana* (American Morel)
- *Morchella norvegienses* (Norwegian Morel - broad PNW distribution)
- *Morchella populiphila* (Cottonwood Morel)
- *Morchella tomentosa* (Gray Morel - burn specialist)
- *Morchella exuberans* (Burn Morel)

### Matsutake Complex
- *Tricholoma magnivelare* (American Matsutake)
- *Tricholoma caligatum* (False Matsutake)

## Completed Work

### Comprehensive Forest Data Extraction ✅
**Status**: **COMPLETED** - We have extracted comprehensive forest data across all Washington FIA plots

**Script**: `scripts/data_extraction/extract_comprehensive_forest_data.py`
**Output**: `outputs/comprehensive_forest_data_wa.geojson` (15MB)

**Results**:
- **12,878 forested plots** across Washington State
- **Comprehensive species composition**: Douglas-fir (74.9% of plots), Western Hemlock (39.4%), Western Red Cedar (26.6%), Hardwoods (32.6%)
- **Forest type classification**: Douglas-fir Dominant (34.0%), Mixed Conifer (17.1%), Fir Species Dominant (12.6%)
- **Topographic data**: Elevation 3-7,900 feet, mean slope 35.0°, aspect distribution
- **Carbon data**: Above-ground carbon by species
- **Diversity metrics**: Species count per plot, basal area by species

### GBIF Mushroom Data Integration ✅
**Status**: **COMPLETED** - Extracted and filtered mushroom observation data

**Script**: `scripts/data_extraction/explore_gbif_mushroom_data.py`
**Outputs**: 
- `outputs/gbif_mushroom_data_raw.csv` (1,826 records)
- `outputs/gbif_mushroom_data_filtered.csv` (977 high-quality records ≤1km uncertainty)
- `outputs/gbif_species_summary.csv`

**Key Findings**:
- **Seasonal patterns**: Peak fruiting in September (228 obs) and April-May (274 combined)
- **Species distribution**: Chanterelles (340 records), Boletes (339), Morels (297), Matsutake (1)
- **Geographic coverage**: Washington, Oregon, Idaho

### Douglas-fir Carbon Raster Pipeline ✅  
**Status**: **COMPLETED** - 6 raster products created

**Scripts**: Multiple scripts in `scripts/analysis/`
**Outputs**: `outputs/rasters/douglas_fir_carbon_wa_*.tif`

**Products**: 500m & 1000m resolution × 1km, 2km, 3km search radii
**Methodology**: Inverse Distance Weighting (IDW) with data-constrained approach
**Recommended**: `douglas_fir_carbon_wa_1000m_r2km.tif` (21.1% coverage, balanced)

## Next Priority Tasks

### 1. Subplot-Level Mapping (Addresses User Question)
**Objective**: Create precise subplot boundary mapping for model training

**Approach**:
- **24-foot radius circles**: Map each subplot as 24ft radius circle (exact FIA protocol)
- **48-foot square pixels**: Alternative approach with 48ft × 48ft pixels (0.053 acres ≈ subplot area)
- **Subplot composition**: Extract tree species data at subplot level for highest spatial precision

**Implementation**:
```python
# Proposed script: scripts/analysis/create_subplot_boundaries.py
# - Query FIADB SUBPLOT table for precise coordinates
# - Create 24ft radius buffers around each subplot center
# - Extract tree composition for each subplot individually
# - Generate high-resolution (3-10m) raster products
```

**Benefits**:
- **Training data precision**: Exact forest composition for ML model training
- **Spatial accuracy**: Match remote sensing pixel scale (10-30m)
- **Validation capability**: Ground-truth data for model validation

### 2. Plot-Level Aggregate Mapping
**Objective**: Create plot-scale habitat suitability mapping

**Approach**:
- **120-foot radius circles**: Approximate total plot area (4 subplots + connecting areas)
- **Aggregate composition**: Plot-level forest type classification
- **Medium resolution**: 30-100m pixels for landscape-scale modeling

### 3. Waterways Data Integration (User Request)
**Priority**: **HIGH** - Critical for mushroom habitat modeling

**Data Sources**:
- **National Hydrography Dataset (NHD)**: Streams, rivers, water bodies
- **USGS Stream Network**: Flow characteristics, seasonal patterns
- **Watershed boundaries**: HUC (Hydrologic Unit Code) boundaries

**Mushroom-Water Associations**:
- **Stream proximity**: Distance to nearest water feature
- **Moisture gradients**: Areas with consistent soil moisture
- **Riparian zones**: Unique forest compositions near water
- **Seasonal flooding**: Areas with spring moisture retention

**Implementation Priority**:
```python
# Proposed script: scripts/data_extraction/extract_waterways_data.py
# - Download NHD data for Washington State
# - Calculate distance-to-water for all FIA plots
# - Create moisture-proximity covariates
# - Analyze mushroom observations vs. water proximity
```

### 4. Environmental Covariate Assembly
**Next Phase**: Integrate comprehensive environmental predictors

**Priority Data**:
- **Climate**: PRISM temperature/precipitation (monthly normals)
- **Topography**: Slope, aspect, topographic wetness index
- **Soils**: SSURGO soil types, pH, drainage class
- **Remote sensing**: NDVI time series, forest type classifications

## Methodological Framework

### Multi-Scale Approach
1. **Subplot level (24ft radius)**: Precise species composition for training
2. **Plot level (120ft radius)**: Forest type classification for regional mapping  
3. **Landscape level (1-10km)**: Broad habitat suitability predictions

### Spatial Resolution Strategy
- **Fine scale (10-30m)**: Subplot-based training data, matches Landsat/Sentinel-2
- **Medium scale (100-500m)**: Plot-based regional mapping
- **Coarse scale (1-5km)**: Climate and broad landscape patterns

### Model Training Philosophy
- **Ground truth**: Use subplot data as precise training points
- **Validation**: Plot-level aggregates for model validation
- **Prediction**: Scale up to landscape using environmental covariates

## Technical Standards

### Coordinate Systems
- **Primary projection**: Albers Equal Area (EPSG:5070) for area-accurate calculations
- **Working projection**: WGS84 (EPSG:4326) for data integration
- **Local projections**: State-specific systems when needed for accuracy

### Data Quality Standards
- **Spatial accuracy**: ≤100m for training data, ≤1km for covariate data
- **Temporal consistency**: Use 2010-2020 data when possible
- **Validation requirements**: Hold-out datasets for all models

### Processing Standards
- **Memory efficiency**: Chunked processing for large rasters
- **Reproducibility**: All parameters documented in scripts
- **Version control**: Date-stamped outputs with methodology notes

## Current Data Inventory

### ✅ Available Data
- **WA FIA Database**: 12,878 forested plots with comprehensive species data
- **OR FIA Database**: 3.3GB database (not yet processed)
- **GBIF Mushroom Data**: 977 high-quality observations across target species
- **NLCD Tree Cover**: 2021 data for Pacific Northwest
- **Douglas-fir Rasters**: 6 products at different resolutions/search radii

### 🔄 In Development
- **Waterways data**: National Hydrography Dataset integration
- **Subplot boundaries**: 24ft radius precise mapping
- **Oregon FIA data**: Extraction and integration with WA data

### 📋 Needed Data
- **Climate data**: PRISM monthly normals
- **Soil data**: SSURGO database for PNW  
- **Elevation**: High-resolution DEM for topographic variables
- **Remote sensing**: Landsat/Sentinel time series for forest monitoring

## Research Questions

### Ecological Questions
1. **Habitat specificity**: Which forest types are most associated with each mushroom species?
2. **Water dependency**: How important is proximity to water features for fruiting?
3. **Seasonal patterns**: Can we predict fruiting timing based on weather/climate?
4. **Scale effects**: Do habitat relationships change at different spatial scales?

### Technical Questions  
1. **Optimal resolution**: What pixel size best captures mushroom habitat relationships?
2. **Model performance**: How well can forest inventory data predict mushroom habitat?
3. **Temporal stability**: Are habitat models stable across different years?
4. **Transferability**: Do models trained in WA work in OR/ID?

## Success Metrics
1. **Data completeness**: >90% spatial coverage of target region with environmental covariates
2. **Model performance**: AUC >0.8 for species habitat suitability models
3. **Spatial validation**: Correspondence with known mushroom collection sites
4. **Expert validation**: Review and approval by professional mycologists
5. **Application success**: Useful predictions for citizen scientists and researchers

---
*Last updated: December 2024*
*Project status: Active development - Comprehensive forest data extraction complete* 