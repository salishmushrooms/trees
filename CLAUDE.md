# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a **Pacific Northwest Forest Habitat Mapping Project** that uses FIA (Forest Inventory and Analysis) survey data to create tree species habitat range and relative abundance visualizations. The project analyzes forest composition data from Washington, Oregon, and Idaho to generate species distribution maps for use in web mapping applications via Mapbox.

## The Team
I am a single developer working as part of a single person company. I prefer simple methods for achieving our goals and want to avoid overengineering or creating complicated solutions with many configuration options that would be unlikely to be used.

We will use this directory to coordinate communications among our related projects
/Users/JJC/communications/tree-climate-mushroom

## Simplify when possible
Coordinates may use a precision of 10m max unless stated otherwise
Mass, area, volume may be in whole numbers
Simplify vectors where possible but ask when level of optimal simplification is ambiguous

## Core Architecture

### Data Pipeline Structure
The project follows a multi-stage data processing pipeline:

1. **Raw Data Sources** (`data/raw/`):
   - FIA forest inventory databases (SQLite): 
     - `trees_SQLite_FIADB_WA.db` (Washington)
     - `trees_SQLite_FIADB_OR.db` (Oregon)
     - `trees_SQLite_FIADB_ID.db` (Idaho)
   - Environmental rasters: `nlcd_tcc_conus_2021_v2021-4.tif` (NLCD tree cover)
   - Hydrography data: National Hydrography Dataset shapefiles

2. **Data Extraction** (`scripts/data_extraction/`):
   - Database interface: `fia_database_interface.py`
   - Comprehensive forest data: `extract_comprehensive_forest_data.py`
   - Mushroom observations: `explore_gbif_mushroom_data.py`
   - Species-specific extraction: `extract_douglas_fir_carbon.py`

3. **Spatial Processing** (`scripts/visualization/`):
   - Species distribution mapping: `create_species_distribution_maps_cached.py`
   - Bulk processing: `bulk_process_staged_species.py`
   - Shell processing: `process_staged_species.sh`

4. **Output Generation** (`outputs/`):
   - GeoJSON vector data: `outputs/geojson/`
   - CSV tabular data: `outputs/csv/`
   - Mapbox tiles: `outputs/mapbox_mbtiles/`
   - Carbon percentile data: `outputs/plot_carbon_percentiles_latest_surveys.geojson`

### Multi-Scale Spatial Analysis Framework

The project uses a three-tier spatial analysis approach:

1. **Subplot Level (24ft radius)**: High-precision training data using exact FIA subplot boundaries
2. **Plot Level (120ft radius)**: Regional forest type classification using plot-level aggregates  
3. **Landscape Level (1-10km)**: Broad habitat suitability using environmental covariates

## Key Commands

### Data Extraction
```bash
# Extract comprehensive forest data from FIA database
cd scripts/data_extraction
python extract_comprehensive_forest_data.py

# Get mushroom observation data from GBIF
python explore_gbif_mushroom_data.py

# Extract species-specific carbon data
python extract_douglas_fir_carbon.py
```

### Regional Tree Cover Mask Setup (Recommended)
```bash
# Create reusable regional masks for faster processing
cd scripts
python create_regional_tcc_mask.py --threshold 5 --resolution 30 --min-area 2.0

# List existing masks
python create_regional_tcc_mask.py --list-existing
```

### Species Distribution Workflow

#### Step 1: Generate EDIT_ME Files
```bash
# Generate initial distribution maps for manual editing in QGIS
cd scripts/visualization
python create_species_distribution_maps_cached.py --species douglas-fir --workflow generate
```

#### Step 2: Manual Editing in QGIS
1. Load `outputs/plot_carbon_percentiles_latest_surveys.geojson` into QGIS
2. Adjust centroid sizes based on carbon percentile values
3. Load `{species}_merged_habitat_{buffer}_EDIT_ME.geojson`
4. Manually draw habitat ranges based on plot data
5. Save edited file as `{species}_merged_habitat_{buffer}_edited.geojson`

#### Step 3: Process Edited Files
```bash
# Process a single edited species with elevation and tree cover constraints
python create_species_distribution_maps_cached.py --species black-cottonwood --workflow staged

# Bulk process all edited species files
cd scripts
python bulk_process_staged_species.py
# OR use the shell script version:
./process_staged_species.sh
```

### Mapbox Tile Generation
```bash
# Create combined mbtiles file with all species layers
cd scripts/visualization
./create_mapbox_tilesets.sh

# Output will be in outputs/mapbox_mbtiles/
# Upload to Mapbox Studio for web mapping integration
```

## Database Schema

### FIA Database Structure
The project uses SQLite FIA databases with key tables:
- `PLOT`: Plot locations and metadata (LAT, LON, ELEV, INVYR)
- `COND`: Condition/forest type data (FORTYPCD, SLOPE, ASPECT, CANOPY_COVER)
- `TREE`: Individual tree records (SPCD, DIA, CARBON_AG)
- `SUBPLOT`: Subplot locations for precise spatial analysis

### Species Code Reference
Key species codes used throughout the codebase:
- Douglas-fir: SPCD 202
- Western Hemlock: SPCD 263  
- Western Red Cedar: SPCD 242
- Sitka Spruce: SPCD 098
- (See `docs/PNW_SPECIES_REFERENCE.md` for complete list)

## Coordinate Systems

- **Primary Analysis**: Albers Equal Area (EPSG:5070) for area-accurate calculations
- **Data Integration**: WGS84 (EPSG:4326) for combining disparate data sources
- **Web Mapping**: Web Mercator (EPSG:3857) for Mapbox tile outputs

## File Naming Conventions

### Scripts
- Pattern: `verb_object_method.py`
- Examples: `extract_douglas_fir_carbon.py`, `create_species_distribution_maps_cached.py`

### Data Outputs
- Pattern: `species_data_region_resolution_method.format`
- Examples: `douglas_fir_carbon_wa_1000m_r2km.tif`, `bigleaf_maple_plots.geojson`

## Development Standards

### Data Quality Requirements
- **Spatial accuracy**: ≤100m for training data, ≤1km for covariate data
- **Temporal consistency**: Use 2010-2020 data when possible
- **Coordinate precision**: Minimum 4 decimal places for lat/lon. Maximum ~10m precision

### Processing Standards
- **Memory efficiency**: Use chunked processing for large rasters (>1GB)
- **Caching**: Cache intermediate results in `cache/` directory
- **Reproducibility**: All processing parameters documented in script headers
- **Regional TCC Masks**: Use shared regional tree cover masks for 90% faster processing

### Environment Requirements
```python
# Core dependencies
import geopandas
import pandas  
import sqlite3
import shapely
import rasterio
import numpy
import scipy
```

## Common Workflows

### Adding New Species Analysis
1. Add species code to `docs/PNW_SPECIES_REFERENCE.md`
2. Update database extraction queries in `extract_comprehensive_forest_data.py`
3. Create species-specific distribution script in `scripts/visualization/`
4. Generate habitat maps using `create_species_distribution_maps_cached.py`

### Processing New Geographic Region
1. Obtain FIA database for target state/region
2. Place in `data/raw/` following naming convention
3. Update database connection paths in extraction scripts
4. Run comprehensive data extraction pipeline
5. Update regional masks in `cache/species_masks/`

## Current Status

### Completed Components
- **FIA Data Integration**: Complete databases for Washington, Oregon, and Idaho
- **Species Distribution Maps**: 18 tree species processed with habitat range mapping:
  - bigleaf maple, black cottonwood, douglas fir, engelmann spruce, grand fir
  - lodgepole pine, mountain hemlock, noble fir, pacific madrone, pacific silver fir
  - ponderosa pine, quaking aspen, red alder, subalpine fir, tanoak
  - western hemlock, western juniper, white fir
- **QGIS Workflow**: Established manual editing pipeline for habitat refinement
- **Mapbox Integration**: Combined mbtiles generation for all species layers
- **Carbon Analysis**: Plot-level carbon percentile calculations across all species

### Active Workflows
- **Manual Habitat Editing**: Ongoing refinement of species ranges using QGIS
- **Elevation/Tree Cover Constraints**: Applied to all processed species
- **Web Map Development**: WordPress integration with Mapbox layer toggles

## Documentation Structure

- `README.md`: Project overview and quick start guide
- `docs/PROJECT_PLAN.md`: Comprehensive project planning and methodology
- `docs/MUSHROOM_HABITAT_MODELING_PLAN.md`: Detailed modeling framework
- `docs/FIA_DATABASE_STRUCTURE_REFERENCE.md`: Database schema and query guide
- `docs/PNW_SPECIES_REFERENCE.md`: Complete species codes and query examples
- `MAPBOX_MBTILES_README.md`: Mapbox tile generation and upload instructions

