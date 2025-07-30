# Pacific Northwest Forest Habitat Mapping Project

## Overview
This project uses FIA (Forest Inventory and Analysis) survey data to create tree species habitat range and relative abundance visualizations across the Pacific Northwest. The system generates species distribution maps that can be displayed as interactive layers in web mapping applications via Mapbox.

**Primary Purpose**: Query FIA databases and build geospatial visualizations of tree species distributions  
**Geographic Scope**: Washington, Oregon, and Idaho  
**Output Products**: Species habitat ranges, relative abundance maps, and Mapbox-ready tile layers  

## Project Organization

```
trees/
├── README.md                       # This file - Project overview
├── CLAUDE.md                       # Complete project guide & workflows
├── docs/                           # Technical documentation
│   ├── FIA_DATABASE_STRUCTURE_REFERENCE.md  # Database schema reference
│   ├── PNW_SPECIES_REFERENCE.md             # Species codes and examples
│   └── QGIS_INTEGRATION_GUIDE.md            # QGIS workflow guide
├── data/raw/                       # Source data
│   ├── trees_SQLite_FIADB_{WA,OR,ID}.db     # FIA databases
│   ├── nlcd_tcc_conus_2021_v2021-4.tif     # NLCD tree cover
│   └── hydrography/                         # NHD water data
├── scripts/                        # Analysis scripts
│   ├── data_extraction/                     # Database queries
│   ├── database/                            # Database utilities
│   ├── visualization/                       # Species mapping
│   └── utilities/                           # Standalone tools
├── archive/                        # Deprecated/experimental scripts
└── outputs/                        # Generated data
    ├── species_distribution/                # Species habitat files
    ├── mapbox_mbtiles/                      # Web mapping tiles
    └── plot_carbon_percentiles_latest_surveys.geojson
```

## Current Project Status

### ✅ Completed Components
- **FIA Database Integration**: Complete databases for Washington, Oregon, and Idaho
- **Species Distribution Maps**: 18 tree species with habitat range mapping
  - bigleaf maple, black cottonwood, douglas fir, engelmann spruce, grand fir
  - lodgepole pine, mountain hemlock, noble fir, pacific madrone, pacific silver fir
  - ponderosa pine, quaking aspen, red alder, subalpine fir, tanoak
  - western hemlock, western juniper, white fir
- **QGIS Workflow**: Established manual editing pipeline for habitat refinement
- **Mapbox Integration**: Combined mbtiles generation for web mapping
- **Carbon Analysis**: Plot-level carbon percentile calculations across all species

### 🔄 Active Development
- **Manual Habitat Refinement**: Ongoing QGIS editing of species ranges
- **Web Map Integration**: WordPress/Mapbox layer toggle implementation
- **Additional Species**: Expanding to include more PNW tree species

## Quick Start

See **`CLAUDE.md`** for complete workflows and detailed implementation guidance.

### Basic Commands
```bash
# Extract forest data
cd scripts/data_extraction
python extract_comprehensive_forest_data.py

# Create species maps
cd scripts/visualization  
python create_species_distribution_maps_cached.py --species douglas-fir --workflow generate

# Generate Mapbox tiles
cd scripts
./create_mapbox_tilesets.sh
```

## Key Files
- **Main Guide**: `CLAUDE.md` - Complete project documentation
- **Species Data**: `outputs/plot_carbon_percentiles_latest_surveys.geojson`  
- **Web Tiles**: `outputs/mapbox_mbtiles/*.mbtiles`
- **Databases**: `data/raw/trees_SQLite_FIADB_{WA,OR,ID}.db`

## Technical Documentation
- `docs/FIA_DATABASE_STRUCTURE_REFERENCE.md` - Database schema reference
- `docs/PNW_SPECIES_REFERENCE.md` - Species codes and examples  
- `docs/QGIS_INTEGRATION_GUIDE.md` - QGIS workflow guide
- `MAPBOX_MBTILES_README.md` - Mapbox tile generation guide

---
*Last updated: December 2024*  
*Project status: Active - 18 species processed, ongoing habitat refinement* 