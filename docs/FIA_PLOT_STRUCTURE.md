# FIA Plot Structure and Data Organization

## Plot Overview
FIA (Forest Inventory and Analysis) plots are the fundamental sampling units used to collect forest data. Each plot represents a specific location where detailed forest measurements are taken.

### Plot Identification
- Each plot has a unique identifier (`PLOT_CN`)
- Plots are organized by:
  - State Code (`STATECD`)
  - Unit Code (`UNITCD`)
  - County Code (`COUNTYCD`)
  - Plot Number (`PLOT`)
  - Inventory Year (`INVYR`)

### Plot Location
- Plots are georeferenced using:
  - Latitude (`LAT`)
  - Longitude (`LON`)
- Coordinates are stored in WGS84 (EPSG:4326) for geographic reference

## Plot Structure

### Plot Design
- Each FIA plot consists of **four subplots** arranged in a fixed cluster pattern
- **Total plot area**: Approximately 1 acre (variable due to subplot arrangement)
- **Subplot spacing**: Center-to-center distance of 120 feet between subplots
- **Plot center**: Located at Subplot 1 (center subplot)

### Subplot Organization and Precise Mapping

#### Subplot Layout
```
    Subplot 2 (North)
         |
         | 120 ft
         |
Subplot 3 --- Subplot 1 (Center) --- [120 ft spacing]
(East)                 
         |
         | 120 ft  
         |
    Subplot 4 (South)
```

#### Subplot Specifications
- **Subplot 1**: Center subplot (Plot center coordinates)
- **Subplot 2**: 120 ft north of center (bearing 0°/360°)
- **Subplot 3**: 120 ft southeast of center (bearing 120°)
- **Subplot 4**: 120 ft southwest of center (bearing 240°)

#### Measurement Areas per Subplot
- **Microplot**: 6.8 ft radius (0.033 acres) - seedlings and saplings
- **Subplot**: **24.0 ft radius** (0.133 acres) - **main tree measurements**
- **Macroplot**: 58.9 ft radius (1.0 acres) - select measurements

### **Precise Spatial Mapping for Mushroom Habitat Modeling**

#### Subplot-Level Mapping (24ft Radius Circles)
**Optimal for**: High-precision forest composition training data

**Specifications**:
- **Circle radius**: 24 feet (7.32 meters)
- **Area coverage**: 0.133 acres (538 m²) per subplot
- **Coordinate precision**: ±10 feet from calculated subplot center
- **Resolution match**: 3-10m pixels (matches Sentinel-2, high-res Landsat)

**Implementation approach**:
```python
# Calculate subplot centers from plot center coordinates
# Subplot 1: plot_lat, plot_lon (center)
# Subplot 2: plot_lat + 120ft*cos(0°), plot_lon + 120ft*sin(0°)  
# Subplot 3: plot_lat + 120ft*cos(120°), plot_lon + 120ft*sin(120°)
# Subplot 4: plot_lat + 120ft*cos(240°), plot_lon + 120ft*sin(240°)

# Create 24ft radius buffers around each subplot center
# Extract precise tree species composition per subplot
# Generate 3-10m resolution training data
```

#### Plot-Level Mapping (120ft Radius Circles)
**Optimal for**: Regional forest type classification

**Specifications**:
- **Circle radius**: 120 feet (36.6 meters) - encompasses subplot cluster
- **Area coverage**: ~1.0 acre (4,047 m²) per plot
- **Coordinate precision**: Plot center coordinates
- **Resolution match**: 30-100m pixels (matches Landsat 8/9)

#### Alternative: Square Pixel Approach
**48ft × 48ft pixels** (14.6m × 14.6m):
- **Area**: 0.053 acres (214 m²) - close to subplot area
- **Advantages**: Aligns with regular grid systems, easier raster processing
- **Disadvantages**: May not capture circular measurement protocol exactly

### Data Collection

#### Tree Measurements (Subplot Level)
- **Individual trees** measured within each 24ft radius subplot
- **Key measurements** include:
  - Species identification (`SPCD`)
  - Diameter at breast height (`DIA`) 
  - Tree height (`HT`)
  - Tree condition and status
  - Above-ground carbon content (`CARBON_AG`)
  - Live/dead status

#### Subplot-Level Data (Available in FIADB)
- **SUBPLOT table**: Contains subplot-specific measurements
- **TREE table**: Individual tree records linked to specific subplots
- **Precise coordinates**: Can be calculated from plot center + subplot geometry

#### Plot-Level Data (Aggregated)
- **Forest type classification** (`FORTYPCD`)
- **Stand characteristics** (age, density, structure)
- **Topographic data** (slope, aspect, elevation)
- **Canopy cover** (`LIVE_CANOPY_CVR_PCT`)
- **Environmental conditions**

## Mushroom Habitat Modeling Applications

### Training Data Precision Levels

#### **Level 1: Subplot Precision (24ft radius)**
- **Use case**: Machine learning model training with maximum precision
- **Spatial scale**: 3-10m pixel resolution
- **Data quality**: Exact species composition, DBH distribution, carbon content
- **Sample size**: ~51,000 subplots across Washington (4 × 12,878 plots)

#### **Level 2: Plot Aggregation (120ft radius)**  
- **Use case**: Regional forest type mapping and validation
- **Spatial scale**: 30-100m pixel resolution
- **Data quality**: Dominant forest type, diversity indices, structural metrics
- **Sample size**: 12,878 plots across Washington

#### **Level 3: Landscape Prediction (1-10km)**
- **Use case**: Broad habitat suitability mapping
- **Spatial scale**: 1-5km pixel resolution  
- **Data quality**: Environmental covariates, climate, topography
- **Coverage**: Complete landscape coverage via environmental modeling

### Advantages of Multi-Scale Approach

1. **High-resolution training**: Subplot data provides precise forest composition for ML training
2. **Medium-resolution validation**: Plot data validates models at management-relevant scales  
3. **Low-resolution prediction**: Environmental covariates enable landscape-scale habitat maps
4. **Scale matching**: Aligns with satellite remote sensing pixel sizes (10m, 30m, 1km)

## Data Extraction Strategy

### Current Status (✅ Completed)
- **Plot-level extraction**: 12,878 plots with aggregated forest composition
- **Species composition**: All major Pacific Northwest tree species
- **Topographic data**: Slope, aspect, elevation for each plot

### Next Steps (🔄 In Development)
- **Subplot-level extraction**: Individual subplot coordinates and tree composition
- **Precise boundaries**: 24ft radius circles for each subplot
- **High-resolution rasters**: 3-10m resolution forest composition maps
- **Validation framework**: Compare subplot vs. plot vs. landscape predictions

---

## **Boletus rex-veris Correlation Analysis Application**

### **Purpose**: Correlate mushroom observations with forest composition using FIA plot data

#### **Spatial Matching Strategy**

**1. Proximity Search**
- Search radius: **2km** around each Boletus observation
- Account for FIA coordinate "fuzziness" (±0.5-1.0 mile uncertainty)
- Weight correlations by distance (closer plots = stronger signal)

**2. Elevation Matching**
- Primary tolerance: **±200 feet** elevation difference
- Secondary analysis: **±500 feet** for broader patterns
- Critical for elevation-sensitive species like Boletus rex-veris

**3. Temporal Considerations**
- Use most recent inventory cycle (2015-2020)
- Account for forest succession since measurement
- Match observation season with plot measurement timing when possible

#### **Data Variables for Analysis**

**Plot-Level Variables** (from PLOT/COND tables):
- `LAT`, `LON`: Fuzzy plot coordinates
- `ELEV`: Plot elevation (feet)
- `FORTYPCD`: Forest type classification
- `LIVE_CANOPY_CVR_PCT`: Live canopy cover percentage
- `STDAGE`: Stand age
- `SLOPE`, `ASPECT`: Topographic variables

**Tree Species Composition** (from TREE table):
- `SPCD`: Species codes for all trees on plot
- `DIA`: Tree diameter (proxy for forest structure)
- `CARBON_AG`: Above-ground carbon by species
- Species presence/absence within 2km radius
- Relative basal area and dominance

**Key Pacific Northwest Species to Analyze**:
```
202 - Douglas-fir (Pseudotsuga menziesii)
263 - Western hemlock (Tsuga heterophylla)
242 - Western red cedar (Thuja plicata)
017 - Pacific silver fir (Abies amabilis)
011 - Noble fir (Abies procera)
122 - Ponderosa pine (Pinus ponderosa)
108 - Lodgepole pine (Pinus contorta)
351 - Red alder (Alnus rubra)
```

#### **Statistical Analysis Framework**

**Correlation Metrics**:
1. **Species Co-occurrence Frequency**: % plots with species X near Boletus
2. **Distance-Weighted Association**: Closer plots weighted higher
3. **Chi-square Tests**: Statistical significance of associations
4. **Logistic Regression**: P(Boletus | Forest composition)

**Environmental Envelope Analysis**:
1. **Elevation Range**: Min, max, optimal for Boletus rex-veris
2. **Canopy Cover Preferences**: Optimal density range
3. **Forest Type Associations**: Which types are preferred/avoided
4. **Negative Habitat Identification**: Where Boletus is never found

#### **Expected Research Questions**

1. **Primary Host Associations**: Which tree species are most strongly correlated with Boletus rex-veris?
2. **Elevation Envelope**: What is the precise elevation range for this species?
3. **Canopy Requirements**: Optimal forest density (not too open, not too closed)?
4. **Negative Indicators**: Which forest types to avoid in surveys?
5. **Succession Preferences**: Young vs. mature forest associations?
6. **Management Impacts**: How do thinning/harvesting affect habitat quality?

#### **Coordinate Uncertainty Handling**

**For Exact Coordinates** (geoprivacy=open):
- Use precise 2km search radius
- Calculate exact elevation from DEM
- High confidence in spatial relationships

**For Obscured Coordinates** (geoprivacy=obscured):
- Expand search radius to 5-10km
- Use elevation bands rather than exact elevation
- Lower weight in statistical analysis
- Reserve for model validation

#### **Output Products**

1. **Species Association Matrix**: Quantified correlations with confidence intervals
2. **Habitat Suitability Maps**: Probability surfaces for Oregon/Washington
3. **Environmental Envelope**: Precise parameter ranges for optimal habitat
4. **Survey Recommendations**: Priority areas for future Boletus searches
5. **Management Guidelines**: Forest practices that enhance/degrade habitat

---

## References
- FIADB User Guide: `data/raw/FIADB User Guide P2_9-1_final.pdf`
- WA FIA Database: `data/raw/trees_SQLite_FIADB_WA.db`
- OR FIA Database: `data/raw/trees_SQLite_FIADB_OR.db`

## Technical Implementation Notes

### Coordinate Calculations
```python
# Convert 120ft spacing to decimal degrees (approximate)
# 120 feet ≈ 36.6 meters
# At 47°N latitude: 1 degree ≈ 111,000m
# 120ft spacing ≈ 0.00033 degrees latitude
# Longitude correction: divide by cos(latitude)

import math
plot_lat = 47.0  # example
plot_lon = -120.0  # example
spacing_degrees = 120 * 0.3048 / 111000  # 120ft to decimal degrees

# Subplot coordinates:
subplot_coords = {
    1: (plot_lat, plot_lon),  # Center
    2: (plot_lat + spacing_degrees, plot_lon),  # North  
    3: (plot_lat + spacing_degrees * math.cos(math.radians(120)), 
        plot_lon + spacing_degrees * math.sin(math.radians(120))),  # SE
    4: (plot_lat + spacing_degrees * math.cos(math.radians(240)),
        plot_lon + spacing_degrees * math.sin(math.radians(240)))   # SW
}
```

### Spatial Accuracy Considerations
- **GPS accuracy**: FIA plot coordinates typically accurate to ±10-30 feet
- **Subplot positioning**: Compass bearings may have ±5° accuracy
- **Recommended buffer**: Include 10-15 foot uncertainty buffer for subplot boundaries
- **Quality control**: Compare calculated subplot positions with any available subplot coordinate data

---
*Last updated: December 2024*
*Focus: Mushroom habitat correlation analysis with Boletus rex-veris priority* 