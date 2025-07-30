# Mushroom Habitat Modeling Project

## Predictive Modeling for Pacific Northwest Mushroom Fruiting

### Project Overview
**Ultimate Goal**: Create predictive models for mushroom fruiting patterns across the Pacific Northwest using forest inventory data, environmental covariates, weather/precip,snow and mushroom observation records.

- Create interesting visuals for users of a map application
- find optimal times and locations for locating mushrooms
- understand more about how weather patterns affect mushroom fruiting
- help people discover new places to look for mushrooms by applying 

**Target Species**:
- 4 Chanterelle species (*Cantharellus* spp.)
- Boletus genus (King boletes, etc.)
- Morel species (*Morchella* spp.) 
- Matsutake complex (*Tricholoma matsutake*, *T. bakamatsutake*, etc.)

---

## 🍄 PRIORITY PROJECT: Boletus rex-veris Habitat Correlation Analysis

### Current Data Status
**iNaturalist Dataset**: `data/raw/iNat/observations-585945.csv`
- **792 total** Boletus rex-veris observations
- **514 with exact coordinates** (65%) - excellent sample size for analysis!
- **277 with obscured coordinates** (35%) - reserve for validation
- **53 in Oregon and Washington** specifically

**FIA Databases Available**:
- `data/raw/trees_SQLite_FIADB_OR.db` - Oregon forest plots
- `data/raw/trees_SQLite_FIADB_WA.db` - Washington forest plots

### Phase 1: Exact Coordinate Analysis (PRIORITY)

#### Step 1: Data Preparation and Filtering
1. **Extract Oregon/Washington Boletus rex-veris observations**
   - Filter for `geoprivacy=open` AND `coordinates_obscured=false` 
   - Focus on WA/OR observations first (53 confirmed)
   - Expand to broader Pacific Northwest region for larger sample

2. **Spatial and Temporal Quality Control**
   - Remove observations with `positional_accuracy > 100m`
   - Filter by observation quality (`quality_grade = research`)
   - Extract elevation from observation coordinates using DEM
   - Seasonal analysis: focus on prime fruiting months (May-July)

#### Step 2: FIA Plot Matching Strategy
1. **Proximity Matching**
   - Find FIA plots within 2km radius of each Boletus observation
   - Weight by distance (closer plots = higher correlation strength)
   - Account for plot coordinate "fuzziness" (±0.5-1.0 mile uncertainty)

2. **Elevation Matching**
   - Match observations to plots within ±200ft elevation range
   - Secondary analysis: ±500ft for broader patterns
   - Document elevation preferences of Boletus rex-veris

3. **Temporal Considerations**
   - Use most recent FIA inventory data (2010-2020 cycles)
   - Account for forest succession since plot measurement

#### Step 3: Forest Composition Analysis

**Primary Variables to Correlate**:

1. **Tree Species Composition** (from TREE table)
   - Species presence/absence within 2km radius
   - Relative basal area by species
   - Above-ground carbon content by species (`CARBON_AG`)
   - Dominant species (>50% basal area)

2. **Forest Structure** (from COND table)
   - Forest type classification (`FORTYPCD`)
   - Canopy cover percentage (`LIVE_CANOPY_CVR_PCT`)
   - Stand age (`STDAGE`)
   - Stand size class (`STDSZCD`)

3. **Habitat Classification** (from REF_HABTYP_DESCRIPTION)
   - Plant association groups (where available)
   - Habitat type correlations

4. **Topographic Variables** (from PLOT/COND tables)
   - Elevation range and preferences
   - Slope and aspect preferences
   - Site productivity class

#### Step 4: Statistical Correlation Analysis

**Positive Correlations** (Species/Conditions favored by Boletus rex-veris):
- Species co-occurrence frequency analysis
- Chi-square tests for species association
- Logistic regression: P(Boletus | Species composition)

**Negative Correlations** (Species/Conditions avoided by Boletus rex-veris):
- Identify forest types where Boletus is never/rarely found
- Exclusion analysis: conditions that prevent Boletus occurrence

**Environmental Envelope Analysis**:
- Elevation range (min, max, optimal)
- Canopy cover preferences (too dense vs. too open)
- Slope and aspect preferences

#### Step 5: Habitat Suitability Modeling

1. **Training Data Preparation**
   - Boletus presence locations (514 exact coordinates)
   - Background/absence locations (systematic sampling)
   - Environmental variables for each location

2. **Model Development**
   - Random Forest classification
   - MaxEnt species distribution modeling
   - Ensemble model combining approaches

3. **Model Validation**
   - Cross-validation with held-out Boletus observations
   - Validation against obscured coordinate observations
   - Expert review by mycologists

### Phase 2: Expanded Analysis (If Phase 1 Successful)

#### Including Obscured Coordinates
If exact coordinate analysis shows strong patterns:
1. **Fuzzy Matching Approach**
   - Use obscured coordinates with larger search radius (5-10km)
   - Weight by uncertainty area
   - Ensemble predictions across possible locations

2. **Habitat Envelope Validation**
   - Test whether obscured observations fall within predicted suitable habitat
   - Refine model parameters based on expanded dataset

#### Comparison with Other Bolete Species
- *Boletus edulis* (if available in dataset)
- Habitat niche differentiation analysis
- Competitive exclusion patterns

### Expected Outcomes

#### Quantitative Results
1. **Species Association Matrix**
   - Tree species positively correlated with Boletus rex-veris
   - Tree species negatively correlated or absent from Boletus habitat
   - Quantified association strengths (correlation coefficients)

2. **Environmental Envelope**
   - Elevation range: Expected 2,000-6,000 feet based on literature
   - Canopy cover: Expected 40-80% (not too dense, not too open)
   - Forest type preferences: Likely Douglas-fir, True fir associations

3. **Habitat Suitability Map**
   - Probability surface across Oregon and Washington
   - Validation statistics (AUC, sensitivity, specificity)
   - Confidence intervals for predictions

#### Biological Insights
1. **Mycorrhizal Associations**
   - Primary host tree species for Boletus rex-veris
   - Secondary host associations
   - Forest successional stage preferences

2. **Habitat Constraints**
   - Minimum canopy cover thresholds
   - Maximum elevation limits
   - Slope and aspect limitations

3. **Management Implications**
   - Forest practices that enhance/degrade Boletus habitat
   - Indicator species for suitable habitat
   - Conservation priority areas

### Questions This Analysis Will Answer

#### Species-Specific Questions
1. **What tree species are most strongly associated with Boletus rex-veris?**
   - Expected: Douglas-fir, True firs, possibly Western Hemlock
   - Quantified through co-occurrence analysis

2. **What forest types should be avoided in Boletus surveys?**
   - Identify "negative habitat" - areas where Boletus is never found
   - Pure hardwood stands, very young forests, clearcuts

3. **What is the optimal elevation range for Boletus rex-veris?**
   - Define precise elevation envelope from data
   - Compare to literature values (2,000-6,000 ft)

4. **How does canopy cover affect Boletus occurrence?**
   - Threshold analysis for optimal canopy density
   - Balance between too open (no mycorrhizal partners) and too dense (insufficient light)

#### Ecological Questions
5. **Do Boletus rex-veris and other mushroom species compete or co-occur?**
   - Species interaction patterns
   - Niche partitioning by elevation, forest type, etc.

6. **How does forest age/succession affect Boletus habitat suitability?**
   - Stand age correlations
   - Successional stage preferences

7. **What role does forest management play in Boletus habitat?**
   - Thinning effects, harvest rotation impacts
   - Natural vs. managed forest preferences

#### Predictive Questions
8. **Can we predict new Boletus locations based on forest composition?**
   - Model accuracy assessment
   - Validation against independent observations

9. **What areas have suitable habitat but no current observations?**
   - Gap analysis for survey prioritization
   - New location recommendations for foragers

10. **How might climate change affect Boletus habitat in the Pacific Northwest?**
    - Elevation range shifts
    - Host tree distribution changes

---

### Multi-Scale Data Framework

#### 0. Other Ideas
Applying Sample Plot Data - more broadly
- could I catergorize some items visually

#### 1. High-Resolution Point Data (FIA Subplots)
**Source**: FIADB subplot-level data
**Spatial Scale**: Exact coordinates, ~0.1 acre subplots
**Temporal**: Specific survey years (5-10 year cycles)
**Data Type**: Precise species composition, DBH, carbon, forest structure

**Key Variables**:
- Tree species composition and dominance
- Stand age and structure
- Canopy closure
- Understory vegetation
- Dead wood/CWD volume
- Soil characteristics (where available)

#### 2. Moderate-Resolution Point Data (FIA Plots) 
**Source**: FIADB plot-level aggregated data
**Spatial Scale**: ~1 acre plots, fuzzy coordinates and sometimes swapped, possibly only swapped when private land.
- I could use fuzzy coordinates with the elevation data to narrow it down
**Temporal**: Survey year averages
**Data Type**: Plot-level forest type classification

**Key Variables**:
- Dominant forest type (Douglas-fir, Western Hemlock, Mixed Conifer, etc.)
- Species diversity indices
- Basal area by species
- Stand density
- Forest productivity metrics

#### 3. Landscape-Scale Predictions (1-10km grid)
**Source**: Environmental covariates + model predictions
**Spatial Scale**: 1km to 10km grid cells
**Temporal**: Multi-year averages or seasonal
**Data Type**: Predicted forest type probabilities

**Environmental Predictors**:

##### Climate/Weather
- Temperature (mean, min, max, seasonality)
- Precipitation (annual, seasonal, timing)
- Humidity/moisture indices
- Growing season length
- Frost dates

##### Topography
- Elevation
- Slope and aspect
- Topographic wetness index
- Heat load index
- Terrain ruggedness
- Distance to ridges/valleys

##### Soils & Geology
- Soil types (SSURGO database)
- Soil pH and nutrient status
- Soil drainage class
- Organic matter content
- Parent material/bedrock geology
- Mineral composition

##### Remote Sensing
- NDVI, EVI (vegetation health)
- Landsat forest type classifications
- Canopy cover percentage (NLCD)
- Leaf Area Index (LAI)
- Surface temperature
- Moisture indices

##### Hydrology
- Distance to streams/rivers
- Watershed characteristics
- Stream flow regimes
- Groundwater proximity

### Forest Type Classification System

#### FIA Forest Types (Primary)
- Douglas-fir (*Pseudotsuga menziesii*)
- Western Hemlock (*Tsuga heterophylla*)
- Western Red Cedar (*Thuja plicata*)
- Noble Fir (*Abies procera*)
- Pacific Silver Fir (*Abies amabilis*)
- Sitka Spruce (*Picea sitchensis*)
- Mixed Conifer
- Hardwood (Alder, Maple, etc.)

#### Ecological Classifications (Secondary)
- Plant Association Groups (PAG)
- Habitat Types
- Potential Vegetation Groups
- Ecoregion classifications

### Mushroom Observation Data

#### GBIF Data Integration
**Source**: Global Biodiversity Information Facility
**Target Species**: iNaturalist, research collections, citizen science
**Data Processing**:
- Spatial accuracy filtering (≤1km uncertainty)
- Temporal filtering (recent decades)
- Taxonomic validation
- Seasonal fruiting pattern analysis

#### Data Structure
```
mushroom_observations:
- species_name
- latitude/longitude
- observation_date
- habitat_description
- associated_tree_species
- elevation
- data_source
- confidence_level
```

### Modeling Framework

#### Stage 1: Forest Type Prediction
**Input**: FIA plots + environmental covariates
**Output**: Forest type probability maps
**Methods**: Random Forest, Gradient Boosting, Neural Networks
**Validation**: Cross-validation with held-out FIA plots

#### Stage 2: Mushroom Habitat Suitability
**Input**: Forest type predictions + environmental variables + GBIF data
**Output**: Species-specific habitat suitability maps
**Methods**: MaxEnt, Random Forest, Ensemble modeling
**Validation**: Cross-validation with GBIF observations

#### Stage 3: Seasonal/Temporal Modeling
**Input**: Climate data + phenology models
**Output**: Fruiting timing predictions
**Methods**: Time series analysis, climate envelope modeling

### Data Sources to Acquire

#### Soil Data
- SSURGO (Soil Survey Geographic Database)
- STATSGO (State Soil Geographic Database) 
- POLARIS soil properties (pH, nutrients, texture)
- SoilGrids global soil information

#### Climate Data
- PRISM (Parameter-elevation Regressions on Independent Slopes Model)
- Daymet daily climate data
- WorldClim bioclimatic variables
- MODIS land surface temperature

#### Remote Sensing
- Landsat time series (30m resolution)
- Sentinel-2 (10m resolution)
- MODIS vegetation indices
- Lidar-derived forest structure (where available)

#### Hydrological
- National Hydrography Dataset (NHD)
- USGS stream gauge data
- Watershed boundary datasets

### Implementation Phases

#### Phase 1: Data Assembly & Forest Type Modeling
1. Extract comprehensive FIA data (all tree species, not just Douglas-fir)
2. Assemble environmental covariate layers
3. Develop forest type classification models
4. Create forest type probability maps

#### Phase 2: Mushroom Data Integration
1. Query GBIF for target mushroom species in Pacific Northwest
2. Clean and validate mushroom observation data
3. Spatial and temporal analysis of fruiting patterns
4. Habitat association analysis

#### Phase 3: Habitat Suitability Modeling
1. Develop species-specific habitat models
2. Integrate forest type predictions with environmental variables
3. Create mushroom habitat suitability maps
4. Validate against independent observation data

#### Phase 4: Temporal/Seasonal Modeling
1. Analyze seasonal fruiting patterns
2. Develop climate-based fruiting timing models
3. Create seasonal prediction maps
4. Integrate with habitat suitability for comprehensive forecasting

### Technical Considerations

#### Spatial Resolution Strategy
- **Fine scale (30m-100m)**: Detailed habitat mapping for specific sites
- **Medium scale (250m-1km)**: Regional habitat assessment
- **Coarse scale (1km-10km)**: Landscape-level patterns and planning

#### Temporal Resolution
- **Annual**: Forest type and general habitat suitability
- **Seasonal**: Fruiting timing predictions
- **Event-based**: Post-precipitation fruiting forecasts

#### Uncertainty Quantification
- Model ensemble approaches
- Bootstrap confidence intervals
- Spatial uncertainty propagation
- Validation uncertainty assessment

### Success Metrics
1. **Forest Type Accuracy**: >80% classification accuracy on held-out FIA plots
2. **Habitat Model Performance**: AUC >0.8 for mushroom presence/absence
3. **Spatial Validation**: Correspondence with known mushroom hotspots
4. **Temporal Validation**: Accurate fruiting timing predictions
5. **Expert Validation**: Review by mycologists and forest ecologists 