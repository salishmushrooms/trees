# FIA Database Structure Reference Guide

## Overview

This document provides a reference for navigating and querying the FIA (Forest Inventory and Analysis) database for forest and mushroom habitat modeling. It focuses on key tables, fields, and relationships most relevant for ecological analysis.

Official Documentation Here
data/raw/FIADB Trees Guide.pdf

**Important Note on Plot Locations**: FIA plot coordinates are **NOT EXACT** due to privacy protection. Coordinates are "fuzzed" (moved up to 0.5-1.0 mile) and up to 20% of private plots are "swapped" with similar plots in the same county to protect landowner privacy.

## Core Data Tables

### 1. PLOT Table - Primary Plot Information
**Purpose**: Main table containing plot-level metadata and location information

**Key Fields for Analysis**:
- `CN`: Unique plot identifier (VARCHAR2(34)) - **Primary key for all joins**
- `STATECD`: State code (NUMBER(4)) - 53 for Washington
- `UNITCD`: Survey unit code (NUMBER(2))
- `COUNTYCD`: County code (NUMBER(3))  
- `PLOT`: Plot number (NUMBER(5))
- `INVYR`: Inventory year (NUMBER(4)) - Survey year
- `LAT`: Latitude (NUMBER(8,6)) - **Fuzzed, not exact**
- `LON`: Longitude (NUMBER(9,6)) - **Fuzzed, not exact**
- `ELEV`: Elevation in feet (NUMBER(5))
- `PLOT_STATUS_CD`: Plot status (NUMBER(1)) - 1=Forested, 2=Nonforested, 3=Noncensus water
- `MEASYEAR`, `MEASMON`, `MEASDAY`: Measurement timing
- `CYCLE`, `SUBCYCLE`: Inventory cycle information

**Usage Notes**:
- Use `PLT_CN = PLOT.CN` for all joins to other tables
- Filter `PLOT_STATUS_CD = 1` for forested plots only
- Coordinates suitable for regional analysis but not precise ground-truthing

### 2. SUBPLOT Table - Individual Subplot Information  
**Purpose**: Each plot has 4 subplots arranged in 120-foot spacing pattern

**Key Fields**:
- `CN`: Unique subplot identifier (VARCHAR2(34))
- `PLT_CN`: Links to PLOT.CN (VARCHAR2(34)) - **Foreign key**
- `SUBP`: Subplot number (NUMBER) - 1,2,3,4 
- `SUBP_STATUS_CD`: Subplot status (NUMBER)
- `MICRCOND`, `SUBPCOND`, `MACRCOND`: Condition assignments

**Subplot Layout**:
```
    Subplot 2 (North, 120ft from center)
         |
Subplot 3 --- Subplot 1 (Center) --- [120 ft spacing]  
(SE, 120°)    (Plot center)
         |
    Subplot 4 (South, 240°)
```

**Measurement Areas**:
- **Microplot**: 6.8 ft radius - trees <5.0" DBH
- **Subplot**: 24.0 ft radius - trees ≥5.0" DBH (**primary measurement area**)
- **Macroplot**: 58.9 ft radius - select measurements

### 3. COND (Condition) Table - Forest Conditions and Habitat Classification
**Purpose**: Describes forest conditions within plots, including forest types and habitat classifications

**Key Fields for Habitat Analysis**:
- `CN`: Unique condition identifier (VARCHAR2(34))
- `PLT_CN`: Links to PLOT.CN (VARCHAR2(34)) - **Foreign key**
- `CONDID`: Condition ID within plot (NUMBER) - Usually 1 for uniform plots

**Forest Classification**:
- `FORTYPCD`: **Forest type code** (NUMBER) - Links to REF_FOREST_TYPE
- `FLDTYPCD`: Field-determined forest type (NUMBER)
- `ADFORCD`: Administrative forest code (NUMBER)

**Habitat Classification** ⭐:
- `HABTYPCD1`: **Primary habitat type code** (VARCHAR2(10)) - Links to REF_HABTYP_DESCRIPTION
- `HABTYPCD1_PUB_CD`: Publication code for habitat type (VARCHAR2(10))
- `HABTYPCD2`: Secondary habitat type (VARCHAR2(10))

**Stand Characteristics**:
- `STDAGE`: Stand age (NUMBER)
- `STDSZCD`: Stand size class (NUMBER)
- `SITECLCD`: Site class (NUMBER)
- `SLOPE`: Slope percentage (NUMBER) 
- `ASPECT`: Aspect in degrees (NUMBER)
- `LIVE_CANOPY_CVR_PCT`: Live canopy cover percentage (NUMBER)

**Usage**: 55% of Washington conditions have habitat type data (16,113/29,073 records)

### 4. TREE Table - Individual Tree Measurements
**Purpose**: Individual tree data within subplots

**Key Fields**:
- `CN`: Unique tree identifier (VARCHAR2(34))
- `PLT_CN`: Links to PLOT.CN (VARCHAR2(34)) - **Foreign key**
- `SUBP`: Subplot number (NUMBER) - Links to SUBPLOT.SUBP
- `SPCD`: **Species code** (NUMBER) - Links to REF_SPECIES
- `DIA`: Diameter at breast height (NUMBER)
- `HT`: Tree height (NUMBER)
- `CARBON_AG`: Above-ground carbon (NUMBER)
- `STATUSCD`: Tree status - 1=Live, 2=Dead

### 5. SUBP_COND Table - Subplot-Condition Links
**Purpose**: Links subplots to their condition classifications

**Key Fields**:
- `PLT_CN`: Plot identifier (VARCHAR2(34))
- `SUBP`: Subplot number (NUMBER)
- `CONDID`: Condition ID (NUMBER)
- `SUBPCOND_PROP`: Proportion of subplot in this condition (NUMBER)

## Reference Tables (Lookup/Classification)

### REF_FOREST_TYPE - Forest Type Classifications
**Purpose**: Standardized forest type descriptions

**Key Pacific Northwest Forest Types** (from Washington data):
```
201 - Douglas-fir (6,393 plots - 44% of WA forests)
301 - Western hemlock (1,672 plots - 11%)  
221 - Ponderosa pine (1,300 plots - 9%)
264 - Pacific silver fir (933 plots - 6%)
304 - Western redcedar (476 plots - 3%)
270 - Mountain hemlock (323 plots - 2%)
305 - Sitka spruce (for coastal areas)
```

**Fields**:
- `VALUE`: Forest type code (NUMBER) - Links to COND.FORTYPCD
- `MEANING`: Forest type description (VARCHAR2(80))
- `TYPGRPCD`: Forest type group code (NUMBER)

### REF_HABTYP_DESCRIPTION - Habitat Type Classifications ⭐
**Purpose**: Detailed plant association and habitat classifications (8,712 records)

**Common Washington Habitat Types**:
```
TSHE/GASH/POMU - western hemlock/salal/western swordfern (261 plots)
TSHE/POMU-OXOR - western hemlock/western swordfern-redwood-sorrel (239 plots)
PSME/PHMA5 - Douglas-fir/mallow ninebark (199 plots)
THPL/CLUN2 - western red cedar/bride's bonnet (124 plots)
```

**Fields**:
- `HABTYPCD`: Habitat type code (VARCHAR2(10)) - Links to COND.HABTYPCD1
- `PUB_CD`: Publication code (VARCHAR2(10)) - Links to COND.HABTYPCD1_PUB_CD  
- `SCIENTIFIC_NAME`: Scientific plant association name (VARCHAR2(115))
- `COMMON_NAME`: Common habitat description (VARCHAR2(255))

### REF_SPECIES - Tree Species Information
**Purpose**: Tree species codes and characteristics

**Key Pacific Northwest Species**:
```
202 - Douglas-fir (Pseudotsuga menziesii)
263 - Western hemlock (Tsuga heterophylla) 
242 - Western red cedar (Thuja plicata)
017 - Pacific silver fir (Abies amabilis)
011 - Noble fir (Abies procera)
098 - Sitka spruce (Picea sitchensis)
```

**Fields**:
- `SPCD`: Species code (NUMBER) - Links to TREE.SPCD
- `COMMON_NAME`: Common name (VARCHAR2(100))
- `SCIENTIFIC_NAME`: Scientific name (VARCHAR2(4000))
- `SFTWD_HRDWD`: Softwood/Hardwood classification (VARCHAR2(1))

## Data Relationships and Query Patterns

### Basic Plot → Condition → Habitat Lookup
```sql
SELECT p.PLOT, p.LAT, p.LON, p.ELEV,
       f.MEANING as FOREST_TYPE,
       h.SCIENTIFIC_NAME as HABITAT_TYPE
FROM PLOT p
JOIN COND c ON p.CN = c.PLT_CN AND c.CONDID = 1
JOIN REF_FOREST_TYPE f ON c.FORTYPCD = f.VALUE
LEFT JOIN REF_HABTYP_DESCRIPTION h ON c.HABTYPCD1 = h.HABTYPCD 
    AND c.HABTYPCD1_PUB_CD = h.PUB_CD
WHERE p.PLOT_STATUS_CD = 1;
```

### Subplot-Level Tree Composition
```sql
SELECT p.PLOT, s.SUBP, sp.COMMON_NAME as SPECIES,
       COUNT(*) as TREE_COUNT,
       AVG(t.DIA) as AVG_DBH,
       SUM(t.CARBON_AG) as TOTAL_CARBON
FROM PLOT p
JOIN SUBPLOT s ON p.CN = s.PLT_CN
JOIN TREE t ON p.CN = t.PLT_CN AND s.SUBP = t.SUBP
JOIN REF_SPECIES sp ON t.SPCD = sp.SPCD
WHERE t.STATUSCD = 1  -- Live trees only
GROUP BY p.CN, s.SUBP, t.SPCD;
```

### Habitat-Forest Type Cross-Classification
```sql
SELECT f.MEANING as FOREST_TYPE,
       h.SCIENTIFIC_NAME as HABITAT_TYPE,
       COUNT(*) as PLOT_COUNT
FROM COND c
JOIN REF_FOREST_TYPE f ON c.FORTYPCD = f.VALUE
JOIN REF_HABTYP_DESCRIPTION h ON c.HABTYPCD1 = h.HABTYPCD 
    AND c.HABTYPCD1_PUB_CD = h.PUB_CD
GROUP BY c.FORTYPCD, c.HABTYPCD1
ORDER BY PLOT_COUNT DESC;
```

## Forest Type Applicability for Broader Analysis

**FORTYPCD Analysis Results**:

✅ **FORTYPCD is HIGHLY APPLICABLE** for regional forest classification because:

1. **Standardized Classification**: Forest types are nationally standardized (201=Douglas-fir everywhere)
2. **Ecologically Meaningful**: Based on dominant tree species composition
3. **Regionally Consistent**: Same codes apply across similar ecosystems
4. **Well-Documented**: Extensive reference documentation available

**Geographic Transferability**:
- **Within Region**: Forest types transfer well within Pacific Northwest (WA → OR → ID)
- **Similar Climates**: Douglas-fir type (201) applies wherever Douglas-fir dominates
- **Elevation Zones**: Mountain hemlock (270), subalpine fir (268) define elevation belts
- **Habitat Associations**: Can predict mushroom habitat based on forest type patterns

**Limitations**:
- Forest types are **generalized** - don't capture local microclimates
- **Age/disturbance** not captured in forest type alone
- **Understory variation** better captured by habitat types (HABTYPCD1)

## Navigation Guide for Additional Information

**For detailed field definitions**: See `data/raw/FIADB Trees Guide.pdf`
- Chapter 2: Individual table descriptions
- Appendix F: Code definitions and valid values
- Section 1.2.3: Plot design and measurement protocols

**Key sections to reference**:
- **Plot Design**: Section 1.2.3 (4-subplot arrangement, measurement radii)
- **Forest Types**: Appendix F (complete forest type code list)
- **Species Codes**: Appendix F (complete species list with scientific names)
- **Habitat Types**: Regional habitat type publications (referenced in PUB_CD)

## Recommended Analysis Workflow

1. **Start with PLOT table** - filter for forested plots, appropriate years
2. **Join to COND** - get forest types and habitat classifications  
3. **Add SUBPLOT/TREE data** - for detailed species composition
4. **Use REF tables** - for human-readable descriptions
5. **Spatial analysis** - remember coordinates are approximate (±0.5-1.0 mile)

**For Mushroom Habitat Modeling**:
- **Plot level**: Forest type + environmental variables → broad habitat suitability
- **Subplot level**: Detailed species composition → precise training data
- **Habitat types**: Plant associations → mushroom microhabitat preferences 