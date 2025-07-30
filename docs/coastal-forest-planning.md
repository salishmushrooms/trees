# Coastal Forest Bioregion Planning Document

## Project Objective
Create a coastal forest bioregion that:
- Extends to the actual coastline/water boundary
- Follows natural topographic boundaries (valleys, ridges)
- Uses elevation constraints to exclude high-elevation areas unsuitable for shore pine
- Has irregular, natural boundaries rather than geometric buffers
- Can be layered under water in Mapbox for natural masking

## Available Data Resources

### 1. Tree Cover Data
**Primary Source**: `/Users/JJC/trees/outputs/mapbox_masks/pnw_tree_cover_90m_mapbox.tif`
- **Resolution**: 90m 
- **Coverage**: Pacific Northwest region
- **Values**: 0-100% tree canopy cover
- **Use Case**: Exclude non-forested areas, define forest boundaries

**Alternative Sources**:
- `pnw_tree_cover_30m_full.tif` (30m resolution, higher detail)
- `pnw_tree_cover_60m_half.tif` (60m resolution)
- `pnw_tree_cover_120m_quarter.tif` (120m resolution, faster processing)

### 2. Digital Elevation Model (DEM)
**Current Sources** (need to verify):
- Cached DEM data in `data/raw/` directory
- NLCD elevation data (referenced in existing scripts)
- Possible USGS 30m DEM

**Elevation Constraints**:
- Shore pine (Pinus contorta var. contorta): 0-500ft elevation
- Coastal transition zone: 500-1000ft elevation
- Maximum coastal forest: <1500ft elevation

### 3. Cached Elevation Masks
**Available**: `cache/species_masks/pnw_elevation_mask_*`
- Various elevation ranges (0-4500ft, 0-5000ft, etc.)
- Pre-processed GeoJSON polygons
- Ready for spatial operations

### 4. Tree Species Database Analysis
**CONFIRMED**: REF_SPECIES table has subspecies fields, but they're NOT populated for Pinus contorta!
- **Primary Table**: `TREE` table with `SPCD` codes
- **Species Reference**: `REF_SPECIES` table
- **Subspecies Fields**: VARIETY and SUBSPECIES exist but are NULL for SPCD 108
- **Key Finding**: **FIA database does NOT distinguish shore pine from lodgepole pine**
- **Implication**: Must use geographic location + elevation to identify shore pine habitat
- **Shore Pine Strategy**: Use coastal proximity (<50km) + low elevation (<500ft) as proxies

### 5. Existing Coastal Species Data
**Processed Data**: `outputs/coastal_analysis/coastal_forest_plots.geojson`
- 2,234 coastal plots at <1000ft elevation
- 418 plots with Sitka spruce
- 52 plots with lodgepole pine (includes shore pine)
- Geographic range: -124.70� to -123.00� longitude

## Technical Approach

### Phase 1: Data Assessment
1. **Subspecies Investigation**:
   ```sql
   SELECT DISTINCT SPCD, COMMON_NAME, SCIENTIFIC_NAME 
   FROM REF_SPECIES 
   WHERE SPCD = 108 OR SCIENTIFIC_NAME LIKE '%contorta%'
   ```

2. **DEM Data Inventory**:
   - Locate and assess available DEM sources
   - Determine resolution and coverage
   - Test elevation extraction at known coastal plot locations

3. **Tree Cover Analysis**:
   - Load 90m Mapbox tree cover raster
   - Analyze coverage along known coastline
   - Determine optimal threshold for coastal forests

### Phase 2: Coastal Boundary Definition
1. **Water Boundary Approach**:
   - Since Mapbox water layer will mask the coastline
   - Create bioregion that extends into water areas
   - Let water layer provide natural coastal boundary

2. **Elevation-Based Inland Boundary**:
   - Use DEM to create elevation contours
   - Apply species-specific elevation limits
   - Create smooth transitions between elevation zones

3. **Tree Cover Integration**:
   - Apply minimum tree cover threshold (5-10%)
   - Remove non-forested areas from bioregion
   - Maintain connectivity of forest patches

### Phase 3: Bioregion Refinement
1. **Topographic Following**:
   - Use DEM to identify ridges and valleys
   - Extend bioregion up coastal valleys
   - Respect natural watershed boundaries

2. **Species Validation**:
   - Overlay known Sitka spruce locations
   - Overlay known shore pine locations  
   - Adjust boundaries to include species concentrations

3. **Connectivity Analysis**:
   - Ensure bioregion forms connected forest corridors
   - Bridge small gaps in forest cover
   - Remove isolated fragments below minimum size

### Phase 4: Output Optimization
1. **Mapbox Integration**:
   - Optimize geometry for web mapping
   - Reduce polygon complexity where appropriate
   - Test layering with water boundaries

2. **Multi-Resolution Outputs**:
   - High-detail version for analysis
   - Simplified version for web mapping
   - Point validation datasets

## Implementation Strategy

### Method A: Raster-Based Approach (Recommended)
```python
# Workflow
1. Load tree cover raster (90m)
2. Load DEM raster
3. Create elevation masks (<500ft, <1000ft, <1500ft)
4. Apply tree cover threshold (>5%)
5. Combine constraints with logical AND
6. Vectorize resulting raster to polygons
7. Apply morphological operations for smoothing
8. Validate against species occurrence data
```

**Advantages**:
- Follows natural boundaries precisely
- Integrates multiple data sources seamlessly
- Produces irregular, realistic boundaries
- Computationally efficient

### Method B: Species-Driven Buffer Approach
```python
# Workflow  
1. Create buffers around coastal species plots
2. Merge overlapping buffers
3. Clip to elevation constraints
4. Clip to tree cover constraints
5. Smooth boundaries
```

**Advantages**:
- Directly based on species occurrences
- Simpler implementation
- Easy to validate

## Data Processing Requirements

### Computing Considerations
- **Memory**: Tree cover raster (~90m) may be large
- **Processing Time**: Raster operations can be slow
- **Storage**: Multiple intermediate outputs needed

### Quality Control Checks
1. **Elevation Validation**: Check that bioregion elevations match species data
2. **Tree Cover Validation**: Verify adequate forest cover in bioregion
3. **Species Coverage**: Ensure bioregion captures majority of coastal species plots
4. **Boundary Quality**: Check for unrealistic geometric artifacts

## Expected Outputs

### Primary Deliverable
- `coastal_forest_bioregion_final.geojson`: Detailed coastal forest boundary
- Irregular boundary following coastline and topography
- Area estimate: 8,000-15,000 km� (refined from current 14,120 km�)

### Supporting Deliverables
- `coastal_forest_simplified.geojson`: Web-optimized version
- `coastal_forest_validation.geojson`: Species plot coverage analysis
- `coastal_forest_metadata.json`: Processing parameters and statistics

### Integration Outputs
- `coastal_forest_90m.tif`: Raster mask for further analysis
- Mapbox-ready styling specifications
- Documentation for water layer integration

## Data Investigation Results

### Tree Cover Raster Analysis ✅
**CONFIRMED**: Tree cover data is excellent for coastal areas!
- **Best raster**: `pnw_tree_cover_30m_full.tif` (30m resolution, most detail)
- **Good alternative**: `pnw_tree_cover_90m_mapbox.tif` (90m, faster processing)
- **Coastal coverage**: All areas show good tree cover data (32-64% mean coverage)
- **Olympic Peninsula**: 58.7% mean tree cover
- **WA/OR Coast**: 43-64% mean tree cover
- **Port Angeles**: 32% mean tree cover (lower due to urban areas)

### Subspecies Investigation ✅  
**CONFIRMED**: No subspecies differentiation, but we have coastal Pinus contorta data!
- **Total coastal Pinus contorta**: 1,528 trees in 101 plots (WA: 502 trees/45 plots, OR: 1,026 trees/56 plots)
- **Elevation range**: 9-500 ft (perfect for shore pine identification)
- **Average elevation**: 132-216 ft
- **Geographic distribution**: All west of -123.0° longitude

### Elevation Data Status ✅
**CONFIRMED**: Excellent elevation data available!
- **DEM Raster**: `outputs/mapbox_masks/pnw_elevation_120m_mapbox.tif` (120m resolution, 208MB)
- **Elevation Masks**: Multiple pre-processed masks available
- **Key masks for coastal**: `pnw_elevation_mask_0_5000ft.geojson`, `pnw_elevation_mask_0_4500ft.geojson`
- **Shore pine constraint**: Can use 0-500ft elevation mask for precise shore pine habitat

## Terminal Commands for Data Investigation

### 1. Database Schema Investigation
```bash
# Check REF_SPECIES table structure and Pinus contorta entries
sqlite3 data/raw/trees_SQLite_FIADB_WA.db << 'EOF'
.schema REF_SPECIES
.mode column
.headers on
SELECT SPCD, COMMON_NAME, GENUS, SPECIES, VARIETY, SUBSPECIES 
FROM REF_SPECIES 
WHERE SPCD = 108 OR GENUS LIKE '%Pinus%' AND SPECIES LIKE '%contorta%'
LIMIT 5;
.quit
EOF
```

### 2. Tree Cover Raster Analysis
```bash
# Check raster properties safely
python3 -c "
import rasterio
import numpy as np

# Check 90m tree cover raster
with rasterio.open('outputs/mapbox_masks/pnw_tree_cover_90m_mapbox.tif') as src:
    print('Tree Cover Raster Info:')
    print(f'Shape: {src.shape}')
    print(f'CRS: {src.crs}')
    print(f'Bounds: {src.bounds}')
    
    # Read coastal area sample (WA coast)
    window = rasterio.windows.from_bounds(-125, 46, -123, 48, src.transform)
    data = src.read(1, window=window)
    print(f'Coastal sample - Min: {data.min()}, Max: {data.max()}, Mean: {data.mean():.1f}')
"
```

### 3. DEM Data Discovery
```bash
# Find elevation/DEM files
find data/raw -name "*elevation*" -o -name "*dem*" -o -name "*elev*" -type f
find . -name "*.tif" | grep -i elev
ls -la data/raw/nlcd_tcc_*
```

### 4. Coastal Species Subspecies Check
```bash
# Check if any coastal plots have subspecies data
sqlite3 data/raw/trees_SQLite_FIADB_WA.db << 'EOF'
SELECT COUNT(*) as total_pinus_contorta,
       COUNT(CASE WHEN r.VARIETY IS NOT NULL THEN 1 END) as with_variety,
       COUNT(CASE WHEN r.SUBSPECIES IS NOT NULL THEN 1 END) as with_subspecies
FROM TREE t
JOIN REF_SPECIES r ON t.SPCD = r.SPCD
WHERE t.SPCD = 108;
.quit
EOF
```

## Implementation Plan

### READY TO IMPLEMENT! 🚀

All data sources confirmed and available:

#### Available Data Assets:
1. **Tree Cover**: `pnw_tree_cover_30m_full.tif` (30m resolution)
2. **Elevation**: `pnw_elevation_120m_mapbox.tif` (120m resolution) 
3. **Elevation Masks**: Pre-processed 0-500ft mask for shore pine
4. **Species Validation**: 101 coastal Pinus contorta plots (9-500ft elevation)

#### Implementation Strategy:
```python
# Raster-based coastal forest creation workflow:
1. Load tree cover raster (30m) - Apply >10% threshold for coastal forests
2. Load elevation raster (120m) - Apply <500ft constraint for shore pine habitat  
3. Create coastal proximity buffer - Extend from known species plots
4. Combine all constraints with logical AND operations
5. Vectorize result to create coastal forest polygon
6. Extend polygon to water boundary for Mapbox integration
7. Validate against 101 known shore pine plot locations
```

#### Expected Output:
- Irregular coastal forest boundary following natural topography
- Extends to water's edge for Mapbox water layer masking
- Captures shore pine habitat (0-500ft elevation)
- Includes Sitka spruce coastal forests
- Area estimate: 8,000-12,000 km² (refined from our previous attempts)

### Next Steps Priority

1. **HIGH PRIORITY**: Create raster-based coastal forest bioregion script
2. **HIGH PRIORITY**: Implement tree cover + elevation constraints  
3. **MEDIUM PRIORITY**: Validate against known coastal species plots
4. **MEDIUM PRIORITY**: Optimize geometry for Mapbox integration
5. **LOW PRIORITY**: Create multiple resolution outputs for different uses

## Research Questions

1. **Subspecies Data**: Does FIA database distinguish Pinus contorta varieties?
2. **Elevation Precision**: What DEM resolution is optimal for coastal analysis?
3. **Tree Cover Threshold**: What minimum canopy cover defines coastal forests?
4. **Boundary Complexity**: How much detail is needed for web mapping?
5. **Seasonal Variation**: Do we need to account for seasonal forest changes?

---

*This document will be updated as data assessment and prototyping proceed.*