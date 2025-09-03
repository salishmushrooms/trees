import sqlite3
import json
import pandas as pd
from shapely.geometry import Point
import geopandas as gpd

# Database paths for WA, OR, ID
databases = {
    'WA': '/Users/JJC/trees/data/raw/trees_SQLite_FIADB_WA.db',
    'OR': '/Users/JJC/trees/data/raw/trees_SQLite_FIADB_OR.db',
    'ID': '/Users/JJC/trees/data/raw/trees_SQLite_FIADB_ID.db'
}

# Query to extract comprehensive forest data for mushroom habitat modeling
# Uses most recent survey data only (latest INVYR) for each plot location
query = """
WITH latest_plots AS (
    SELECT p.CN as PLOT_CN,
           p.STATECD,
           p.UNITCD, 
           p.COUNTYCD,
           p.PLOT,
           p.LAT,
           p.LON,
           p.INVYR,
           p.ELEV,
           p.PLOT_STATUS_CD,
           ROW_NUMBER() OVER (
               PARTITION BY p.STATECD, p.UNITCD, p.COUNTYCD, p.PLOT 
               ORDER BY p.INVYR DESC
           ) as rn
    FROM PLOT p
    WHERE p.LAT IS NOT NULL 
        AND p.LON IS NOT NULL
        AND p.PLOT_STATUS_CD = 1  -- Forested plots only
)
SELECT 
    lp.PLOT_CN,
    lp.STATECD,
    lp.UNITCD, 
    lp.COUNTYCD,
    lp.PLOT,
    lp.LAT,
    lp.LON,
    lp.INVYR,
    lp.ELEV,
    
    -- Plot-level forest type information
    c.FORTYPCD,
    c.FLDTYPCD,
    c.STDAGE,
    c.ADFORCD,
    c.SLOPE,
    c.ASPECT,
    c.LIVE_CANOPY_CVR_PCT as CANOPY_COVER,
    
    -- Tree species composition (aggregated by plot)
    COUNT(t.CN) as TOTAL_TREES,
    
    -- Douglas-fir metrics (Pseudotsuga menziesii) - SPCD 202
    SUM(CASE WHEN t.SPCD = 202 THEN 1 ELSE 0 END) as DOUGLAS_FIR_COUNT,
    ROUND(SUM(CASE WHEN t.SPCD = 202 THEN t.CARBON_AG ELSE 0 END), 0) as DOUGLAS_FIR_CARBON,
    ROUND(SUM(CASE WHEN t.SPCD = 202 THEN t.CARBON_BG ELSE 0 END), 0) as DOUGLAS_FIR_CARBON_BG,
    ROUND(SUM(CASE WHEN t.SPCD = 202 THEN t.DIA * t.DIA ELSE 0 END)) as DOUGLAS_FIR_BA,
    ROUND(SUM(CASE WHEN t.SPCD = 202 THEN POWER((0.3 * t.DIA + 1.2), 2) * 3.14159 ELSE 0 END)) as DOUGLAS_FIR_CROWN_AREA,
    AVG(CASE WHEN t.SPCD = 202 AND t.UNCRCD IS NOT NULL THEN t.UNCRCD END) as DOUGLAS_FIR_AVG_UNCRCD,
    AVG(CASE WHEN t.SPCD = 202 AND t.CR IS NOT NULL THEN t.CR END) as DOUGLAS_FIR_AVG_CR,
    
    -- Western Hemlock (Tsuga heterophylla) - SPCD 263
    SUM(CASE WHEN t.SPCD = 263 THEN 1 ELSE 0 END) as WESTERN_HEMLOCK_COUNT,
    ROUND(SUM(CASE WHEN t.SPCD = 263 THEN t.CARBON_AG ELSE 0 END), 0) as WESTERN_HEMLOCK_CARBON,
    ROUND(SUM(CASE WHEN t.SPCD = 263 THEN t.CARBON_BG ELSE 0 END), 0) as WESTERN_HEMLOCK_CARBON_BG,
    ROUND(SUM(CASE WHEN t.SPCD = 263 THEN t.DIA * t.DIA ELSE 0 END)) as WESTERN_HEMLOCK_BA,
    ROUND(SUM(CASE WHEN t.SPCD = 263 THEN POWER((0.3 * t.DIA + 1.2), 2) * 3.14159 ELSE 0 END)) as WESTERN_HEMLOCK_CROWN_AREA,
    AVG(CASE WHEN t.SPCD = 263 AND t.UNCRCD IS NOT NULL THEN t.UNCRCD END) as WESTERN_HEMLOCK_AVG_UNCRCD,
    
    -- Western Red Cedar (Thuja plicata) - SPCD 242
    SUM(CASE WHEN t.SPCD = 242 THEN 1 ELSE 0 END) as WESTERN_RED_CEDAR_COUNT,
    ROUND(SUM(CASE WHEN t.SPCD = 242 THEN t.CARBON_AG ELSE 0 END), 0) as WESTERN_RED_CEDAR_CARBON,
    ROUND(SUM(CASE WHEN t.SPCD = 242 THEN t.CARBON_BG ELSE 0 END), 0) as WESTERN_RED_CEDAR_CARBON_BG,
    ROUND(SUM(CASE WHEN t.SPCD = 242 THEN t.DIA * t.DIA ELSE 0 END)) as WESTERN_RED_CEDAR_BA,
    ROUND(SUM(CASE WHEN t.SPCD = 242 THEN POWER((0.3 * t.DIA + 1.2), 2) * 3.14159 ELSE 0 END)) as WESTERN_RED_CEDAR_CROWN_AREA,
    AVG(CASE WHEN t.SPCD = 242 AND t.UNCRCD IS NOT NULL THEN t.UNCRCD END) as WESTERN_RED_CEDAR_AVG_UNCRCD,
    
    -- Sitka Spruce (Picea sitchensis) - SPCD 098
    SUM(CASE WHEN t.SPCD = 98 THEN 1 ELSE 0 END) as SITKA_SPRUCE_COUNT,
    ROUND(SUM(CASE WHEN t.SPCD = 98 THEN t.CARBON_AG ELSE 0 END), 0) as SITKA_SPRUCE_CARBON,
    ROUND(SUM(CASE WHEN t.SPCD = 98 THEN POWER((0.25 * t.DIA + 1.0), 2) * 3.14159 ELSE 0 END)) as SITKA_SPRUCE_CROWN_AREA,
    
    -- Noble Fir (Abies procera) - SPCD 011
    SUM(CASE WHEN t.SPCD = 11 THEN 1 ELSE 0 END) as NOBLE_FIR_COUNT,
    ROUND(SUM(CASE WHEN t.SPCD = 11 THEN t.CARBON_AG ELSE 0 END), 0) as NOBLE_FIR_CARBON,
    ROUND(SUM(CASE WHEN t.SPCD = 11 THEN POWER((0.28 * t.DIA + 1.1), 2) * 3.14159 ELSE 0 END)) as NOBLE_FIR_CROWN_AREA,
    
    -- Pacific Silver Fir (Abies amabilis) - SPCD 017
    SUM(CASE WHEN t.SPCD = 17 THEN 1 ELSE 0 END) as PACIFIC_SILVER_FIR_COUNT,
    ROUND(SUM(CASE WHEN t.SPCD = 17 THEN t.CARBON_AG ELSE 0 END), 0) as PACIFIC_SILVER_FIR_CARBON,
    ROUND(SUM(CASE WHEN t.SPCD = 17 THEN POWER((0.28 * t.DIA + 1.1), 2) * 3.14159 ELSE 0 END)) as PACIFIC_SILVER_FIR_CROWN_AREA,
    
    
    -- Individual hardwood species for mushroom habitat analysis
    -- Bigleaf Maple (Acer macrophyllum) - SPCD 312
    SUM(CASE WHEN t.SPCD = 312 THEN 1 ELSE 0 END) as BIGLEAF_MAPLE_COUNT,
    ROUND(SUM(CASE WHEN t.SPCD = 312 THEN t.CARBON_AG ELSE 0 END), 0) as BIGLEAF_MAPLE_CARBON,
    ROUND(SUM(CASE WHEN t.SPCD = 312 THEN POWER((0.4 * t.DIA + 0.5), 2) * 3.14159 ELSE 0 END)) as BIGLEAF_MAPLE_CROWN_AREA,
    
    -- Red Alder (Alnus rubra) - SPCD 351
    SUM(CASE WHEN t.SPCD = 351 THEN 1 ELSE 0 END) as RED_ALDER_COUNT,
    ROUND(SUM(CASE WHEN t.SPCD = 351 THEN t.CARBON_AG ELSE 0 END), 0) as RED_ALDER_CARBON,
    ROUND(SUM(CASE WHEN t.SPCD = 351 THEN POWER((0.4 * t.DIA + 0.5), 2) * 3.14159 ELSE 0 END)) as RED_ALDER_CROWN_AREA,
    
    -- Black Cottonwood (Populus trichocarpa) - SPCD 746
    SUM(CASE WHEN t.SPCD = 746 THEN 1 ELSE 0 END) as BLACK_COTTONWOOD_COUNT,
    ROUND(SUM(CASE WHEN t.SPCD = 746 THEN t.CARBON_AG ELSE 0 END), 0) as BLACK_COTTONWOOD_CARBON,
    ROUND(SUM(CASE WHEN t.SPCD = 746 THEN POWER((0.4 * t.DIA + 0.5), 2) * 3.14159 ELSE 0 END)) as BLACK_COTTONWOOD_CROWN_AREA,
    
    -- Hardwood percentage (for forest type classification)
    SUM(CASE WHEN t.SPGRPCD >= 40 THEN 1 ELSE 0 END) as HARDWOOD_COUNT,
    
    -- Total metrics
    ROUND(SUM(t.CARBON_AG), 0) as TOTAL_CARBON_AG,
    ROUND(SUM(t.CARBON_BG), 0) as TOTAL_CARBON_BG,
    ROUND(SUM(t.DIA * t.DIA)) as TOTAL_BASAL_AREA,
    ROUND(SUM(POWER((0.3 * t.DIA + 1.2), 2) * 3.14159)) as TOTAL_CROWN_AREA,
    AVG(t.DIA) as MEAN_DBH,
    MAX(t.DIA) as MAX_DBH,
    AVG(t.UNCRCD) as MEAN_UNCRCD,
    AVG(t.CR) as MEAN_CR,
    AVG(t.ACTUALHT) as MEAN_HEIGHT,
    
    -- Species diversity
    COUNT(DISTINCT t.SPCD) as SPECIES_COUNT

FROM latest_plots lp
LEFT JOIN COND c ON lp.PLOT_CN = c.PLT_CN AND c.CONDID = 1  -- Main condition
LEFT JOIN TREE t ON lp.PLOT_CN = t.PLT_CN
WHERE lp.rn = 1  -- Only most recent survey for each plot location
GROUP BY lp.PLOT_CN, lp.STATECD, lp.UNITCD, lp.COUNTYCD, lp.PLOT, lp.LAT, lp.LON, lp.INVYR, lp.ELEV,
         c.FORTYPCD, c.FLDTYPCD, c.STDAGE, c.ADFORCD, c.SLOPE, c.ASPECT, c.LIVE_CANOPY_CVR_PCT
HAVING TOTAL_TREES > 0  -- Only plots with trees
ORDER BY TOTAL_CARBON_AG DESC;
"""

print("Executing comprehensive forest data query for WA, OR, and ID...")
print("This may take several minutes due to complex aggregations...")

# Execute query for each state and combine results
all_dataframes = []

for state, db_path in databases.items():
    print(f"\nProcessing {state} database...")
    try:
        conn = sqlite3.connect(db_path)
        df_state = pd.read_sql_query(query, conn)
        conn.close()
        print(f"  Extracted {len(df_state):,} plots from {state}")
        all_dataframes.append(df_state)
    except Exception as e:
        print(f"  Warning: Could not process {state} database: {e}")
        continue

# Combine all state dataframes
if all_dataframes:
    df = pd.concat(all_dataframes, ignore_index=True)
    print(f"\nCombined data: {len(df):,} total forested plots from {len(all_dataframes)} states")
else:
    raise Exception("No data could be extracted from any database")

# Calculate relative composition
df['DOUGLAS_FIR_PERCENT'] = (df['DOUGLAS_FIR_COUNT'] / df['TOTAL_TREES'] * 100).round(2)
df['WESTERN_HEMLOCK_PERCENT'] = (df['WESTERN_HEMLOCK_COUNT'] / df['TOTAL_TREES'] * 100).round(2)
df['WESTERN_RED_CEDAR_PERCENT'] = (df['WESTERN_RED_CEDAR_COUNT'] / df['TOTAL_TREES'] * 100).round(2)
df['SITKA_SPRUCE_PERCENT'] = (df['SITKA_SPRUCE_COUNT'] / df['TOTAL_TREES'] * 100).round(2)
df['NOBLE_FIR_PERCENT'] = (df['NOBLE_FIR_COUNT'] / df['TOTAL_TREES'] * 100).round(2)
df['PACIFIC_SILVER_FIR_PERCENT'] = (df['PACIFIC_SILVER_FIR_COUNT'] / df['TOTAL_TREES'] * 100).round(2)
df['BIGLEAF_MAPLE_PERCENT'] = (df['BIGLEAF_MAPLE_COUNT'] / df['TOTAL_TREES'] * 100).round(2)
df['RED_ALDER_PERCENT'] = (df['RED_ALDER_COUNT'] / df['TOTAL_TREES'] * 100).round(2)
df['BLACK_COTTONWOOD_PERCENT'] = (df['BLACK_COTTONWOOD_COUNT'] / df['TOTAL_TREES'] * 100).round(2)
df['HARDWOOD_PERCENT'] = (df['HARDWOOD_COUNT'] / df['TOTAL_TREES'] * 100).round(2)

# Define forest type based on dominant species
def classify_forest_type(row):
    if row['DOUGLAS_FIR_PERCENT'] >= 50:
        return 'Douglas-fir Dominant'
    elif row['WESTERN_HEMLOCK_PERCENT'] >= 50:
        return 'Western Hemlock Dominant'
    elif row['WESTERN_RED_CEDAR_PERCENT'] >= 50:
        return 'Western Red Cedar Dominant'
    elif row['HARDWOOD_PERCENT'] >= 50:
        return 'Hardwood Dominant'
    elif (row['DOUGLAS_FIR_PERCENT'] + row['WESTERN_HEMLOCK_PERCENT']) >= 50:
        return 'Douglas-fir/Western Hemlock Mixed'
    else:
        return 'Mixed Conifer'

df['FOREST_TYPE_PREDICTED'] = df.apply(classify_forest_type, axis=1)

# Create GeoDataFrame
geometry = [Point(xy) for xy in zip(df['LON'], df['LAT'])]
gdf = gpd.GeoDataFrame(df, geometry=geometry, crs='EPSG:4326')

# Export to GeoJSON
output_file = 'comprehensive_forest_data_pnw.geojson'
gdf.to_file(output_file, driver='GeoJSON')

print(f"Data exported to {output_file}")

# Display summary statistics
print(f"\nSummary Statistics:")
print(f"Total forested plots: {len(df):,}")
print(f"Elevation range: {df['ELEV'].min():.0f} to {df['ELEV'].max():.0f} feet")
print(f"Mean trees per plot: {df['TOTAL_TREES'].mean():.1f}")
print(f"Mean species per plot: {df['SPECIES_COUNT'].mean():.1f}")

print(f"\nForest Type Distribution:")
forest_type_counts = df['FOREST_TYPE_PREDICTED'].value_counts()
for forest_type, count in forest_type_counts.items():
    pct = 100 * count / len(df)
    print(f"  {forest_type}: {count:,} plots ({pct:.1f}%)")

print(f"\nSpecies Composition Summary:")
print(f"Plots with Douglas-fir: {(df['DOUGLAS_FIR_COUNT'] > 0).sum():,} ({100*(df['DOUGLAS_FIR_COUNT'] > 0).sum()/len(df):.2f}%)")
print(f"Plots with Western Hemlock: {(df['WESTERN_HEMLOCK_COUNT'] > 0).sum():,} ({100*(df['WESTERN_HEMLOCK_COUNT'] > 0).sum()/len(df):.2f}%)")
print(f"Plots with Western Red Cedar: {(df['WESTERN_RED_CEDAR_COUNT'] > 0).sum():,} ({100*(df['WESTERN_RED_CEDAR_COUNT'] > 0).sum()/len(df):.2f}%)")
print(f"Plots with Sitka Spruce: {(df['SITKA_SPRUCE_COUNT'] > 0).sum():,} ({100*(df['SITKA_SPRUCE_COUNT'] > 0).sum()/len(df):.2f}%)")
print(f"Plots with Noble Fir: {(df['NOBLE_FIR_COUNT'] > 0).sum():,} ({100*(df['NOBLE_FIR_COUNT'] > 0).sum()/len(df):.2f}%)")
print(f"Plots with Pacific Silver Fir: {(df['PACIFIC_SILVER_FIR_COUNT'] > 0).sum():,} ({100*(df['PACIFIC_SILVER_FIR_COUNT'] > 0).sum()/len(df):.2f}%)")
print(f"Plots with Bigleaf Maple: {(df['BIGLEAF_MAPLE_COUNT'] > 0).sum():,} ({100*(df['BIGLEAF_MAPLE_COUNT'] > 0).sum()/len(df):.2f}%)")
print(f"Plots with Red Alder: {(df['RED_ALDER_COUNT'] > 0).sum():,} ({100*(df['RED_ALDER_COUNT'] > 0).sum()/len(df):.2f}%)")
print(f"Plots with Black Cottonwood: {(df['BLACK_COTTONWOOD_COUNT'] > 0).sum():,} ({100*(df['BLACK_COTTONWOOD_COUNT'] > 0).sum()/len(df):.2f}%)")
print(f"Plots with Hardwoods (any): {(df['HARDWOOD_COUNT'] > 0).sum():,} ({100*(df['HARDWOOD_COUNT'] > 0).sum()/len(df):.2f}%)")

# Show some example records
print(f"\nExample High-Diversity Plots:")
high_diversity = df.nlargest(5, 'SPECIES_COUNT')[['PLOT_CN', 'LAT', 'LON', 'ELEV', 'SPECIES_COUNT', 'FOREST_TYPE_PREDICTED', 'TOTAL_TREES']]
print(high_diversity.to_string(index=False))

# Show slope/aspect distribution for mushroom habitat analysis
print(f"\nTopographic Summary (for mushroom habitat):")
print(f"Slope range: {df['SLOPE'].min():.0f}° to {df['SLOPE'].max():.0f}°")
print(f"Mean slope: {df['SLOPE'].mean():.1f}°")

# Aspect bins (for interpretability)
def assign_aspect(aspect_deg):
    if pd.isna(aspect_deg):
        return 'Unknown'
    elif aspect_deg <= 45 or aspect_deg > 315:
        return 'North'
    elif 45 < aspect_deg <= 135:
        return 'East'
    elif 135 < aspect_deg <= 225:
        return 'South'
    else:
        return 'West'

df['ASPECT_CARDINAL'] = df['ASPECT'].apply(assign_aspect)
aspect_counts = df['ASPECT_CARDINAL'].value_counts()
print(f"\nAspect distribution:")
for aspect, count in aspect_counts.items():
    pct = 100 * count / len(df)
    print(f"  {aspect}: {count:,} plots ({pct:.1f}%)") 