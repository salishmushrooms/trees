import sqlite3
import json
import pandas as pd
from shapely.geometry import Point
import geopandas as gpd

# Connect to the Washington FIA database
conn = sqlite3.connect('trees_SQLite_FIADB_WA.db')

# Query to extract comprehensive forest data for mushroom habitat modeling
query = """
SELECT 
    p.CN as PLOT_CN,
    p.STATECD,
    p.UNITCD, 
    p.COUNTYCD,
    p.PLOT,
    p.LAT,
    p.LON,
    p.INVYR,
    p.ELEV,
    
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
    
    -- Douglas-fir metrics
    SUM(CASE WHEN t.SPCD = 202 THEN 1 ELSE 0 END) as DOUGLAS_FIR_COUNT,
    SUM(CASE WHEN t.SPCD = 202 THEN t.CARBON_AG ELSE 0 END) as DOUGLAS_FIR_CARBON,
    SUM(CASE WHEN t.SPCD = 202 THEN t.DIA * t.DIA ELSE 0 END) as DOUGLAS_FIR_BA,
    
    -- Western Hemlock (Tsuga heterophylla) - SPCD 263
    SUM(CASE WHEN t.SPCD = 263 THEN 1 ELSE 0 END) as HEMLOCK_COUNT,
    SUM(CASE WHEN t.SPCD = 263 THEN t.CARBON_AG ELSE 0 END) as HEMLOCK_CARBON,
    SUM(CASE WHEN t.SPCD = 263 THEN t.DIA * t.DIA ELSE 0 END) as HEMLOCK_BA,
    
    -- Western Red Cedar (Thuja plicata) - SPCD 242
    SUM(CASE WHEN t.SPCD = 242 THEN 1 ELSE 0 END) as CEDAR_COUNT,
    SUM(CASE WHEN t.SPCD = 242 THEN t.CARBON_AG ELSE 0 END) as CEDAR_CARBON,
    SUM(CASE WHEN t.SPCD = 242 THEN t.DIA * t.DIA ELSE 0 END) as CEDAR_BA,
    
    -- Sitka Spruce (Picea sitchensis) - SPCD 098
    SUM(CASE WHEN t.SPCD = 98 THEN 1 ELSE 0 END) as SITKA_SPRUCE_COUNT,
    SUM(CASE WHEN t.SPCD = 98 THEN t.CARBON_AG ELSE 0 END) as SITKA_SPRUCE_CARBON,
    
    -- Noble Fir (Abies procera) - SPCD 011
    SUM(CASE WHEN t.SPCD = 11 THEN 1 ELSE 0 END) as NOBLE_FIR_COUNT,
    SUM(CASE WHEN t.SPCD = 11 THEN t.CARBON_AG ELSE 0 END) as NOBLE_FIR_CARBON,
    
    -- Pacific Silver Fir (Abies amabilis) - SPCD 017
    SUM(CASE WHEN t.SPCD = 17 THEN 1 ELSE 0 END) as SILVER_FIR_COUNT,
    SUM(CASE WHEN t.SPCD = 17 THEN t.CARBON_AG ELSE 0 END) as SILVER_FIR_CARBON,
    
    -- All firs combined
    SUM(CASE WHEN t.SPCD BETWEEN 10 AND 19 THEN 1 ELSE 0 END) as ALL_FIR_COUNT,
    SUM(CASE WHEN t.SPCD BETWEEN 10 AND 19 THEN t.CARBON_AG ELSE 0 END) as ALL_FIR_CARBON,
    
    -- All pines 
    SUM(CASE WHEN t.SPCD BETWEEN 100 AND 199 THEN 1 ELSE 0 END) as ALL_PINE_COUNT,
    SUM(CASE WHEN t.SPCD BETWEEN 100 AND 199 THEN t.CARBON_AG ELSE 0 END) as ALL_PINE_CARBON,
    
    -- Hardwoods
    SUM(CASE WHEN t.SPGRPCD >= 40 THEN 1 ELSE 0 END) as HARDWOOD_COUNT,
    SUM(CASE WHEN t.SPGRPCD >= 40 THEN t.CARBON_AG ELSE 0 END) as HARDWOOD_CARBON,
    
    -- Total metrics
    SUM(t.CARBON_AG) as TOTAL_CARBON_AG,
    SUM(t.DIA * t.DIA) as TOTAL_BASAL_AREA,
    AVG(t.DIA) as MEAN_DBH,
    MAX(t.DIA) as MAX_DBH,
    
    -- Species diversity
    COUNT(DISTINCT t.SPCD) as SPECIES_COUNT

FROM PLOT p
LEFT JOIN COND c ON p.CN = c.PLT_CN AND c.CONDID = 1  -- Main condition
LEFT JOIN TREE t ON p.CN = t.PLT_CN
WHERE p.LAT IS NOT NULL 
    AND p.LON IS NOT NULL
    AND p.PLOT_STATUS_CD = 1  -- Forested plots only
GROUP BY p.CN, p.STATECD, p.UNITCD, p.COUNTYCD, p.PLOT, p.LAT, p.LON, p.INVYR, p.ELEV,
         c.FORTYPCD, c.FLDTYPCD, c.STDAGE, c.ADFORCD, c.SLOPE, c.ASPECT, c.LIVE_CANOPY_CVR_PCT
HAVING TOTAL_TREES > 0  -- Only plots with trees
ORDER BY TOTAL_CARBON_AG DESC;
"""

print("Executing comprehensive forest data query...")
print("This may take a few minutes due to complex aggregations...")

df = pd.read_sql_query(query, conn)
conn.close()

print(f"Extracted data for {len(df)} forested plots")

# Calculate relative composition
df['DOUGLAS_FIR_PERCENT'] = (df['DOUGLAS_FIR_COUNT'] / df['TOTAL_TREES'] * 100).round(1)
df['HEMLOCK_PERCENT'] = (df['HEMLOCK_COUNT'] / df['TOTAL_TREES'] * 100).round(1)
df['CEDAR_PERCENT'] = (df['CEDAR_COUNT'] / df['TOTAL_TREES'] * 100).round(1)
df['HARDWOOD_PERCENT'] = (df['HARDWOOD_COUNT'] / df['TOTAL_TREES'] * 100).round(1)

# Define forest type based on dominant species
def classify_forest_type(row):
    if row['DOUGLAS_FIR_PERCENT'] >= 50:
        return 'Douglas-fir Dominant'
    elif row['HEMLOCK_PERCENT'] >= 50:
        return 'Western Hemlock Dominant'
    elif row['CEDAR_PERCENT'] >= 50:
        return 'Western Red Cedar Dominant'
    elif row['HARDWOOD_PERCENT'] >= 50:
        return 'Hardwood Dominant'
    elif (row['DOUGLAS_FIR_PERCENT'] + row['HEMLOCK_PERCENT']) >= 50:
        return 'Douglas-fir/Hemlock Mixed'
    elif row['ALL_FIR_COUNT'] >= row['TOTAL_TREES'] * 0.5:
        return 'Fir Species Dominant'
    elif row['ALL_PINE_COUNT'] >= row['TOTAL_TREES'] * 0.5:
        return 'Pine Species Dominant'
    else:
        return 'Mixed Conifer'

df['FOREST_TYPE_PREDICTED'] = df.apply(classify_forest_type, axis=1)

# Create GeoDataFrame
geometry = [Point(xy) for xy in zip(df['LON'], df['LAT'])]
gdf = gpd.GeoDataFrame(df, geometry=geometry, crs='EPSG:4326')

# Export to GeoJSON
output_file = 'comprehensive_forest_data_wa.geojson'
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
print(f"Plots with Douglas-fir: {(df['DOUGLAS_FIR_COUNT'] > 0).sum():,} ({100*(df['DOUGLAS_FIR_COUNT'] > 0).sum()/len(df):.1f}%)")
print(f"Plots with Western Hemlock: {(df['HEMLOCK_COUNT'] > 0).sum():,} ({100*(df['HEMLOCK_COUNT'] > 0).sum()/len(df):.1f}%)")
print(f"Plots with Western Red Cedar: {(df['CEDAR_COUNT'] > 0).sum():,} ({100*(df['CEDAR_COUNT'] > 0).sum()/len(df):.1f}%)")
print(f"Plots with Hardwoods: {(df['HARDWOOD_COUNT'] > 0).sum():,} ({100*(df['HARDWOOD_COUNT'] > 0).sum()/len(df):.1f}%)")

# Show some example records
print(f"\nExample High-Diversity Plots:")
high_diversity = df.nlargest(5, 'SPECIES_COUNT')[['PLOT_CN', 'LAT', 'LON', 'ELEV', 'SPECIES_COUNT', 'FOREST_TYPE_PREDICTED', 'TOTAL_TREES']]
print(high_diversity.to_string(index=False))

# Show slope/aspect distribution for mushroom habitat analysis
print(f"\nTopographic Summary (for mushroom habitat):")
print(f"Slope range: {df['SLOPE'].min():.0f}° to {df['SLOPE'].max():.0f}°")
print(f"Mean slope: {df['SLOPE'].mean():.1f}°")

# Aspect bins (for interpretability)
df['ASPECT_CARDINAL'] = pd.cut(df['ASPECT'], 
                               bins=[0, 45, 135, 225, 315, 360], 
                               labels=['North', 'East', 'South', 'West', 'North'],
                               include_lowest=True)
aspect_counts = df['ASPECT_CARDINAL'].value_counts()
print(f"\nAspect distribution:")
for aspect, count in aspect_counts.items():
    pct = 100 * count / len(df)
    print(f"  {aspect}: {count:,} plots ({pct:.1f}%)") 