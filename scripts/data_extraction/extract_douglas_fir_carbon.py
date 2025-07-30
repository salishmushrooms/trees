import sqlite3
import json
import pandas as pd
from shapely.geometry import Point
import geopandas as gpd

# Connect to the Washington FIA database
conn = sqlite3.connect('trees_SQLite_FIADB_WA.db')

# Query to extract Douglas-fir carbon data aggregated by plot
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
    COUNT(t.CN) as DOUGLAS_FIR_COUNT,
    SUM(t.CARBON_AG) as TOTAL_CARBON_AG,
    AVG(t.CARBON_AG) as AVG_CARBON_AG,
    MIN(t.CARBON_AG) as MIN_CARBON_AG,
    MAX(t.CARBON_AG) as MAX_CARBON_AG
FROM PLOT p
INNER JOIN TREE t ON p.CN = t.PLT_CN
WHERE t.SPCD = 202 
    AND t.CARBON_AG IS NOT NULL 
    AND p.LAT IS NOT NULL 
    AND p.LON IS NOT NULL
GROUP BY p.CN, p.STATECD, p.UNITCD, p.COUNTYCD, p.PLOT, p.LAT, p.LON, p.INVYR
ORDER BY TOTAL_CARBON_AG DESC;
"""

print("Executing query to extract Douglas-fir carbon data...")
df = pd.read_sql_query(query, conn)
conn.close()

print(f"Extracted data for {len(df)} plots with Douglas-fir trees")
print(f"Total carbon range: {df['TOTAL_CARBON_AG'].min():.2f} to {df['TOTAL_CARBON_AG'].max():.2f}")
print(f"Average carbon per plot: {df['TOTAL_CARBON_AG'].mean():.2f}")

# Create GeoDataFrame
geometry = [Point(xy) for xy in zip(df['LON'], df['LAT'])]
gdf = gpd.GeoDataFrame(df, geometry=geometry, crs='EPSG:4326')

# Export to GeoJSON
output_file = 'douglas_fir_carbon_points_wa.geojson'
gdf.to_file(output_file, driver='GeoJSON')

print(f"Data exported to {output_file}")
print(f"Coordinate bounds:")
print(f"  Latitude: {df['LAT'].min():.6f} to {df['LAT'].max():.6f}")
print(f"  Longitude: {df['LON'].min():.6f} to {df['LON'].max():.6f}")

# Display summary statistics
print("\nSummary Statistics:")
print(df[['DOUGLAS_FIR_COUNT', 'TOTAL_CARBON_AG', 'AVG_CARBON_AG']].describe())

# Show top 10 plots by carbon
print("\nTop 10 plots by total above-ground carbon:")
print(df[['PLOT_CN', 'COUNTYCD', 'LAT', 'LON', 'DOUGLAS_FIR_COUNT', 'TOTAL_CARBON_AG']].head(10)) 