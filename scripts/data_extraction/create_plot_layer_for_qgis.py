import sqlite3
import pandas as pd
from shapely.geometry import Point
import geopandas as gpd

# Database paths for WA, OR, ID
databases = {
    'WA': '/Users/JJC/trees/data/raw/trees_SQLite_FIADB_WA.db',
    'OR': '/Users/JJC/trees/data/raw/trees_SQLite_FIADB_OR.db',
    'ID': '/Users/JJC/trees/data/raw/trees_SQLite_FIADB_ID.db'
}

# Lightweight query - just plot identifiers and basic metadata
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
        AND p.PLOT_STATUS_CD = 1
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

    -- Basic plot characteristics
    c.FORTYPCD,
    c.SLOPE,
    c.ASPECT,
    c.LIVE_CANOPY_CVR_PCT as CANOPY_COVER,

    -- Tree count for validation
    COUNT(t.CN) as TOTAL_TREES,

    -- Create unique plot identifier for lookup
    lp.STATECD || '_' || lp.UNITCD || '_' || lp.COUNTYCD || '_' || lp.PLOT || '_' || lp.INVYR as PLOT_ID

FROM latest_plots lp
LEFT JOIN COND c ON lp.PLOT_CN = c.PLT_CN AND c.CONDID = 1
LEFT JOIN TREE t ON lp.PLOT_CN = t.PLT_CN
WHERE lp.rn = 1
GROUP BY lp.PLOT_CN, lp.STATECD, lp.UNITCD, lp.COUNTYCD, lp.PLOT, lp.LAT, lp.LON,
         lp.INVYR, lp.ELEV, c.FORTYPCD, c.SLOPE, c.ASPECT, c.LIVE_CANOPY_CVR_PCT
HAVING TOTAL_TREES > 0
ORDER BY lp.STATECD, lp.UNITCD, lp.COUNTYCD, lp.PLOT;
"""

print("Creating lightweight plot layer for QGIS...")
print("This extracts plot locations and identifiers without species data for faster loading")

# Execute query for each state and combine results
all_dataframes = []
state_counts = {}

for state, db_path in databases.items():
    print(f"\nProcessing {state} database...")
    try:
        conn = sqlite3.connect(db_path)
        df_state = pd.read_sql_query(query, conn)
        conn.close()

        # Add state identifier for easy filtering in QGIS
        df_state['STATE'] = state

        print(f"  Extracted {len(df_state):,} plots from {state}")
        state_counts[state] = len(df_state)
        all_dataframes.append(df_state)
    except Exception as e:
        print(f"  Warning: Could not process {state} database: {e}")
        continue

# Combine all state dataframes
if all_dataframes:
    df = pd.concat(all_dataframes, ignore_index=True)
    print(f"\nCombined data: {len(df):,} total forested plots")

    # Show state breakdown
    for state, count in state_counts.items():
        pct = 100 * count / len(df)
        print(f"  {state}: {count:,} plots ({pct:.1f}%)")

else:
    raise Exception("No data could be extracted from any database")

# Create GeoDataFrame
geometry = [Point(xy) for xy in zip(df['LON'], df['LAT'])]
gdf = gpd.GeoDataFrame(df, geometry=geometry, crs='EPSG:4326')

# Export to GeoJSON for QGIS
output_file = '/Users/JJC/trees/outputs/fia_plots_for_qgis.geojson'
gdf.to_file(output_file, driver='GeoJSON')

print(f"\nLightweight plot layer exported to: {output_file}")
print(f"File size should be much smaller than comprehensive extraction")

# Display summary
print(f"\nLayer Statistics:")
print(f"Total plots: {len(df):,}")
print(f"Elevation range: {df['ELEV'].min():.0f} to {df['ELEV'].max():.0f} feet")
print(f"Tree count range: {df['TOTAL_TREES'].min():.0f} to {df['TOTAL_TREES'].max():.0f} trees per plot")
print(f"Mean trees per plot: {df['TOTAL_TREES'].mean():.1f}")

print(f"\nAttributes in QGIS layer:")
for col in df.columns:
    if col not in ['geometry']:
        print(f"  {col}")

print(f"\nNext steps:")
print(f"1. Load {output_file} into QGIS")
print(f"2. Use selection tools to select plots of interest")
print(f"3. Run the species lookup script with selected PLOT_IDs")