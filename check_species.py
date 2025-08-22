import pandas as pd
import geopandas as gpd

# Load the data
gdf = gpd.read_file('/Users/JJC/communications/tree-climate-mushroom/data/plot_data_with_crown_areas.geojson')

# Check species columns
count_columns = [col for col in gdf.columns if '_COUNT' in col and col != 'SPECIES_COUNT']
print('Species with count data:')
for col in count_columns:
    species_name = col.replace('_COUNT', '')
    total_plots = (gdf[col] > 0).sum()
    max_count = gdf[col].max()
    print(f'  {species_name}: {total_plots} plots, max {max_count} trees')

print(f'\nTotal plots: {len(gdf)}')

# Check which species have crown area data
crown_columns = [col for col in gdf.columns if '_CROWN_AREA' in col]
print('\nSpecies with crown area data:')
for col in crown_columns:
    species_name = col.replace('_CROWN_AREA', '')
    plots_with_crown = (gdf[col] > 0).sum()
    print(f'  {species_name}: {plots_with_crown} plots')