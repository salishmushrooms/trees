#!/usr/bin/env python3
"""
Create Species Distribution Maps based on Species Density/Dominance

This approach calculates the percentage that each species represents within 
each plot (e.g., Douglas Fir carbon / Total plot carbon) and creates
distribution areas starting from low density (5%) up to high dominance (80%+).

This is more ecologically meaningful for habitat mapping than absolute values.
"""

import sqlite3
import pandas as pd
import geopandas as gpd
import numpy as np
from pathlib import Path
from shapely.geometry import Point
from shapely.ops import unary_union
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

# Configuration
DB_PATH = "data/raw/trees_SQLite_FIADB_WA.db"
OUTPUT_DIR = Path("outputs/species_distribution")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Buffer parameters
BASE_BUFFER_KM = 1.5       # Base buffer radius in km
MAX_BUFFER_KM = 4.0        # Maximum buffer radius in km  
MIN_BUFFER_KM = 0.8        # Minimum buffer radius in km

def get_species_density_data(species_code, species_name):
    """Get species data with density percentage per plot"""
    print(f"\n🌲 Processing {species_name} (code: {species_code})")
    
    query = """
    WITH latest_plots AS (
        SELECT 
            STATECD, UNITCD, COUNTYCD, PLOT,
            MAX(INVYR) as LATEST_INVYR
        FROM PLOT 
        WHERE PLOT_STATUS_CD = 1
            AND LAT IS NOT NULL 
            AND LON IS NOT NULL
        GROUP BY STATECD, UNITCD, COUNTYCD, PLOT
    ),
    plot_totals AS (
        SELECT 
            p.CN as PLOT_CN,
            p.LAT,
            p.LON,
            p.ELEV as ELEVATION_FT,
            SUM(t.CARBON_AG) as TOTAL_PLOT_CARBON,
            COUNT(t.CN) as TOTAL_PLOT_TREES
        FROM PLOT p
        INNER JOIN latest_plots lp ON (
            p.STATECD = lp.STATECD AND 
            p.UNITCD = lp.UNITCD AND 
            p.COUNTYCD = lp.COUNTYCD AND 
            p.PLOT = lp.PLOT AND 
            p.INVYR = lp.LATEST_INVYR
        )
        INNER JOIN TREE t ON p.CN = t.PLT_CN
        WHERE t.STATUSCD = 1  -- Live trees only
            AND t.CARBON_AG IS NOT NULL
            AND t.CARBON_AG > 0
        GROUP BY p.CN, p.LAT, p.LON, p.ELEV
    ),
    species_data AS (
        SELECT 
            p.CN as PLOT_CN,
            p.LAT,
            p.LON,
            p.ELEV as ELEVATION_FT,
            COUNT(t.CN) as SPECIES_TREE_COUNT,
            SUM(t.CARBON_AG) as SPECIES_CARBON_AG
        FROM PLOT p
        INNER JOIN latest_plots lp ON (
            p.STATECD = lp.STATECD AND 
            p.UNITCD = lp.UNITCD AND 
            p.COUNTYCD = lp.COUNTYCD AND 
            p.PLOT = lp.PLOT AND 
            p.INVYR = lp.LATEST_INVYR
        )
        INNER JOIN TREE t ON p.CN = t.PLT_CN
        WHERE t.SPCD = ?
            AND t.STATUSCD = 1  -- Live trees only
            AND t.CARBON_AG IS NOT NULL
            AND t.CARBON_AG > 0
        GROUP BY p.CN, p.LAT, p.LON, p.ELEV
    )
    SELECT 
        pt.LAT,
        pt.LON,
        pt.ELEVATION_FT,
        sd.SPECIES_TREE_COUNT,
        sd.SPECIES_CARBON_AG,
        pt.TOTAL_PLOT_CARBON,
        pt.TOTAL_PLOT_TREES,
        ROUND((sd.SPECIES_CARBON_AG * 100.0 / pt.TOTAL_PLOT_CARBON), 2) as SPECIES_DENSITY_PCT,
        ROUND((sd.SPECIES_TREE_COUNT * 100.0 / pt.TOTAL_PLOT_TREES), 2) as SPECIES_COUNT_PCT
    FROM plot_totals pt
    INNER JOIN species_data sd ON pt.PLOT_CN = sd.PLOT_CN
    WHERE pt.TOTAL_PLOT_CARBON > 0
        AND (sd.SPECIES_CARBON_AG * 100.0 / pt.TOTAL_PLOT_CARBON) >= 5.0  -- Minimum 5% density
    ORDER BY SPECIES_DENSITY_PCT DESC;
    """
    
    conn = sqlite3.connect(DB_PATH)
    df = pd.read_sql_query(query, conn, params=[species_code])
    conn.close()
    
    print(f"  📍 Found {len(df):,} plots with ≥5% {species_name}")
    print(f"  🗺️  Extent: Lat {df.LAT.min():.3f}-{df.LAT.max():.3f}, Lon {df.LON.min():.3f}-{df.LON.max():.3f}")
    print(f"  📊 Density range: {df['SPECIES_DENSITY_PCT'].min():.1f}%-{df['SPECIES_DENSITY_PCT'].max():.1f}%")
    print(f"  📊 Mean density: {df['SPECIES_DENSITY_PCT'].mean():.1f}%")
    
    return df

def create_water_mask():
    """Create simple water body polygons for Washington state"""
    water_bodies = [
        # Puget Sound (simplified)
        [(-123.2, 47.0), (-122.2, 47.0), (-122.2, 48.5), (-123.2, 48.5)],
        # Lake Washington area  
        [(-122.3, 47.5), (-122.2, 47.5), (-122.2, 47.7), (-122.3, 47.7)],
        # Columbia River (much more realistic - only along southern border)
        [(-124.0, 45.5), (-121.0, 45.5), (-121.0, 46.3), (-124.0, 46.3)],  # Western portion
        [(-121.0, 45.5), (-117.0, 45.5), (-117.0, 46.1), (-121.0, 46.1)],  # Eastern portion - smaller
    ]
    
    from shapely.geometry import Polygon
    water_polygons = []
    for coords in water_bodies:
        if len(coords) >= 3:
            water_polygons.append(Polygon(coords))
    
    return water_polygons

def calculate_buffer_size(density_pct, tree_count):
    """Calculate buffer size based on species density percentage and tree count"""
    # Scale buffer size based on density percentage (5% = small, 80%+ = large)
    density_factor = min(density_pct / 80.0, 1.0)  # Cap at 80% for max buffer
    
    # Tree count factor (log scale)
    tree_factor = min(np.log10(tree_count + 1) / 3.0, 1.0)  # Cap at 1000 trees
    
    # Combine factors (weight density more heavily for habitat mapping)
    combined_score = 0.8 * density_factor + 0.2 * tree_factor
    
    # Scale to buffer range
    buffer_km = MIN_BUFFER_KM + (MAX_BUFFER_KM - MIN_BUFFER_KM) * combined_score
    
    return buffer_km

def create_density_distribution(species_df, species_name):
    """Create distribution using species density thresholds"""
    print(f"  🎯 Creating density-based distribution for {species_name}...")
    
    if len(species_df) < 3:
        print(f"  ⚠️  Insufficient data for {species_name} (need ≥3 plots)")
        return None
    
    # Define ecologically meaningful density thresholds
    density_thresholds = [
        (5, "Low presence (5-15%)"),
        (15, "Moderate presence (15-30%)"), 
        (30, "Common (30-50%)"),
        (50, "Abundant (50-70%)"),
        (70, "Dominant (70%+)")
    ]
    
    print(f"  📊 Density distribution:")
    for threshold, desc in density_thresholds:
        count = len(species_df[species_df['SPECIES_DENSITY_PCT'] >= threshold])
        print(f"    {desc}: {count:,} plots")
    
    # Create points and calculate buffer sizes
    points_with_buffers = []
    
    print(f"  🔧 Calculating buffers for {len(species_df)} plots...")
    
    for idx, row in tqdm(species_df.iterrows(), total=len(species_df), desc="  Processing plots"):
        lat, lon = row['LAT'], row['LON']
        density_pct = row['SPECIES_DENSITY_PCT']
        tree_count = row['SPECIES_TREE_COUNT']
        
        # Calculate buffer size based on density and tree count
        buffer_km = calculate_buffer_size(density_pct, tree_count)
        
        # Convert km to degrees (approximate)
        buffer_deg = buffer_km / 111.0
        
        # Create point and buffer
        point = Point(lon, lat)
        buffer_geom = point.buffer(buffer_deg)
        
        # Determine density category
        density_category = "low_5"
        for threshold, desc in density_thresholds:
            if density_pct >= threshold:
                density_category = f"density_{threshold}"
        
        points_with_buffers.append({
            'geometry': buffer_geom,
            'point_geom': point,
            'species': species_name,
            'density_pct': density_pct,
            'tree_count': tree_count,
            'buffer_km': buffer_km,
            'density_category': density_category,
            'latitude': lat,
            'longitude': lon,
            'elevation_ft': row.get('ELEVATION_FT', 0)
        })
    
    # Create GeoDataFrame
    print(f"  📐 Creating GeoDataFrame with {len(points_with_buffers)} buffered plots...")
    gdf = gpd.GeoDataFrame(points_with_buffers, crs='EPSG:4326')
    
    # Remove areas over water - DISABLED FOR NOW
    # water_polygons = create_water_mask()
    # if water_polygons:
    #     print(f"  🌊 Removing water areas...")
    #     water_union = unary_union(water_polygons)
    #     
    #     # Remove water areas from each buffer
    #     for idx in tqdm(gdf.index, desc="  Water masking"):
    #         buffer_geom = gdf.loc[idx, 'geometry']
    #         try:
    #             # Remove intersection with water
    #             land_geom = buffer_geom.difference(water_union)
    #             gdf.loc[idx, 'geometry'] = land_geom
    #         except Exception:
    #             # Keep original if difference fails
    #             continue
    
    print(f"  ⚠️  Water masking disabled - showing complete species distribution")
    
    # Create merged polygons for each density category
    print(f"  🔗 Creating merged distribution areas...")
    merged_polygons = []
    
    # Process in order from low to high density for proper layering
    category_order = ['density_5', 'density_15', 'density_30', 'density_50', 'density_70']
    
    for category in category_order:
        category_gdf = gdf[gdf['density_category'] == category]
        
        if len(category_gdf) > 0:
            threshold = int(category.split('_')[1])
            desc = next(desc for thresh, desc in density_thresholds if thresh == threshold)
            
            # Union all buffers in this category
            try:
                merged_geom = unary_union(category_gdf['geometry'].tolist())
                
                # Handle both single polygons and multipolygons
                if hasattr(merged_geom, 'geoms'):
                    # MultiPolygon - add each polygon separately
                    for geom in merged_geom.geoms:
                        if geom.area > 0.0005:  # Filter very small fragments
                            merged_polygons.append({
                                'geometry': geom,
                                'species': species_name,
                                'density_threshold': threshold,
                                'density_category': category,
                                'area_sq_deg': geom.area,
                                'level_desc': desc,
                                'opacity_suggested': min(0.3 + (threshold / 100.0), 0.8)  # 30-80% opacity
                            })
                else:
                    # Single Polygon
                    if merged_geom.area > 0.0005:
                        merged_polygons.append({
                            'geometry': merged_geom,
                            'species': species_name,
                            'density_threshold': threshold,
                            'density_category': category,
                            'area_sq_deg': merged_geom.area,
                            'level_desc': desc,
                            'opacity_suggested': min(0.3 + (threshold / 100.0), 0.8)  # 30-80% opacity
                        })
                        
            except Exception as e:
                print(f"    ⚠️  Failed to merge {category}: {e}")
                continue
    
    if not merged_polygons:
        print(f"  ⚠️  No merged polygons created for {species_name}")
        return None
    
    # Create final GeoDataFrame with merged areas
    result_gdf = gpd.GeoDataFrame(merged_polygons, crs='EPSG:4326')
    
    # Sort by density threshold (lowest first for proper layering in maps)
    result_gdf = result_gdf.sort_values('density_threshold').reset_index(drop=True)
    
    print(f"  ✅ Created {len(result_gdf)} distribution areas")
    total_area = result_gdf['area_sq_deg'].sum()
    print(f"  📏 Total distribution area: {total_area:.3f} square degrees")
    
    return result_gdf

def main():
    """Main execution function"""
    print("🌲 CREATING SPECIES DENSITY DISTRIBUTION MAPS")
    print("="*60)
    
    # Test with both species to show different patterns
    species_list = [
        (122, "Ponderosa Pine"),     # More constrained, interesting density patterns
        (202, "Douglas Fir"),        # Ubiquitous but density varies
    ]
    
    for species_code, species_name in species_list:
        try:
            # Get species density data (≥5% threshold)
            species_df = get_species_density_data(species_code, species_name)
            
            if len(species_df) == 0:
                print(f"  ⚠️  No plots found with ≥5% {species_name}")
                continue
            
            # Create density-based distribution
            species_gdf = create_density_distribution(species_df, species_name)
            
            if species_gdf is None:
                continue
            
            # Save to file
            output_file = OUTPUT_DIR / f"{species_name.lower().replace(' ', '_')}_density_distribution_wa.geojson"
            species_gdf.to_file(output_file, driver='GeoJSON')
            print(f"  💾 Saved: {output_file}")
            
            # Show extent check
            bounds = species_gdf.total_bounds
            print(f"  🗺️  Final extent: Lat {bounds[1]:.3f}-{bounds[3]:.3f}, Lon {bounds[0]:.3f}-{bounds[2]:.3f}")
            
        except Exception as e:
            print(f"  ❌ Error processing {species_name}: {e}")
            import traceback
            traceback.print_exc()
    
    print(f"\n✅ Species density distribution mapping complete!")
    print(f"\nFiles created in: {OUTPUT_DIR}")
    print(f"\nVisualization recommendations:")
    print(f"  1. Load GeoJSON in QGIS or MapBox")
    print(f"  2. Style by 'density_threshold' field")
    print(f"  3. Layer order: 5% (bottom, light) → 70% (top, dark)")
    print(f"  4. Use 'opacity_suggested' field for transparency")
    print(f"  5. Color scheme: Light green (5%) → Dark green (70%+)")
    print(f"  6. All areas ≥5% species presence will be shown")

if __name__ == "__main__":
    main() 