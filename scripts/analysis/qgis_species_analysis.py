# QGIS Python Console Script - Species Summary for Text Fields
# Copy and paste this entire script into the QGIS Python console
# Outputs simple format for copying into QGIS text fields

import sqlite3
from collections import defaultdict

# Database paths
databases = {
    'WA': '/Users/JJC/trees/data/raw/trees_SQLite_FIADB_WA.db',
    'OR': '/Users/JJC/trees/data/raw/trees_SQLite_FIADB_OR.db',
    'ID': '/Users/JJC/trees/data/raw/trees_SQLite_FIADB_ID.db'
}

# Species lookup will be loaded dynamically from database

def capitalize_species_name(species_name):
    """Apply proper capitalization to species names for iOS app compatibility"""
    if not species_name:
        return species_name

    # Convert to title case, but handle special cases
    capitalized = species_name.title()

    # Handle special words that should remain lowercase
    lowercase_words = ['spp.', 'or', 'and', 'of']
    words = capitalized.split()

    for i, word in enumerate(words):
        # Keep first word capitalized, handle others based on rules
        if i > 0 and word.lower() in lowercase_words:
            words[i] = word.lower()
        # Handle hyphenated words
        elif '-' in word:
            parts = word.split('-')
            # Capitalize each part of hyphenated words
            words[i] = '-'.join(part.capitalize() for part in parts)

    return ' '.join(words)

def is_manual_data(trees_field):
    """
    Check if tree data is manually entered based on 'nearby' indicator
    Manual data format: "...nearby sample plots)"
    Automated data format: "...sample plots)" (without 'nearby')
    """
    if not trees_field:
        return False

    trees_str = str(trees_field).strip()

    # Manual data contains "nearby sample plots"
    return "nearby sample plots)" in trees_str

def parse_plot_id(plot_id):
    """Parse PLOT_ID string into components"""
    parts = plot_id.split('_')
    return {
        'STATECD': int(parts[0]), 'UNITCD': int(parts[1]),
        'COUNTYCD': int(parts[2]), 'PLOT': int(parts[3]), 'INVYR': int(parts[4])
    }

def get_state_from_statecd(statecd):
    """Map state code to state abbreviation"""
    return {41: 'OR', 53: 'WA', 16: 'ID'}.get(statecd)

def load_species_mapping():
    """Load species code to common name mapping from first available database"""
    species_query = """
    SELECT SPCD, COMMON_NAME
    FROM REF_SPECIES
    WHERE COMMON_NAME IS NOT NULL
    ORDER BY SPCD;
    """

    for state, db_path in databases.items():
        try:
            conn = sqlite3.connect(db_path)
            cursor = conn.cursor()
            cursor.execute(species_query)
            results = cursor.fetchall()
            conn.close()

            species_map = {}
            for spcd, common_name in results:
                if common_name:
                    # Apply proper capitalization
                    species_map[spcd] = capitalize_species_name(common_name.strip())

            print(f"Loaded {len(species_map)} species mappings from {state} database")
            return species_map

        except Exception as e:
            print(f"Warning: Could not load species mapping from {state}: {e}")
            continue

    # Fallback to basic mapping if database lookup fails
    print("Warning: Using fallback species mapping")
    return {
        202: 'Douglas-fir', 263: 'Western Hemlock', 242: 'Western Red Cedar',
        122: 'Ponderosa Pine', 22: 'Grand Fir', 312: 'Bigleaf Maple',
        351: 'Red Alder', 17: 'Pacific Silver Fir'
    }

def analyze_selected_plots():
    """Main analysis function"""
    # Get the active layer
    layer = iface.activeLayer()

    if layer is None:
        print("ERROR: No active layer! Please click on the plot layer first.")
        return

    if layer.selectedFeatureCount() == 0:
        print("ERROR: No features selected! Please select some plots first.")
        return

    # Load species mapping from database
    print("Loading species names from database...")
    species_names = load_species_mapping()

    # Get selected plot IDs
    selected_features = layer.selectedFeatures()
    plot_ids = [str(f['PLOT_ID']) for f in selected_features if f['PLOT_ID']]

    if not plot_ids:
        print("ERROR: No valid PLOT_IDs found in selected features!")
        return

    print(f"Analyzing {len(plot_ids)} selected plots...")

    # Aggregate results across all plots
    species_totals = defaultdict(lambda: {
        'tree_count': 0, 'total_carbon': 0, 'total_basal_area': 0,
        'diameter_sum': 0, 'diameter_count': 0, 'plot_count': 0
    })

    total_plots_processed = 0
    total_trees = 0

    # Process each plot
    for i, plot_id in enumerate(plot_ids):
        if (i + 1) % 10 == 0:
            print(f"  Processing plot {i + 1}/{len(plot_ids)}...")

        try:
            plot_info = parse_plot_id(plot_id)
            state = get_state_from_statecd(plot_info['STATECD'])

            if not state or state not in databases:
                continue

            # Query for this plot
            query = """
            WITH plot_data AS (
                SELECT p.CN as PLOT_CN
                FROM PLOT p
                WHERE p.STATECD = ? AND p.UNITCD = ? AND p.COUNTYCD = ?
                    AND p.PLOT = ? AND p.INVYR = ? AND p.PLOT_STATUS_CD = 1
            )
            SELECT t.SPCD, COUNT(*) as TREE_COUNT,
                   ROUND(SUM(t.CARBON_AG), 0) as TOTAL_CARBON,
                   ROUND(SUM(t.DIA * t.DIA)) as BASAL_AREA,
                   SUM(t.DIA) as DIA_SUM, COUNT(t.DIA) as DIA_COUNT
            FROM plot_data pd
            JOIN TREE t ON pd.PLOT_CN = t.PLT_CN
            WHERE t.SPCD IS NOT NULL
            GROUP BY t.SPCD;
            """

            conn = sqlite3.connect(databases[state])
            cursor = conn.cursor()
            cursor.execute(query, [
                plot_info['STATECD'], plot_info['UNITCD'], plot_info['COUNTYCD'],
                plot_info['PLOT'], plot_info['INVYR']
            ])

            results = cursor.fetchall()
            conn.close()

            if results:
                total_plots_processed += 1
                plot_species = set()

                for spcd, tree_count, carbon, basal_area, dia_sum, dia_count in results:
                    species_name = species_names.get(spcd, f'Other Species')
                    # Species name already capitalized from database lookup

                    species_totals[species_name]['tree_count'] += tree_count
                    species_totals[species_name]['total_carbon'] += carbon or 0
                    species_totals[species_name]['total_basal_area'] += basal_area or 0
                    species_totals[species_name]['diameter_sum'] += dia_sum or 0
                    species_totals[species_name]['diameter_count'] += dia_count or 0

                    if species_name not in plot_species:
                        species_totals[species_name]['plot_count'] += 1
                        plot_species.add(species_name)

                    total_trees += tree_count

        except Exception as e:
            print(f"    Warning: Error processing plot {plot_id}: {e}")
            continue

    if not species_totals:
        print("No tree data found for selected plots!")
        return

    # Convert to sorted results
    results = []
    for species_name, data in species_totals.items():
        avg_diameter = data['diameter_sum'] / data['diameter_count'] if data['diameter_count'] > 0 else 0
        results.append({
            'Species': species_name,
            'Tree_Count': data['tree_count'],
            'Plot_Count': data['plot_count'],
            'Total_Carbon': data['total_carbon'],
            'Total_Basal_Area': data['total_basal_area'],
            'Avg_Diameter': round(avg_diameter, 1)
        })

    # Sort by tree count (most abundant first)
    results.sort(key=lambda x: x['Tree_Count'], reverse=True)

    # Filter for species >5% and create copy-paste format
    significant_species = []
    for result in results:
        tree_pct = 100 * result['Tree_Count'] / total_trees
        if tree_pct >= 5.0:
            significant_species.append((result['Species'], tree_pct))

    print(f"\n{'='*60}")
    print(f"SPECIES SUMMARY: {total_plots_processed} plots, {total_trees:,} trees")
    print(f"{'='*60}")

    if significant_species:
        print(f"\nSpecies >5% (iOS app format):")
        print("-" * 50)

        # Format in iOS app style (comma-separated with proper capitalization)
        species_list = []
        for species, pct in significant_species:
            # Use proper capitalization (already applied via capitalize_species_name)
            species_list.append(f"{species} ({pct:.0f}%)")

        # Create iOS-compatible format
        ios_format = ", ".join(species_list)
        ios_format += f" (Analysis based on {total_plots_processed} nearby sample plots)"

        print(f"\n" + "-" * 50)
        print("COPY THIS TEXT (iOS App Format):")
        print("-" * 50)
        print(ios_format)
        print("-" * 50)

    else:
        print("\nNo species found with >5% abundance")

# Run the analysis
analyze_selected_plots()