#!/usr/bin/env python3
"""
Species Lookup Tool for QGIS-Selected Plots

This script takes plot IDs (from QGIS selections) and dynamically queries
the FIA databases to return the top 10 species by tree count for those plots.

Usage:
1. Select features in QGIS
2. Export selected PLOT_IDs to a text file or CSV
3. Run this script with the plot IDs
4. Get species composition results

Example:
    python species_lookup_tool.py --plot-ids "41_6_33_123_2019,53_2_5_456_2020"
    python species_lookup_tool.py --file selected_plots.txt
"""

import sqlite3
import pandas as pd
import argparse
import sys
from pathlib import Path

# Database paths
databases = {
    'WA': '/Users/JJC/trees/data/raw/trees_SQLite_FIADB_WA.db',
    'OR': '/Users/JJC/trees/data/raw/trees_SQLite_FIADB_OR.db',
    'ID': '/Users/JJC/trees/data/raw/trees_SQLite_FIADB_ID.db'
}

# Species name lookup (SPCD to common name)
species_names = {
    11: 'Noble Fir',
    17: 'Pacific Silver Fir',
    19: 'Subalpine Fir',
    15: 'White Fir',
    22: 'Grand Fir',
    98: 'Sitka Spruce',
    93: 'Engelmann Spruce',
    108: 'Lodgepole Pine',
    122: 'Ponderosa Pine',
    202: 'Douglas-fir',
    242: 'Western Red Cedar',
    263: 'Western Hemlock',
    264: 'Mountain Hemlock',
    312: 'Bigleaf Maple',
    351: 'Red Alder',
    746: 'Black Cottonwood',
    748: 'Quaking Aspen',
    805: 'Pacific Madrone',
    839: 'Tanoak',
    64: 'Western Juniper',
    475: 'Oregon White Oak',
    981: 'Paper Birch'
}

def parse_plot_id(plot_id):
    """Parse PLOT_ID string into components for database lookup"""
    parts = plot_id.split('_')
    if len(parts) != 5:
        raise ValueError(f"Invalid PLOT_ID format: {plot_id}. Expected: STATECD_UNITCD_COUNTYCD_PLOT_INVYR")

    return {
        'STATECD': int(parts[0]),
        'UNITCD': int(parts[1]),
        'COUNTYCD': int(parts[2]),
        'PLOT': int(parts[3]),
        'INVYR': int(parts[4])
    }

def get_state_from_statecd(statecd):
    """Map state code to state abbreviation"""
    state_map = {41: 'OR', 53: 'WA', 16: 'ID'}
    return state_map.get(statecd)

def lookup_species_for_plot(plot_id):
    """Look up all species and counts for a single plot"""
    try:
        plot_info = parse_plot_id(plot_id)
        state = get_state_from_statecd(plot_info['STATECD'])

        if not state or state not in databases:
            print(f"Warning: Unknown state for plot {plot_id}")
            return None

        # Query to get all species for this specific plot
        query = """
        WITH plot_data AS (
            SELECT p.CN as PLOT_CN
            FROM PLOT p
            WHERE p.STATECD = ?
                AND p.UNITCD = ?
                AND p.COUNTYCD = ?
                AND p.PLOT = ?
                AND p.INVYR = ?
                AND p.PLOT_STATUS_CD = 1
        )
        SELECT
            t.SPCD,
            COUNT(*) as TREE_COUNT,
            ROUND(SUM(t.CARBON_AG), 0) as TOTAL_CARBON,
            ROUND(AVG(t.DIA), 1) as AVG_DIAMETER,
            ROUND(SUM(t.DIA * t.DIA)) as BASAL_AREA
        FROM plot_data pd
        JOIN TREE t ON pd.PLOT_CN = t.PLT_CN
        WHERE t.SPCD IS NOT NULL
        GROUP BY t.SPCD
        ORDER BY TREE_COUNT DESC;
        """

        # Execute query
        conn = sqlite3.connect(databases[state])
        df = pd.read_sql_query(query, conn, params=[
            plot_info['STATECD'],
            plot_info['UNITCD'],
            plot_info['COUNTYCD'],
            plot_info['PLOT'],
            plot_info['INVYR']
        ])
        conn.close()

        if df.empty:
            print(f"Warning: No trees found for plot {plot_id}")
            return None

        # Add species names
        df['SPECIES_NAME'] = df['SPCD'].map(species_names)
        df['SPECIES_NAME'] = df['SPECIES_NAME'].fillna(f'Unknown Species')

        # Add plot info
        df['PLOT_ID'] = plot_id
        df['STATE'] = state

        return df

    except Exception as e:
        print(f"Error processing plot {plot_id}: {e}")
        return None

def lookup_multiple_plots(plot_ids, top_n=10):
    """Look up species data for multiple plots"""
    all_results = []

    print(f"Looking up species data for {len(plot_ids)} plots...")

    for i, plot_id in enumerate(plot_ids, 1):
        print(f"Processing plot {i}/{len(plot_ids)}: {plot_id}")

        result = lookup_species_for_plot(plot_id.strip())
        if result is not None and not result.empty:
            # Take top N species for this plot
            top_species = result.head(top_n).copy()
            all_results.append(top_species)

    if not all_results:
        print("No data found for any plots")
        return None

    # Combine all results
    combined_df = pd.concat(all_results, ignore_index=True)

    return combined_df

def summarize_results(df):
    """Create summary statistics from lookup results"""
    if df is None or df.empty:
        return

    print(f"\n{'='*60}")
    print(f"SPECIES LOOKUP RESULTS")
    print(f"{'='*60}")

    total_plots = df['PLOT_ID'].nunique()
    total_trees = df['TREE_COUNT'].sum()

    print(f"Plots analyzed: {total_plots}")
    print(f"Total trees: {total_trees:,}")

    # Species frequency across plots
    print(f"\nTop Species by Frequency Across Plots:")
    species_freq = df.groupby('SPECIES_NAME')['PLOT_ID'].nunique().sort_values(ascending=False)
    for species, count in species_freq.head(10).items():
        pct = 100 * count / total_plots
        print(f"  {species}: {count}/{total_plots} plots ({pct:.1f}%)")

    # Species by total tree count
    print(f"\nTop Species by Total Tree Count:")
    species_trees = df.groupby('SPECIES_NAME')['TREE_COUNT'].sum().sort_values(ascending=False)
    for species, count in species_trees.head(10).items():
        pct = 100 * count / total_trees
        print(f"  {species}: {count:,} trees ({pct:.1f}%)")

    # State breakdown
    print(f"\nPlots by State:")
    state_counts = df.groupby('STATE')['PLOT_ID'].nunique()
    for state, count in state_counts.items():
        print(f"  {state}: {count} plots")

def main():
    parser = argparse.ArgumentParser(
        description='Look up top species for FIA plots selected in QGIS',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single plot
  python species_lookup_tool.py --plot-ids "41_6_33_123_2019"

  # Multiple plots (comma-separated)
  python species_lookup_tool.py --plot-ids "41_6_33_123_2019,53_2_5_456_2020"

  # From file (one PLOT_ID per line)
  python species_lookup_tool.py --file selected_plots.txt

  # Export results to CSV
  python species_lookup_tool.py --file selected_plots.txt --output results.csv
        """
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--plot-ids', type=str,
                      help='Comma-separated list of PLOT_IDs')
    group.add_argument('--file', type=str,
                      help='File containing PLOT_IDs (one per line)')

    parser.add_argument('--top-n', type=int, default=10,
                       help='Number of top species to return per plot (default: 10)')
    parser.add_argument('--output', type=str,
                       help='Output CSV file for results')

    args = parser.parse_args()

    # Get plot IDs
    if args.plot_ids:
        plot_ids = [pid.strip() for pid in args.plot_ids.split(',')]
    else:
        try:
            with open(args.file, 'r') as f:
                plot_ids = [line.strip() for line in f if line.strip()]
        except FileNotFoundError:
            print(f"Error: File {args.file} not found")
            sys.exit(1)

    if not plot_ids:
        print("Error: No plot IDs provided")
        sys.exit(1)

    # Validate plot IDs
    valid_plot_ids = []
    for plot_id in plot_ids:
        try:
            parse_plot_id(plot_id)
            valid_plot_ids.append(plot_id)
        except ValueError as e:
            print(f"Warning: Skipping invalid plot ID: {e}")

    if not valid_plot_ids:
        print("Error: No valid plot IDs found")
        sys.exit(1)

    # Look up species data
    results_df = lookup_multiple_plots(valid_plot_ids, args.top_n)

    if results_df is not None:
        # Display summary
        summarize_results(results_df)

        # Export if requested
        if args.output:
            results_df.to_csv(args.output, index=False)
            print(f"\nDetailed results exported to: {args.output}")

        # Show sample data
        print(f"\nSample Results (first 20 rows):")
        print(results_df[['PLOT_ID', 'SPECIES_NAME', 'TREE_COUNT', 'TOTAL_CARBON', 'AVG_DIAMETER']].head(20).to_string(index=False))

    else:
        print("No results found")
        sys.exit(1)

if __name__ == '__main__':
    main()