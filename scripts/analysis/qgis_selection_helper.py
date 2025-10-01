"""
QGIS Selection Helper

This script can be run from the QGIS Python console to export selected
PLOT_IDs to a text file for use with the species lookup tool.

Instructions for QGIS:
1. Load fia_plots_for_qgis.geojson
2. Select features using any selection tool
3. Open QGIS Python console (Ctrl+Alt+P)
4. Copy and paste this script into the console
5. Run the script
6. It will save selected PLOT_IDs to a text file

Alternative: Save as .py file and run with execfile() in QGIS console
"""

# Get the active layer (should be the plot layer)
layer = iface.activeLayer()

if layer is None:
    print("No active layer! Please select the plot layer first.")
else:
    # Get selected features
    selected_features = layer.selectedFeatures()

    if len(selected_features) == 0:
        print("No features selected! Please select some plots first.")
    else:
        # Extract PLOT_IDs
        plot_ids = []
        for feature in selected_features:
            plot_id = feature['PLOT_ID']
            if plot_id:
                plot_ids.append(plot_id)

        if plot_ids:
            # Save to file
            output_file = '/Users/JJC/trees/outputs/selected_plot_ids.txt'
            with open(output_file, 'w') as f:
                for plot_id in plot_ids:
                    f.write(f"{plot_id}\n")

            print(f"Exported {len(plot_ids)} PLOT_IDs to: {output_file}")
            print("Sample PLOT_IDs:")
            for i, plot_id in enumerate(plot_ids[:5]):
                print(f"  {i+1}. {plot_id}")

            if len(plot_ids) > 5:
                print(f"  ... and {len(plot_ids) - 5} more")

            print(f"\nNext step:")
            print(f"python species_lookup_tool.py --file {output_file} --output species_results.csv")
        else:
            print("No valid PLOT_IDs found in selected features!")