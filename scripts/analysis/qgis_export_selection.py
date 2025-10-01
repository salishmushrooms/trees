# QGIS Python Console Script
# Copy and paste this entire script into the QGIS Python console
# Make sure you have selected some features first!

# Get the active layer
layer = iface.activeLayer()

if layer is None:
    print("ERROR: No active layer! Please click on the plot layer first.")
elif layer.selectedFeatureCount() == 0:
    print("ERROR: No features selected! Please select some plots first using any selection tool.")
else:
    # Get selected features
    selected_features = layer.selectedFeatures()
    plot_ids = []

    # Extract PLOT_IDs from selected features
    for feature in selected_features:
        plot_id = feature['PLOT_ID']
        if plot_id:
            plot_ids.append(str(plot_id))

    if plot_ids:
        # Save to file
        output_file = '/Users/JJC/trees/outputs/selected_plot_ids.txt'
        with open(output_file, 'w') as f:
            for plot_id in plot_ids:
                f.write(f"{plot_id}\n")

        print(f"SUCCESS: Exported {len(plot_ids)} PLOT_IDs to: {output_file}")
        print("\nSample PLOT_IDs:")
        for i, plot_id in enumerate(plot_ids[:5]):
            print(f"  {i+1}. {plot_id}")

        if len(plot_ids) > 5:
            print(f"  ... and {len(plot_ids) - 5} more")

        print(f"\nNext step (run in terminal):")
        print(f"cd /Users/JJC/trees")
        print(f"python scripts/analysis/species_lookup_tool.py --file outputs/selected_plot_ids.txt --output outputs/species_results.csv")
    else:
        print("ERROR: No valid PLOT_IDs found in selected features!")