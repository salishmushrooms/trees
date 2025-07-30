"""
QGIS Python Script to sum carbon values for all species in selected features.
Run this in the QGIS Python Console.
"""

# Get the active layer
layer = iface.activeLayer()

if not layer:
    print("No active layer selected!")
else:
    # Get selected features
    selected = layer.selectedFeatures()
    
    if not selected:
        print("No features selected!")
    else:
        # Find all carbon columns
        carbon_fields = [field.name() for field in layer.fields() 
                        if field.name().endswith('_AG_CARBON') 
                        and not field.name().endswith('_x') 
                        and not field.name().endswith('_y')]
        
        # Calculate sums
        sums = {}
        for field in carbon_fields:
            total = sum(f[field] for f in selected if f[field] is not None)
            if total > 0:  # Only show species with carbon
                species = field.replace('_AG_CARBON', '').replace('_', ' ').title()
                sums[species] = total
        
        # Display results
        print(f"\nCarbon Summary for {len(selected)} selected plots:")
        print("-" * 50)
        
        # Sort by carbon amount
        for species, carbon in sorted(sums.items(), key=lambda x: x[1], reverse=True):
            print(f"{species:.<30} {carbon:,.0f} units")
        
        print("-" * 50)
        print(f"TOTAL CARBON:.................  {sum(sums.values()):,.0f} units")
        
        # Optional: Copy summary to clipboard
        import json
        QApplication.clipboard().setText(json.dumps(sums, indent=2))