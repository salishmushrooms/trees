#!/bin/bash

# Script to create Mapbox-ready vector tilesets from tree species range GeoJSON files
# Creates a single tileset with each species as a separate layer for easy toggle functionality

set -e

echo "🌲 Creating Mapbox Vector Tileset for Tree Species Ranges"
echo "========================================================="

# Configuration
INPUT_DIR="outputs/species_distribution"
OUTPUT_DIR="outputs/mapbox"
TILESET_NAME="pnw_tree_species_ranges"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Check if tippecanoe is installed
if ! command -v tippecanoe &> /dev/null; then
    echo "❌ Error: tippecanoe is not installed"
    echo "Install with: brew install tippecanoe (macOS) or see https://github.com/mapbox/tippecanoe"
    exit 1
fi

echo "📁 Found range files:"
find "$INPUT_DIR" -name "*-range.geojson" -type f | sort

# Verify we have range files
range_files=($(find "$INPUT_DIR" -name "*-range.geojson" -type f))
if [ ${#range_files[@]} -eq 0 ]; then
    echo "❌ Error: No *-range.geojson files found in $INPUT_DIR"
    exit 1
fi

echo ""
echo "🔄 Creating tileset with individual species layers..."

# Build tippecanoe arguments
TIPPECANOE_ARGS=(
    --output="$OUTPUT_DIR/${TILESET_NAME}.mbtiles"
    --force
    --maximum-zoom=12
    --minimum-zoom=0
    --drop-densest-as-needed
    --extend-zooms-if-still-dropping
    --buffer=5
    --simplification=10
    --detect-shared-borders
    --coalesce-densest-as-needed
    --attribution="Forest Inventory and Analysis (FIA) Database"
)

# Add each range file as a separate layer
species_count=0
for file in "${range_files[@]}"; do
    if [[ -f "$file" ]]; then
        # Extract species name from filename (e.g., douglas-fir-range.geojson -> douglas-fir)
        species_name=$(basename "$file" | sed 's/-range\.geojson$//')
        echo "  📊 Adding layer: $species_name"
        TIPPECANOE_ARGS+=("--layer=$species_name:$file")
        ((species_count++))
    fi
done

echo ""
echo "🚀 Running tippecanoe with $species_count species layers..."
tippecanoe "${TIPPECANOE_ARGS[@]}"

echo "✅ Created: $OUTPUT_DIR/${TILESET_NAME}.mbtiles"

# Display tileset information
echo ""
echo "📊 Tileset Information:"
echo "======================="
echo "File: $OUTPUT_DIR/${TILESET_NAME}.mbtiles"
ls -lh "$OUTPUT_DIR/${TILESET_NAME}.mbtiles"

# List layers if tile-join is available
if command -v tile-join &> /dev/null; then
    echo ""
    echo "Layers in tileset:"
    tile-join --list-layers "$OUTPUT_DIR/${TILESET_NAME}.mbtiles" 2>/dev/null | sed 's/^/  - /' || echo "  (Layer listing not available)"
fi

echo ""
echo "🎯 Next Steps:"
echo "=============="
echo "1. Upload the .mbtiles file to your Mapbox account:"
echo "   - Go to https://studio.mapbox.com/"
echo "   - Navigate to Tilesets section"
echo "   - Click 'New tileset' and upload: $OUTPUT_DIR/${TILESET_NAME}.mbtiles"
echo ""
echo "2. Get your tileset URL (will be something like):"
echo "   mapbox://your-username.${TILESET_NAME}"
echo ""
echo "3. In your web application:"
echo "   - Each species will be a separate layer you can toggle"
echo "   - Layer names match the species identifiers:"
for file in "${range_files[@]}"; do
    if [[ -f "$file" ]]; then
        species_name=$(basename "$file" | sed 's/-range\.geojson$//')
        echo "     * $species_name"
    fi
done
echo ""
echo "4. Use Mapbox GL JS to add layers like:"
echo "   map.addLayer({"
echo "     id: 'douglas-fir-fill',"
echo "     type: 'fill',"
echo "     source: 'species-ranges',"
echo "     'source-layer': 'douglas-fir',"
echo "     paint: { 'fill-color': '#228B22', 'fill-opacity': 0.6 }"
echo "   });"
echo ""
echo "💡 For a complete web implementation example, see:"
echo "   https://docs.mapbox.com/mapbox-gl-js/example/toggle-layers/" 