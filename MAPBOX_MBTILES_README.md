# Tree Habitat Mbtiles Creation

This directory contains scripts to process PNW tree habitat GeoJSON files into Mapbox-compatible mbtiles using tippecanoe.

## Overview

The scripts process 17 tree species habitat files (`*_edited_cached.geojson`) with a total size of ~234MB and create a single optimized mbtiles file for web mapping.

### Species Included

- Bigleaf Maple
- Black Cottonwood
- Douglas Fir  
- Engelmann Spruce
- Grand Fir
- Lodgepole Pine
- Mountain Hemlock
- Noble Fir
- Pacific Madrone
- Pacific Silver Fir
- Ponderosa Pine
- Quaking Aspen
- Red Alder
- Subalpine Fir
- Tanoak
- Western Hemlock
- Western Juniper
- White Fir

## Scripts Available

### 1. Python Script: `create_tree_habitat_mbtiles.py`
- More detailed output and error handling
- Cross-platform compatibility
- Creates detailed metadata JSON

### 2. Bash Script: `create_tree_habitat_mbtiles.sh`
- Faster execution
- Colorized output
- Unix/macOS optimized

## Prerequisites

Install tippecanoe:
```bash
# macOS
brew install tippecanoe

# Ubuntu/Debian
sudo apt-get install tippecanoe

# From source
git clone https://github.com/mapbox/tippecanoe.git
cd tippecanoe && make && sudo make install
```

## Usage

### Python Version
```bash
python3 create_tree_habitat_mbtiles.py
```

### Bash Version
```bash
./create_tree_habitat_mbtiles.sh
```

## Output

Both scripts create:
- **Mbtiles file**: `outputs/mapbox/pnw_tree_habitat_distributions.mbtiles`
- **Metadata**: `outputs/mapbox/pnw_tree_habitat_distributions_metadata.json`

## Tippecanoe Settings

The scripts use optimized settings for habitat polygon data:

- **Zoom levels**: 0-12 (global to detailed regional view)
- **Simplification**: 10 (balance between accuracy and file size)
- **Buffer**: 5 pixels (smooth tile boundaries)
- **Drop densest**: Automatic handling of high-density areas
- **No limits**: Preserves all features and detail

## Layer Names

Each species becomes a separate layer with names like:
- `douglas_fir`
- `ponderosa_pine` 
- `western_hemlock`
- etc.

## File Size Expectations

- **Input**: ~234MB (17 GeoJSON files)
- **Expected output**: ~30-60MB (compressed mbtiles)
- **Compression ratio**: ~4-8x smaller

## Processing Time

Expect 2-10 minutes depending on your system:
- **Fast systems**: 2-4 minutes
- **Average systems**: 4-8 minutes  
- **Slower systems**: 8-15 minutes

## Troubleshooting

### tippecanoe not found
```bash
brew install tippecanoe
# or check PATH: echo $PATH
```

### Memory issues
If processing fails due to memory:
- Close other applications
- Try processing smaller subsets first
- Consider increasing swap space

### Large output file
If the mbtiles file is larger than expected:
- Check if source files have excessive detail
- Consider increasing simplification parameter
- Verify polygon topology is clean

## Using the Mbtiles

### Upload to Mapbox
1. Go to [Mapbox Studio](https://studio.mapbox.com/)
2. Navigate to Tilesets
3. Upload the `.mbtiles` file
4. Use in your map styles

### Local Development
```javascript
// Add to your Mapbox GL JS map
map.addSource('tree-habitat', {
    'type': 'vector',
    'url': 'mapbox://your-username.tileset-id'
});

// Add layers for each species
map.addLayer({
    'id': 'douglas-fir',
    'type': 'fill',
    'source': 'tree-habitat',
    'source-layer': 'douglas_fir',
    'paint': {
        'fill-color': '#2d5016',
        'fill-opacity': 0.7
    }
});
```

## Customization

### Adjust Zoom Levels
Edit the scripts to change `--minimum-zoom` and `--maximum-zoom` values.

### Change Output Location
Modify the `OUTPUT_FILE` variable in either script.

### Add More Species
The scripts automatically detect any `*_edited_cached.geojson` files in the input directory.

## Metadata

The generated metadata includes:
- Processing timestamp
- Input/output file sizes
- Processing time
- List of species and layer names
- Compression statistics 