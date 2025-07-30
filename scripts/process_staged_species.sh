#!/bin/bash

# Process staged species distribution files
echo "🌲 Processing staged species distribution files..."

# Array of species with their codes and names
declare -A species_codes=(
    ["bigleaf_maple"]="312:Bigleaf maple"
    ["engelmann_spruce"]="93:Engelmann spruce"
    ["noble_fir"]="22:Noble fir"
    ["pacific_madrone"]="361:Pacific madrone"
    ["quaking_aspen"]="746:Quaking aspen"
    ["subalpine_fir"]="19:Subalpine fir"
    ["tanoak"]="631:Tanoak"
    ["western_juniper"]="64:Western juniper"
    ["white_fir"]="15:White fir"
)

# Find all staged files
staged_files=(outputs/species_distribution/*_edited.geojson)

if [ ${#staged_files[@]} -eq 0 ]; then
    echo "❌ No staged files found!"
    exit 1
fi

echo "📁 Found ${#staged_files[@]} staged files"

# Process each file
for file in "${staged_files[@]}"; do
    if [[ ! -f "$file" ]]; then
        echo "⚠️  File not found: $file"
        continue
    fi
    
    # Extract species name from filename
    filename=$(basename "$file")
    species_slug="${filename%_merged_habitat_regional_habitat_2km_radius_edited.geojson}"
    
    echo ""
    echo "==============================================="
    echo "🔄 Processing: $species_slug"
    echo "📁 File: $filename"
    
    # Look up species code and name
    if [[ -z "${species_codes[$species_slug]}" ]]; then
        echo "⚠️  Unknown species: $species_slug"
        continue
    fi
    
    IFS=':' read -r species_code species_name <<< "${species_codes[$species_slug]}"
    
    echo "🏷️  Species: $species_name (code: $species_code)"
    echo "==============================================="
    
    # Run the processing command
    echo "🚀 Starting processing..."
    CUSTOM_SPECIES_CODE="$species_code" CUSTOM_SPECIES_NAME="$species_name" \
        python scripts/visualization/create_species_distribution_maps_cached.py \
        --workflow staged --buffer-size medium
    
    exit_code=$?
    
    if [ $exit_code -eq 0 ]; then
        echo "✅ SUCCESS: $species_name"
    else
        echo "❌ FAILED: $species_name (exit code: $exit_code)"
        echo "❓ Continue with next species? (y/n)"
        read -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            echo "🛑 Stopping processing"
            exit 1
        fi
    fi
done

echo ""
echo "🎉 Bulk processing complete!" 