#!/usr/bin/env python3
"""
Convert tree species geojson files to optimized MBTiles for Mapbox
Handles file size optimization and zoom level configuration
"""

import os
import subprocess
import json
from pathlib import Path
import argparse
import shutil

# Configuration for different species and file types
SPECIES_CONFIGS = {
    # For continuous distribution files (large polygons)
    'continuous_distribution': {
        'min_zoom': 4,
        'max_zoom': 12,
        'simplification': 2.0,  # More aggressive simplification for large files
        'layer_name': 'species_range',
        'file_pattern': '*continuous_distribution*.geojson'
    },
    
    # For merged habitat files (medium complexity)
    'merged_habitat': {
        'min_zoom': 6,
        'max_zoom': 14,
        'simplification': 1.0,
        'layer_name': 'species_habitat',
        'file_pattern': '*merged_habitat*edited.geojson'
    },
    
    # For individual buffers (highest detail)
    'individual_buffers': {
        'min_zoom': 8,
        'max_zoom': 16,
        'simplification': 0.5,
        'layer_name': 'species_plots',
        'file_pattern': '*individual_buffers*.geojson'
    },
    
    # For plot points (point data)
    'plots_points': {
        'min_zoom': 6,
        'max_zoom': 16,
        'simplification': 0,  # No simplification for points
        'layer_name': 'species_points',
        'file_pattern': '*plots*.geojson'
    }
}

def check_dependencies():
    """Check if required tools are installed"""
    required_tools = ['tippecanoe']
    
    for tool in required_tools:
        if not shutil.which(tool):
            print(f"❌ {tool} is not installed or not in PATH")
            print(f"Please install {tool}:")
            if tool == 'tippecanoe':
                print("  macOS: brew install tippecanoe")
                print("  Ubuntu: sudo apt-get install tippecanoe")
            return False
    
    return True

def get_file_size_mb(filepath):
    """Get file size in MB"""
    return os.path.getsize(filepath) / 1024 / 1024

def estimate_mbtiles_size(geojson_file, config):
    """Estimate final MBTiles size based on geojson file and configuration"""
    geojson_size_mb = get_file_size_mb(geojson_file)
    
    # Rough estimation factors based on tippecanoe compression and zoom levels
    zoom_factor = (config['max_zoom'] - config['min_zoom'] + 1) / 10
    simplification_factor = 1.0 / (config['simplification'] + 0.1)
    compression_factor = 0.3  # tippecanoe typically achieves ~30% of original size
    
    estimated_size = geojson_size_mb * zoom_factor * simplification_factor * compression_factor
    return estimated_size

def adjust_config_for_size(config, estimated_size_mb, target_size_mb=300):
    """Adjust configuration to meet target file size"""
    if estimated_size_mb <= target_size_mb:
        return config
    
    # Calculate adjustment factor
    adjustment_factor = estimated_size_mb / target_size_mb
    
    # Adjust parameters
    adjusted_config = config.copy()
    
    if adjustment_factor > 2:
        # Significant reduction needed
        adjusted_config['max_zoom'] = min(config['max_zoom'], config['min_zoom'] + 6)
        adjusted_config['simplification'] = max(config['simplification'] * 2, 5.0)
    elif adjustment_factor > 1.5:
        # Moderate reduction needed
        adjusted_config['max_zoom'] = min(config['max_zoom'], config['min_zoom'] + 8)
        adjusted_config['simplification'] = config['simplification'] * 1.5
    else:
        # Minor reduction needed
        adjusted_config['simplification'] = config['simplification'] * 1.2
    
    return adjusted_config

def extract_species_name(filename):
    """Extract species name from filename"""
    # Remove common prefixes and suffixes
    name = filename.replace('_continuous_distribution_pnw_elevation_tcc_constrained_from_edited_cached_optimized', '')
    name = name.replace('_continuous_distribution_pnw_elevation_tcc_constrained_from_edited_cached', '')
    name = name.replace('_merged_habitat_regional_habitat_2km_radius_edited', '')
    name = name.replace('_individual_buffers_regional_habitat_2km_radius', '')
    name = name.replace('_plots_regional_habitat_2km_radius_points', '')
    name = name.replace('.geojson', '')
    
    # Clean up species names
    name = name.replace('_', ' ').replace('-', ' ')
    return name.title()

def create_combined_mbtiles(geojson_files, output_file, config, dry_run=False):
    """Create a single MBTiles file with multiple layers from geojson files"""
    
    if not geojson_files:
        print("❌ No geojson files provided")
        return False
    
    print(f"\n🔄 Creating combined MBTiles with {len(geojson_files)} species layers")
    
    # Build tippecanoe command
    cmd = [
        'tippecanoe',
        '--output', str(output_file),
        '--minimum-zoom', str(config['min_zoom']),
        '--maximum-zoom', str(config['max_zoom']),
        '--simplification', str(config['simplification']),
        '--force',  # Overwrite existing files
        '--no-feature-limit',  # Don't limit features per tile
        '--no-tile-size-limit',  # Don't limit tile size
        '--drop-densest-as-needed',  # Drop features to meet size limits
        '--extend-zooms-if-still-dropping',  # Extend zoom if still dropping
        '--name', 'PNW Species Distribution Combined',
        '--description', 'Combined tree species distribution data for Pacific Northwest',
        '--progress-interval', '10',  # Show progress every 10 seconds
        '--quiet'  # Reduce verbose output except for progress
    ]
    
    # Add each geojson file as a separate layer
    layer_info = []
    total_size_mb = 0
    
    for geojson_file in sorted(geojson_files):
        species_name = extract_species_name(geojson_file.stem)
        layer_name = f"{species_name.lower().replace(' ', '_').replace('-', '_')}_range"
        file_size_mb = get_file_size_mb(geojson_file)
        total_size_mb += file_size_mb
        
        layer_info.append({
            'file': geojson_file,
            'species': species_name,
            'layer_name': layer_name,
            'size_mb': file_size_mb
        })
        
        print(f"  📄 {geojson_file.name}")
        print(f"      Species: {species_name}")
        print(f"      Layer: {layer_name}")
        print(f"      Size: {file_size_mb:.1f} MB")
        
        cmd.extend(['--layer', f"{layer_name}:{geojson_file}"])
    
    print(f"\nTotal input size: {total_size_mb:.1f} MB")
    
    # Estimate output size
    estimated_size = total_size_mb * 0.3  # Rough compression estimate
    print(f"Estimated output size: {estimated_size:.1f} MB")
    
    # Check if this might be too large to process efficiently
    if total_size_mb > 500:
        print(f"⚠️  Warning: Large input size ({total_size_mb:.1f} MB) may take a very long time to process!")
        print(f"   Consider processing smaller batches or reducing file sizes first")
        
        response = input("Continue anyway? (y/N): ").lower().strip()
        if response not in ['y', 'yes']:
            print("Aborted by user")
            return False
    
    if estimated_size > 300:
        print(f"⚠️  Warning: Estimated size exceeds 300MB limit!")
        print(f"   Consider reducing max zoom or increasing simplification")
    
    if dry_run:
        print(f"\n🔍 Dry run - would create: {output_file}")
        print(f"Command: {' '.join(cmd[:10])} ... [+{len(cmd)-10} more args]")
        return True
    
    print(f"\n🔄 Creating combined tileset: {output_file.name}")
    print(f"⏳ This may take several minutes for large datasets...")
    print(f"💡 Tip: Use Ctrl+C to cancel if it takes too long")
    
    try:
        # Run with real-time output instead of capturing
        print(f"\n📝 Tippecanoe output:")
        result = subprocess.run(cmd, text=True, check=True)
        
        if output_file.exists():
            output_size_mb = get_file_size_mb(output_file)
            print(f"✅ Combined tileset created successfully!")
            print(f"   Output size: {output_size_mb:.1f} MB")
            print(f"   Compression ratio: {(total_size_mb/output_size_mb):.1f}x")
            print(f"   Layers: {len(layer_info)}")
            
            print(f"\n📋 Layer names for Mapbox GL JS:")
            for info in layer_info:
                print(f"   '{info['layer_name']}' - {info['species']}")
            
            return True
        else:
            print(f"❌ Combined tileset was not created")
            return False
            
    except subprocess.CalledProcessError as e:
        print(f"❌ Error creating combined MBTiles: {e}")
        print(f"   Command: {' '.join(cmd[:10])} ... [truncated]")
        if e.stdout:
            print(f"   Stdout: {e.stdout}")
        if e.stderr:
            print(f"   Stderr: {e.stderr}")
        return False



def auto_detect_file_type(filename):
    """Auto-detect file type based on filename patterns"""
    filename_lower = filename.lower()
    
    for config_type, config in SPECIES_CONFIGS.items():
        pattern = config['file_pattern'].replace('*', '').replace('.geojson', '')
        if pattern in filename_lower:
            return config_type
    
    # Default fallback
    if 'continuous' in filename_lower:
        return 'continuous_distribution'
    elif 'merged' in filename_lower or 'habitat' in filename_lower:
        return 'merged_habitat'
    elif 'buffer' in filename_lower:
        return 'individual_buffers'
    elif 'plot' in filename_lower:
        return 'plots_points'
    else:
        return 'merged_habitat'  # Default

def create_usage_examples(output_dir):
    """Create usage examples for Mapbox"""
    examples = {
        "mapbox_gl_js_examples": {
            "basic_species_layer": {
                "description": "Basic species distribution layer",
                "code": {
                    "source": {
                        "type": "vector",
                        "url": "mapbox://your-username.your-tileset-id"
                    },
                    "layer": {
                        "id": "douglas-fir-distribution",
                        "type": "fill",
                        "source": "douglas-fir-source",
                        "source-layer": "douglas_fir_range",
                        "paint": {
                            "fill-color": "#2E8B57",
                            "fill-opacity": 0.6,
                            "fill-outline-color": "#1C5F3A"
                        }
                    }
                }
            },
            
            "multiple_species_layers": {
                "description": "Add multiple species as separate layers with unique names",
                "examples": [
                    {
                        "species": "Douglas Fir",
                        "source_layer": "douglas_fir_range",
                        "color": "#2E8B57"
                    },
                    {
                        "species": "Western Hemlock", 
                        "source_layer": "western_hemlock_range",
                        "color": "#228B22"
                    },
                    {
                        "species": "Bigleaf Maple",
                        "source_layer": "bigleaf_maple_range", 
                        "color": "#FF6347"
                    }
                ],
                "note": "Each species gets its own unique source-layer name based on the species name"
            },
            
            "zoom_based_visibility": {
                "description": "Show/hide based on zoom level",
                "code": {
                    "layout": {
                        "visibility": "visible"
                    },
                    "paint": {
                        "fill-opacity": [
                            "interpolate", ["linear"], ["zoom"],
                            6, 0.3,
                            10, 0.7,
                            14, 0.9
                        ]
                    }
                }
            }
        },
        
        "mapbox_studio_examples": {
            "species_habitat_style": {
                "filter": ["==", ["get", "habitat_type"], "primary"],
                "paint": {
                    "fill-color": "#2E8B57",
                    "fill-opacity": ["interpolate", ["linear"], ["zoom"], 6, 0.4, 12, 0.8]
                }
            }
        },
        

        
        "combined_tileset_usage": {
            "description": "Using the combined multi-layer tileset",
            "single_source_multiple_layers": {
                "source": {
                    "type": "vector",
                    "url": "mapbox://your-username.your-combined-tileset-id"
                },
                "layers": [
                    {
                        "id": "douglas-fir-layer",
                        "source": "pnw-species-source",
                        "source-layer": "douglas_fir_range",
                        "type": "fill",
                        "paint": {"fill-color": "#2E8B57", "fill-opacity": 0.7}
                    },
                    {
                        "id": "bigleaf-maple-layer",
                        "source": "pnw-species-source", 
                        "source-layer": "bigleaf_maple_range",
                        "type": "fill",
                        "paint": {"fill-color": "#FF6347", "fill-opacity": 0.7}
                    }
                ]
            },
            "layer_control": {
                "toggle_visibility": "map.setLayoutProperty('douglas-fir-layer', 'visibility', 'visible')",
                "change_opacity": "map.setPaintProperty('douglas-fir-layer', 'fill-opacity', 0.5)"
            }
        },
        
        "upload_instructions": {
            "steps": [
                "1. Go to Mapbox Studio > Tilesets",
                "2. Click 'New tileset'", 
                "3. Upload the single combined .mbtiles file",
                "4. One tileset contains all species as separate layers",
                "5. Reference unique layer names in your map styles",
                "6. Toggle each species layer independently"
            ],
            "tips": [
                "Single tileset upload is faster than multiple uploads",
                "Each species gets its own unique layer name",
                "Use zoom levels 6-10 for regional viewing",
                "Higher zoom levels (12-16) for detailed inspection",
                "Perfect balance of organization and control"
            ]
        }
    }
    
    output_file = Path(output_dir) / "mapbox_usage_examples.json"
    with open(output_file, 'w') as f:
        json.dump(examples, f, indent=2)
    
    print(f"📝 Created usage examples: {output_file}")

def is_valid_file_pattern(filename):
    """Check if file matches the required pattern"""
    return (filename.endswith('tcc_constrained_from_edited_cached.geojson') or 
            filename.endswith('tcc_constrained_from_edited_cached_optimized.geojson'))

def create_individual_mbtiles(geojson_file, output_dir, config, dry_run=False):
    """Create individual mbtiles file for a single species"""
    
    geojson_path = Path(geojson_file)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Extract species info
    species_name = extract_species_name(geojson_path.stem)
    layer_name = f"{species_name.lower().replace(' ', '_').replace('-', '_')}_range"
    
    # Create output filename
    output_filename = f"{species_name.lower().replace(' ', '_').replace('-', '_')}.mbtiles"
    output_path = output_dir / output_filename
    
    file_size_mb = get_file_size_mb(geojson_path)
    
    print(f"\n📁 Processing: {species_name}")
    print(f"   Input: {geojson_path.name}")
    print(f"   Size: {file_size_mb:.1f} MB")
    print(f"   Layer: {layer_name}")
    
    if dry_run:
        print(f"   🔍 Would create: {output_path}")
        return True, output_path
    
    # Build tippecanoe command
    cmd = [
        'tippecanoe',
        '--output', str(output_path),
        '--layer', layer_name,
        '--minimum-zoom', str(config['min_zoom']),
        '--maximum-zoom', str(config['max_zoom']),
        '--simplification', str(config['simplification']),
        '--force',
        '--no-feature-limit',
        '--no-tile-size-limit',
        '--drop-densest-as-needed',
        '--extend-zooms-if-still-dropping',
        '--name', f"{species_name} Distribution",
        '--description', f"Species distribution data for {species_name} in Pacific Northwest",
        '--quiet',
        str(geojson_path)
    ]
    
    try:
        result = subprocess.run(cmd, text=True, check=True)
        
        if output_path.exists():
            output_size_mb = get_file_size_mb(output_path)
            print(f"   ✅ Created: {output_size_mb:.1f} MB")
            return True, output_path
        else:
            print(f"   ❌ Failed to create")
            return False, None
            
    except subprocess.CalledProcessError as e:
        print(f"   ❌ Error: {e}")
        return False, None

def combine_individual_mbtiles(mbtiles_dir, output_file, dry_run=False):
    """Combine individual mbtiles files using tile-join"""
    
    mbtiles_dir = Path(mbtiles_dir)
    mbtiles_files = list(mbtiles_dir.glob('*.mbtiles'))
    
    # Filter out combined files
    mbtiles_files = [f for f in mbtiles_files if 'combined' not in f.name.lower()]
    
    if not mbtiles_files:
        print(f"❌ No individual mbtiles files found in {mbtiles_dir}")
        return False
    
    print(f"\n🔗 Combining {len(mbtiles_files)} individual mbtiles files")
    
    total_size_mb = 0
    for mbtiles_file in sorted(mbtiles_files):
        size_mb = get_file_size_mb(mbtiles_file)
        total_size_mb += size_mb
        print(f"  📄 {mbtiles_file.name} ({size_mb:.1f} MB)")
    
    print(f"\nTotal mbtiles size: {total_size_mb:.1f} MB")
    
    if dry_run:
        print(f"🔍 Would combine into: {output_file}")
        return True
    
    # Build tile-join command
    cmd = ['tile-join', '-f', '-o', str(output_file)]
    cmd.extend([str(f) for f in sorted(mbtiles_files)])
    
    print(f"\n🔄 Combining into: {output_file.name}")
    
    try:
        result = subprocess.run(cmd, text=True, check=True)
        
        if output_file.exists():
            output_size_mb = get_file_size_mb(output_file)
            print(f"✅ Combined tileset created: {output_size_mb:.1f} MB")
            return True
        else:
            print(f"❌ Combined tileset was not created")
            return False
            
    except subprocess.CalledProcessError as e:
        print(f"❌ Error combining: {e}")
        if e.stderr:
            print(f"   Error details: {e.stderr}")
        return False

def main():
    parser = argparse.ArgumentParser(description='Convert species geojson files to Mapbox MBTiles')
    parser.add_argument('input', help='Input geojson file or directory')
    parser.add_argument('-o', '--output', default='outputs/mapbox_mbtiles/pnw_tree_species_combined.mbtiles', 
                       help='Output mbtiles file path (default: outputs/mapbox_mbtiles/pnw_species_combined.mbtiles)')
    parser.add_argument('-t', '--type', choices=list(SPECIES_CONFIGS.keys()),
                       help='File type configuration (auto-detected if not specified)', default='continuous_distribution')
    parser.add_argument('--dry-run', action='store_true',
                       help='Show what would be processed without creating files')
    parser.add_argument('--max-size', type=int, default=300,
                       help='Maximum file size in MB (default: 300)')
    parser.add_argument('--list-configs', action='store_true',
                       help='List available configuration types')
    parser.add_argument('--fast-mode', action='store_true',
                       help='Use faster settings for large datasets (lower quality but faster processing)')
    parser.add_argument('--ultra-fast', action='store_true',
                       help='Use ultra-fast settings for testing (very low quality, very fast)')
    parser.add_argument('--two-step', action='store_true',
                       help='Process in two steps: 1) Create individual mbtiles, 2) Combine them (faster and more reliable)')
    parser.add_argument('--step', choices=['individual', 'combine'], 
                       help='Run only specific step of two-step process')
    
    args = parser.parse_args()
    
    if args.list_configs:
        print("Available configuration types:")
        for config_type, config in SPECIES_CONFIGS.items():
            print(f"  {config_type}:")
            print(f"    Pattern: {config['file_pattern']}")
            print(f"    Zoom: {config['min_zoom']}-{config['max_zoom']}")
            print(f"    Layer: {config['layer_name']}")
        return
    
    # Check dependencies
    if not args.dry_run:
        need_tile_join = args.two_step or args.step == 'combine'
        if need_tile_join:
            if not check_dependencies() or not shutil.which('tile-join'):
                print("❌ tile-join is required for two-step process")
                print("Install with: brew install tippecanoe (tile-join is included)")
                return
        else:
            if not check_dependencies():
                print("Please install missing dependencies and try again.")
                return
    
    input_path = Path(args.input)
    output_path = Path(args.output)
    
    # Create output directory
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Get configuration
    config = SPECIES_CONFIGS[args.type].copy()
    
    # Apply fast mode optimizations if requested
    if args.ultra_fast:
        print("⚡ Ultra-fast mode enabled - minimum quality for maximum speed")
        config['max_zoom'] = 6  # Very low max zoom
        config['min_zoom'] = 4  # Keep min zoom reasonable
        config['simplification'] = 20.0  # Very aggressive simplification
        print(f"   Adjusted zoom: {config['min_zoom']}-{config['max_zoom']}")
        print(f"   Adjusted simplification: {config['simplification']}")
    elif args.fast_mode:
        print("🚀 Fast mode enabled - optimizing for speed over quality")
        config['max_zoom'] = min(config['max_zoom'], 8)  # Even lower max zoom
        config['simplification'] = max(config['simplification'], 10.0)  # Even more aggressive simplification
        print(f"   Adjusted zoom: {config['min_zoom']}-{config['max_zoom']}")
        print(f"   Adjusted simplification: {config['simplification']}")
    
    # Process files
    if input_path.is_file():
        # Single file - check if it matches the required pattern
        if not is_valid_file_pattern(input_path.name):
            print(f"❌ Skipping {input_path.name} - does not match required pattern '*tcc_constrained_from_edited_cached.geojson'")
            return
        
        # Create combined mbtiles with single file
        geojson_files = [input_path]
        
    elif input_path.is_dir():
        # Directory of files - filter to only valid pattern files
        all_geojson_files = list(input_path.glob('*.geojson'))
        geojson_files = [f for f in all_geojson_files if is_valid_file_pattern(f.name)]
        
        if not all_geojson_files:
            print(f"No geojson files found in {input_path}")
            return
        
        if not geojson_files:
            print(f"Found {len(all_geojson_files)} geojson files, but none match the required patterns:")
            print("  - '*tcc_constrained_from_edited_cached.geojson'")
            print("  - '*tcc_constrained_from_edited_cached_optimized.geojson'")
            print("Skipped files:")
            for f in all_geojson_files:
                print(f"  - {f.name}")
            return
        
        print(f"Found {len(all_geojson_files)} geojson files, processing {len(geojson_files)} that match required patterns")
        
        if len(geojson_files) < len(all_geojson_files):
            skipped_files = [f for f in all_geojson_files if not is_valid_file_pattern(f.name)]
            print(f"Skipped {len(skipped_files)} files that don't match pattern:")
            for f in skipped_files:
                print(f"  - {f.name}")
            print()
    else:
        print(f"Input path does not exist: {input_path}")
        return
    
    # Handle two-step process
    if args.two_step or args.step:
        individual_dir = output_path.parent / "individual"
        
        # Step 1: Create individual mbtiles files
        if args.step != 'combine':
            print("🚀 Step 1: Creating individual mbtiles files")
            print(f"Individual files will be saved to: {individual_dir}")
            
            successes = 0
            failed_files = []
            
            for geojson_file in geojson_files:
                success, mbtiles_path = create_individual_mbtiles(geojson_file, individual_dir, config, args.dry_run)
                if success:
                    successes += 1
                else:
                    failed_files.append(geojson_file.name)
            
            print(f"\n📊 Step 1 Results: {successes}/{len(geojson_files)} files processed successfully")
            if failed_files:
                print(f"❌ Failed files: {', '.join(failed_files)}")
            
            if args.step == 'individual':
                print("✅ Individual processing complete!")
                return
        
        # Step 2: Combine individual mbtiles files
        if args.step != 'individual':
            print("\n🔗 Step 2: Combining individual mbtiles files")
            
            if not individual_dir.exists():
                print(f"❌ Individual mbtiles directory not found: {individual_dir}")
                print("Run step 1 first or check the directory path")
                return
            
            success = combine_individual_mbtiles(individual_dir, output_path, args.dry_run)
            
            if not success:
                print("❌ Failed to combine mbtiles files")
                exit(1)
    else:
        # Original single-step process
        success = create_combined_mbtiles(geojson_files, output_path, config, args.dry_run)
        
        if not success:
            print("❌ Failed to create combined mbtiles file")
            exit(1)
    
    # Create usage examples
    if not args.dry_run:
        create_usage_examples(output_path.parent)
        
        print(f"\n📋 Next steps:")
        print(f"1. Upload combined mbtiles file to Mapbox Studio: {output_path}")
        print(f"2. All species will be available as separate layers in one tileset")
        print(f"3. Use unique layer names in your map styles (see above output)")
        print(f"4. See mapbox_usage_examples.json for code examples")

if __name__ == "__main__":
    main() 