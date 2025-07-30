#!/usr/bin/env python3
"""
Regional TCC Mask Usage Examples

This script demonstrates how to use the enhanced regional TCC mask functionality
for creating species distribution maps with different tree canopy cover requirements.

Usage examples for common scenarios:
"""

import subprocess
import sys
from pathlib import Path

def run_command(cmd, description):
    """Run a command with nice output formatting"""
    print(f"\n🚀 {description}")
    print(f"Command: {' '.join(cmd)}")
    print("-" * 60)
    
    try:
        result = subprocess.run(cmd, check=True)
        print(f"✅ Command completed successfully")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ Command failed with exit code {e.returncode}")
        return False
    except Exception as e:
        print(f"❌ Command failed: {e}")
        return False

def main():
    print("🌲 REGIONAL TCC MASK USAGE EXAMPLES")
    print("="*60)
    print()
    
    examples = [
        {
            'title': "1. List Available TCC Masks",
            'description': "See what regional TCC masks are already available",
            'commands': [
                ["python", "scripts/create_regional_tcc_mask.py", "--list-existing"],
                ["python", "scripts/visualization/create_species_distribution_maps_cached.py", "--list-tcc-masks"]
            ]
        },
        
        {
            'title': "2. Create Common TCC Mask Variants (Batch)",
            'description': "Create all commonly needed TCC mask variants at once",
            'commands': [
                ["python", "scripts/create_regional_tcc_mask.py", "--create-batch"]
            ]
        },
        
        {
            'title': "3. Find Best Available Mask",
            'description': "Find the best existing mask for specific parameters",
            'commands': [
                ["python", "scripts/create_regional_tcc_mask.py", "--find-best-match", "5", "30", "2.0"],
                ["python", "scripts/create_regional_tcc_mask.py", "--find-best-match", "10"],
                ["python", "scripts/create_regional_tcc_mask.py", "--find-best-match", "20", "60"]
            ]
        },
        
        {
            'title': "4. Auto-Create Missing Mask",
            'description': "Automatically create a specific mask if it doesn't exist",
            'commands': [
                ["python", "scripts/create_regional_tcc_mask.py", "--auto-create", "15", "30", "1.0"]
            ]
        },
        
        {
            'title': "5. Species Maps with Intelligent TCC Selection",
            'description': "Create species maps with automatic best-match TCC mask selection",
            'commands': [
                ["python", "scripts/visualization/create_species_distribution_maps_cached.py", 
                 "--species", "ponderosa-pine", "--find-best-tcc-mask"],
                ["python", "scripts/visualization/create_species_distribution_maps_cached.py", 
                 "--species", "douglas-fir", "--auto-create-tcc-mask", "--tcc-threshold", "10"]
            ]
        },
        
        {
            'title': "6. Species Maps with Specific TCC Masks",
            'description': "Use a specific regional TCC mask file",
            'commands': [
                ["python", "scripts/visualization/create_species_distribution_maps_cached.py", 
                 "--species", "ponderosa-pine", 
                 "--regional-tcc-mask", "cache/species_masks/pnw_tcc_mask_5pct_30m_2.0sqkm.geojson"]
            ]
        }
    ]
    
    print("📖 AVAILABLE EXAMPLES:")
    for i, example in enumerate(examples, 1):
        print(f"\n{example['title']}")
        print(f"   {example['description']}")
    
    print(f"\n" + "="*60)
    print("💡 WORKFLOW RECOMMENDATIONS:")
    print()
    
    print("📋 FOR FIRST-TIME SETUP:")
    print("1. Check what masks exist: --list-existing")
    print("2. Create common variants: --create-batch")
    print("3. Use species scripts with --find-best-tcc-mask")
    print()
    
    print("🔧 FOR SPECIFIC REQUIREMENTS:")
    print("1. Find best match: --find-best-match THRESHOLD [RESOLUTION] [MIN_AREA]")
    print("2. Auto-create if needed: --auto-create THRESHOLD RESOLUTION MIN_AREA")
    print("3. Use specific mask: --regional-tcc-mask PATH/TO/MASK.geojson")
    print()
    
    print("🎯 FOR DIFFERENT SPECIES NEEDS:")
    print("• Very permissive (1% threshold): Good for open woodland species")
    print("• Standard (5% threshold): Balanced approach for most forest species")
    print("• Moderate (10% threshold): Forest interior species")
    print("• Strict (20% threshold): Dense forest species only")
    print()
    
    print("📏 RESOLUTION TRADE-OFFS:")
    print("• 30m resolution: Maximum detail, slower processing (~15 min)")
    print("• 60m resolution: Good detail, faster processing (~8 min)")
    print("• 120m resolution: Coarse detail, quick processing (~4 min)")
    print()
    
    # Interactive mode
    while True:
        print("="*60)
        choice = input("Enter example number to run (1-6), 'all' to run all, or 'q' to quit: ").strip().lower()
        
        if choice == 'q':
            break
        elif choice == 'all':
            for example in examples:
                print(f"\n{'='*60}")
                print(f"RUNNING: {example['title']}")
                print(f"{'='*60}")
                for cmd in example['commands']:
                    run_command(cmd, f"Executing: {' '.join(cmd)}")
        elif choice.isdigit() and 1 <= int(choice) <= len(examples):
            example = examples[int(choice) - 1]
            print(f"\n{'='*60}")
            print(f"RUNNING: {example['title']}")
            print(f"{'='*60}")
            for cmd in example['commands']:
                run_command(cmd, f"Executing: {' '.join(cmd)}")
        else:
            print("❌ Invalid choice. Please enter 1-6, 'all', or 'q'")

if __name__ == "__main__":
    main() 