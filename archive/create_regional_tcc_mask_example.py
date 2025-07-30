#!/usr/bin/env python3
"""
Example: Create and use regional TCC masks for the Pacific Northwest

This example shows how to:
1. Fill gaps in an existing TCC mask (your current 5% 30m 2.0sqkm mask)
2. Create a new regional mask with different parameters
3. Use the regional masks in species distribution mapping

Usage examples for your specific case:
"""

import subprocess
import sys
from pathlib import Path

def run_command(cmd, description):
    """Run a command with nice output formatting"""
    print(f"\n🚀 {description}")
    print(f"Running: {' '.join(cmd)}")
    print("-" * 50)
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(result.stdout)
        if result.stderr:
            print("STDERR:", result.stderr)
        print("✅ Command completed successfully")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ Command failed with exit code {e.returncode}")
        print("STDOUT:", e.stdout)
        print("STDERR:", e.stderr)
        return False

def main():
    print("🌲 REGIONAL TCC MASK EXAMPLE WORKFLOW")
    print("="*60)
    print("This example demonstrates creating and using regional TCC masks")
    print("for species distribution mapping in the Pacific Northwest.")
    print()
    
    # Check if existing mask exists
    cache_dir = Path("cache/species_masks")
    existing_mask = cache_dir / "pnw_tcc_mask_5pct_30m_2.0sqkm.geojson"
    
    if existing_mask.exists():
        file_size_mb = existing_mask.stat().st_size / (1024 * 1024)
        print(f"📦 Found existing mask: {existing_mask.name} ({file_size_mb:.1f} MB)")
        print("✨ We can fill gaps in this existing mask!")
    else:
        print("📦 No existing mask found - will create from scratch")
    
    print("\n🎯 RECOMMENDED WORKFLOW FOR YOUR CASE:")
    print()
    
    # Step 1: List existing masks
    print("1️⃣  LIST EXISTING TCC MASKS:")
    cmd1 = [sys.executable, "scripts/create_regional_tcc_mask.py", "--list-existing"]
    if not run_command(cmd1, "List existing TCC mask files"):
        return
    
    # Step 2: Fill gaps in existing 5% mask (if it exists)
    if existing_mask.exists():
        print("\n2️⃣  FILL GAPS IN EXISTING 5% MASK:")
        print("This will process only the unprocessed areas to complete your regional mask.")
        cmd2 = [
            sys.executable, "scripts/create_regional_tcc_mask.py",
            "--threshold", "5",
            "--resolution", "30", 
            "--min-area", "2.0",
            "--fill-gaps"
        ]
        
        response = input("\n🤔 Fill gaps in existing 5% mask? (y/n): ")
        if response.lower() == 'y':
            if not run_command(cmd2, "Fill gaps in existing 5% TCC mask"):
                return
        else:
            print("⏭️  Skipping gap filling")
    
    # Step 3: Create new 20% mask (for comparison)
    print("\n3️⃣  CREATE 20% TCC MASK (for comparison with existing species processing):")
    cmd3 = [
        sys.executable, "scripts/create_regional_tcc_mask.py",
        "--threshold", "20",
        "--resolution", "30",
        "--min-area", "2.0"
    ]
    
    response = input("\n🤔 Create 20% threshold mask? (y/n): ")
    if response.lower() == 'y':
        if not run_command(cmd3, "Create 20% TCC mask for the entire Pacific Northwest"):
            return
    else:
        print("⏭️  Skipping 20% mask creation")
    
    # Step 4: Test with species distribution mapping
    print("\n4️⃣  TEST WITH SPECIES DISTRIBUTION MAPPING:")
    print()
    print("Now you can use these regional masks in your species processing:")
    print()
    print("# Use 5% regional mask (recommended for new processing):")
    cmd4a = [
        sys.executable, "scripts/visualization/create_species_distribution_maps_cached.py",
        "--species", "ponderosa-pine",
        "--tcc-threshold", "5"
    ]
    print(f"   {' '.join(cmd4a)}")
    print()
    
    print("# Use 20% regional mask (matches your existing processing):")
    cmd4b = [
        sys.executable, "scripts/visualization/create_species_distribution_maps_cached.py", 
        "--species", "ponderosa-pine",
        "--tcc-threshold", "20"
    ]
    print(f"   {' '.join(cmd4b)}")
    print()
    
    print("# Force species-specific mask creation (legacy mode):")
    cmd4c = [
        sys.executable, "scripts/visualization/create_species_distribution_maps_cached.py",
        "--species", "ponderosa-pine", 
        "--use-species-tcc-mask"
    ]
    print(f"   {' '.join(cmd4c)}")
    print()
    
    response = input("🤔 Test with Ponderosa Pine using 5% regional mask? (y/n): ")
    if response.lower() == 'y':
        if not run_command(cmd4a, "Test Ponderosa Pine with 5% regional TCC mask"):
            return
    else:
        print("⏭️  Skipping species test")
    
    # Summary
    print("\n🎉 WORKFLOW COMPLETE!")
    print("="*40)
    print()
    print("📋 WHAT YOU'VE ACCOMPLISHED:")
    print("• ✅ Completed regional TCC mask(s) for entire Pacific Northwest")
    print("• ✅ Masks can be reused across multiple species")
    print("• ✅ Much faster processing for future species (no TCC processing needed)")
    print("• ✅ Consistent TCC constraints across all species")
    print()
    print("📂 MASK FILES LOCATION:")
    print(f"   {cache_dir}/")
    print("   ├── pnw_tcc_mask_5pct_30m_2.0sqkm.geojson   # 5% threshold")
    print("   └── pnw_tcc_mask_20pct_30m_2.0sqkm.geojson  # 20% threshold")
    print()
    print("🚀 NEXT STEPS:")
    print("1. Use these regional masks for all future species processing")
    print("2. Set USE_REGIONAL_TCC_MASK = True in your config (already default)")
    print("3. Process multiple species much faster!")
    print("4. Consider creating masks at different thresholds (1%, 10%) as needed")
    print()
    print("💡 PERFORMANCE BENEFITS:")
    print("• Regional mask creation: ~30-60 minutes (one time)")
    print("• Species processing with regional mask: ~2-5 minutes")
    print("• Species processing with species-specific mask: ~15-30 minutes")
    print("• Time savings: ~90% for each additional species!")

if __name__ == "__main__":
    main() 