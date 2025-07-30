#!/usr/bin/env python3
"""
Bulk Process Staged Species Distribution Files

This script processes all *_edited.geojson files in the outputs/species_distribution directory
by running the visualization script for each staged species.

It extracts species names from filenames and maps them to appropriate species codes.
"""

import subprocess
import sys
from pathlib import Path
import time
import os
import re

# Species name to code mapping (based on FIA database codes)
SPECIES_MAPPING = {
    "bigleaf_maple": {"code": 312, "name": "Bigleaf maple"},
    "engelmann_spruce": {"code": 93, "name": "Engelmann spruce"},
    "noble_fir": {"code": 22, "name": "Noble fir"},
    "pacific_madrone": {"code": 361, "name": "Pacific madrone"},
    "quaking_aspen": {"code": 746, "name": "Quaking aspen"},
    "subalpine_fir": {"code": 19, "name": "Subalpine fir"},
    "tanoak": {"code": 631, "name": "Tanoak"},
    "western_juniper": {"code": 64, "name": "Western juniper"},
    "white_fir": {"code": 15, "name": "White fir"}
}

def find_staged_files():
    """Find all staged *_edited.geojson files"""
    species_dir = Path("outputs/species_distribution")
    staged_files = list(species_dir.glob("*_edited.geojson"))
    return staged_files

def extract_species_slug(filename):
    """Extract species slug from filename"""
    # Pattern: species_name_merged_habitat_regional_habitat_2km_radius_edited.geojson
    pattern = r"(.+)_merged_habitat_regional_habitat_2km_radius_edited\.geojson"
    match = re.match(pattern, filename.name)
    if match:
        return match.group(1)
    return None

def run_species_processing(species_slug, species_code, species_name):
    """Run the visualization script for a specific species"""
    print(f"\n{'='*60}")
    print(f"🌲 PROCESSING: {species_name} (code: {species_code})")
    print(f"📁 Species slug: {species_slug}")
    print(f"{'='*60}")
    
    # Set environment variables
    env = os.environ.copy()
    env['CUSTOM_SPECIES_CODE'] = str(species_code)
    env['CUSTOM_SPECIES_NAME'] = species_name
    
    # Run the command
    cmd = [
        sys.executable,
        "scripts/visualization/create_species_distribution_maps_cached.py",
        "--workflow", "staged",
        "--buffer-size", "medium"
    ]
    
    try:
        start_time = time.time()
        
        print(f"🚀 Starting processing...")
        print(f"📋 Command: CUSTOM_SPECIES_CODE={species_code} CUSTOM_SPECIES_NAME=\"{species_name}\" python {' '.join(cmd[1:])}")
        
        result = subprocess.run(
            cmd,
            env=env,
            capture_output=True,
            text=True,
            timeout=1800  # 30 minute timeout
        )
        
        elapsed_time = time.time() - start_time
        
        if result.returncode == 0:
            print(f"✅ SUCCESS! Completed in {elapsed_time:.1f} seconds")
            if result.stdout:
                print("📤 Output:")
                print(result.stdout[-1000:])  # Last 1000 chars
        else:
            print(f"❌ FAILED! (exit code: {result.returncode})")
            if result.stderr:
                print("🔥 Error output:")
                print(result.stderr[-1000:])  # Last 1000 chars
            return False
            
    except subprocess.TimeoutExpired:
        print(f"⏰ TIMEOUT! Processing took longer than 30 minutes")
        return False
    except Exception as e:
        print(f"💥 EXCEPTION: {e}")
        return False
    
    return True

def main():
    """Main processing function"""
    print("🌲 BULK SPECIES DISTRIBUTION PROCESSING")
    print("="*60)
    print("🎯 Processing staged *_edited.geojson files")
    print("📋 Running with --workflow staged --buffer-size medium")
    
    # Find all staged files
    staged_files = find_staged_files()
    
    if not staged_files:
        print("❌ No staged files found!")
        print("💡 Make sure you have *_edited.geojson files in outputs/species_distribution/")
        return
    
    print(f"\n📁 Found {len(staged_files)} staged files:")
    for file in staged_files:
        print(f"   • {file.name}")
    
    # Process each file
    results = []
    total_files = len(staged_files)
    
    for i, file in enumerate(staged_files, 1):
        print(f"\n🔄 PROCESSING {i}/{total_files}")
        
        # Extract species slug from filename
        species_slug = extract_species_slug(file)
        
        if not species_slug:
            print(f"⚠️  Could not extract species slug from: {file.name}")
            results.append({"file": file.name, "status": "skipped", "reason": "filename_parse_error"})
            continue
        
        # Look up species code and name
        if species_slug not in SPECIES_MAPPING:
            print(f"⚠️  Unknown species slug: {species_slug}")
            print(f"💡 Available species: {list(SPECIES_MAPPING.keys())}")
            results.append({"file": file.name, "status": "skipped", "reason": "unknown_species"})
            continue
        
        species_info = SPECIES_MAPPING[species_slug]
        species_code = species_info["code"]
        species_name = species_info["name"]
        
        # Process the species
        success = run_species_processing(species_slug, species_code, species_name)
        
        if success:
            results.append({"file": file.name, "status": "success", "species": species_name})
        else:
            results.append({"file": file.name, "status": "failed", "species": species_name})
    
    # Print summary
    print(f"\n{'='*60}")
    print("📊 PROCESSING SUMMARY")
    print(f"{'='*60}")
    
    successful = [r for r in results if r["status"] == "success"]
    failed = [r for r in results if r["status"] == "failed"]
    skipped = [r for r in results if r["status"] == "skipped"]
    
    print(f"✅ Successful: {len(successful)}/{total_files}")
    for result in successful:
        print(f"   • {result['species']}")
    
    if failed:
        print(f"\n❌ Failed: {len(failed)}")
        for result in failed:
            print(f"   • {result['species']}")
    
    if skipped:
        print(f"\n⚠️  Skipped: {len(skipped)}")
        for result in skipped:
            print(f"   • {result['file']} ({result['reason']})")
    
    print(f"\n🎯 Overall success rate: {len(successful)}/{total_files} ({len(successful)/total_files*100:.1f}%)")

if __name__ == "__main__":
    main() 