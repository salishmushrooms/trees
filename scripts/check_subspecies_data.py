#!/usr/bin/env python3
"""
Check FIA database for subspecies information
"""

import sqlite3
from pathlib import Path

def check_database(db_path):
    """Check database for subspecies information"""
    print(f"\n{'='*60}")
    print(f"Checking: {db_path}")
    print('='*60)
    
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        
        # Check REF_SPECIES table structure
        print("\n1. REF_SPECIES Table Structure:")
        print("-" * 40)
        cursor.execute("PRAGMA table_info(REF_SPECIES)")
        columns = cursor.fetchall()
        for col in columns:
            print(f"  {col[1]:20} {col[2]}")
        
        # Check Pinus contorta entries
        print("\n2. Pinus contorta (SPCD 108) Details:")
        print("-" * 40)
        cursor.execute("""
            SELECT SPCD, COMMON_NAME, GENUS, SPECIES, VARIETY, SUBSPECIES, SCIENTIFIC_NAME
            FROM REF_SPECIES 
            WHERE SPCD = 108
        """)
        
        results = cursor.fetchall()
        for row in results:
            print(f"  SPCD: {row[0]}")
            print(f"  Common Name: {row[1]}")
            print(f"  Genus: {row[2]}")
            print(f"  Species: {row[3]}")
            print(f"  Variety: {row[4] if row[4] else 'NULL'}")
            print(f"  Subspecies: {row[5] if row[5] else 'NULL'}")
            print(f"  Scientific Name: {row[6]}")
        
        # Check if any Pinus species have subspecies data
        print("\n3. All Pinus Species with Subspecies Data:")
        print("-" * 40)
        cursor.execute("""
            SELECT SPCD, COMMON_NAME, GENUS, SPECIES, VARIETY, SUBSPECIES
            FROM REF_SPECIES 
            WHERE GENUS = 'Pinus' AND (VARIETY IS NOT NULL OR SUBSPECIES IS NOT NULL)
            ORDER BY SPCD
        """)
        
        results = cursor.fetchall()
        if results:
            for row in results:
                print(f"  SPCD {row[0]}: {row[1]} - Variety: {row[4]}, Subspecies: {row[5]}")
        else:
            print("  No Pinus species have subspecies data")
        
        # Count Pinus contorta trees in coastal areas
        print("\n4. Pinus contorta Trees in Coastal Areas:")
        print("-" * 40)
        
        # Check if PLOT table exists and get coastal plots
        cursor.execute("""
            SELECT COUNT(*) as total_trees,
                   COUNT(DISTINCT p.CN) as unique_plots,
                   MIN(p.ELEV) as min_elevation,
                   MAX(p.ELEV) as max_elevation,
                   AVG(p.ELEV) as avg_elevation
            FROM TREE t
            JOIN PLOT p ON t.PLT_CN = p.CN
            WHERE t.SPCD = 108 
            AND p.LON < -123.0 
            AND p.ELEV < 1000
            AND p.PLOT_STATUS_CD = 1
        """)
        
        result = cursor.fetchone()
        if result[0] > 0:
            print(f"  Total coastal Pinus contorta trees: {result[0]}")
            print(f"  In {result[1]} unique plots")
            print(f"  Elevation range: {result[2]:.0f} - {result[3]:.0f} ft")
            print(f"  Average elevation: {result[4]:.0f} ft")
        else:
            print("  No coastal Pinus contorta trees found")
        
        conn.close()
        
    except Exception as e:
        print(f"Error: {e}")

def main():
    """Check all available FIA databases"""
    print("FIA Database Subspecies Analysis")
    print("================================")
    
    # List of databases to check
    databases = [
        "data/raw/trees_SQLite_FIADB_WA.db",
        "data/raw/trees_SQLite_FIADB_OR.db",
        "data/raw/trees_SQLite_FIADB_ID.db"
    ]
    
    for db_path in databases:
        if Path(db_path).exists():
            check_database(db_path)
        else:
            print(f"\nSkipping {db_path} - file not found")
    
    print(f"\n{'='*60}")
    print("SUMMARY:")
    print("="*60)
    print("- FIA database has VARIETY and SUBSPECIES fields")
    print("- These fields are NOT populated for Pinus contorta (SPCD 108)")
    print("- Shore pine vs. lodgepole pine must be distinguished by location")
    print("- Recommendation: Use coastal proximity + elevation as proxies")

if __name__ == "__main__":
    main()