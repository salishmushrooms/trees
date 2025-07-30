"""
GBIF Mushroom Data Exploration
Explore mushroom observations for target species in Pacific Northwest
"""

import requests
import pandas as pd
import json
from datetime import datetime
import time

# Target mushroom species for Pacific Northwest habitat modeling
TARGET_SPECIES = {
    'chanterelles': [
        'Cantharellus formosus', 
        'Cantharellus cascadensis',
        'Cantharellus formosus',
        'Cantharellus subalpinus'
    ],
    'boletes': [
        'Boletus edulis',
        'Boletus rex-veris',
        'Boletus barrowsii',
        'Boletus regineus'
    ],
    'morels': [ 
        'Morchella esculenta',
        'Morchella elata',
        'Morchella importuna', #easiest to associate with weather conditions and many observations in city
        'Morchella snyderi',
        'Morchella tridentina',
        'Morchella americana', #fairly reliable identifications and habitat distribution
        'morchella norvegienses', #possibly more difficult to identify, broad habitat distribution in PNW
        'Morchella populiphila',
        'Morchella tomentosa', #reliable ID and found in burn habitats
        'Morchella exuberans', #burn morels one of them
    ],
    'matsutake': [
        'Tricholoma magnivelare', #consider adding section matsutake if possible
        'Tricholoma caligatum' #false matsutake or false positive
    ]
}

# Pacific Northwest bounding box (approximate)
# Washington and Oregon
PNW_BOUNDS = {
    'min_lat': 42.0,  # Southern Oregon
    'max_lat': 49.0,  # Northern Washington
    'min_lon': -125.0, # Pacific Coast
    'max_lon': -116.0  # Eastern Washington/Oregon
}

def query_gbif_species(species_name, limit=300):
    """
    Query GBIF for occurrences of a specific species in Pacific Northwest
    """
    base_url = "https://api.gbif.org/v1/occurrence/search"
    
    params = {
        'scientificName': species_name,
        'decimalLatitude': f"{PNW_BOUNDS['min_lat']},{PNW_BOUNDS['max_lat']}",
        'decimalLongitude': f"{PNW_BOUNDS['min_lon']},{PNW_BOUNDS['max_lon']}",
        'hasCoordinate': True,
        'hasGeospatialIssue': False,
        'limit': limit,
        'country': 'US',  # United States only
        'year': '1990,2024'  # Recent decades only
    }
    
    try:
        response = requests.get(base_url, params=params)
        response.raise_for_status()
        return response.json()
    except requests.RequestException as e:
        print(f"Error querying GBIF for {species_name}: {e}")
        return None

def extract_occurrence_data(gbif_response, species_name):
    """
    Extract relevant fields from GBIF response
    """
    if not gbif_response or 'results' not in gbif_response:
        return []
    
    records = []
    for record in gbif_response['results']:
        # Extract key fields
        occurrence = {
            'species_name': species_name,
            'gbif_id': record.get('key'),
            'latitude': record.get('decimalLatitude'),
            'longitude': record.get('decimalLongitude'),
            'coordinate_uncertainty': record.get('coordinateUncertaintyInMeters'),
            'elevation': record.get('elevation'),
            'date': record.get('eventDate'),
            'year': record.get('year'),
            'month': record.get('month'),
            'day': record.get('day'),
            'dataset': record.get('datasetName'),
            'institution': record.get('institutionCode'),
            'collector': record.get('recordedBy'),
            'country': record.get('country'),
            'state_province': record.get('stateProvince'),
            'locality': record.get('locality'),
            'habitat': record.get('habitat'),
            'basis_of_record': record.get('basisOfRecord'),
            'occurrence_status': record.get('occurrenceStatus'),
            'taxon_key': record.get('taxonKey'),
            'scientific_name': record.get('scientificName'),
            'identification_qualifier': record.get('identificationQualifier')
        }
        records.append(occurrence)
    
    return records

def main():
    """
    Main function to explore GBIF data for all target species
    """
    print("Exploring GBIF mushroom occurrence data for Pacific Northwest...")
    print(f"Search area: {PNW_BOUNDS['min_lat']:.1f}°N to {PNW_BOUNDS['max_lat']:.1f}°N, "
          f"{PNW_BOUNDS['min_lon']:.1f}°W to {PNW_BOUNDS['max_lon']:.1f}°W")
    print("=" * 70)
    
    all_records = []
    species_summary = []
    
    for group_name, species_list in TARGET_SPECIES.items():
        print(f"\n{group_name.upper()}:")
        print("-" * 40)
        
        group_records = []
        for species in species_list:
            print(f"Querying {species}...")
            
            # Query GBIF
            gbif_data = query_gbif_species(species)
            
            if gbif_data:
                records = extract_occurrence_data(gbif_data, species)
                group_records.extend(records)
                all_records.extend(records)
                
                # Summary stats
                total_found = gbif_data.get('count', 0)
                records_returned = len(records)
                
                print(f"  Total occurrences: {total_found:,}")
                print(f"  Records with coordinates: {records_returned:,}")
                
                if records:
                    # Quick stats
                    df_temp = pd.DataFrame(records)
                    year_range = f"{df_temp['year'].min()}-{df_temp['year'].max()}"
                    elevation_range = ""
                    if df_temp['elevation'].notna().any():
                        elev_min = df_temp['elevation'].min()
                        elev_max = df_temp['elevation'].max()
                        elevation_range = f", elevation: {elev_min}-{elev_max}m"
                    
                    print(f"  Year range: {year_range}{elevation_range}")
                
                species_summary.append({
                    'group': group_name,
                    'species': species,
                    'total_occurrences': total_found,
                    'coordinate_records': records_returned
                })
            
            # Be nice to GBIF API
            time.sleep(0.5)
    
    # Create comprehensive DataFrame
    if all_records:
        df_all = pd.DataFrame(all_records)
        
        # Clean and process data
        df_all['coordinate_uncertainty'] = pd.to_numeric(df_all['coordinate_uncertainty'], errors='coerce')
        df_all['elevation'] = pd.to_numeric(df_all['elevation'], errors='coerce')
        df_all['latitude'] = pd.to_numeric(df_all['latitude'], errors='coerce')
        df_all['longitude'] = pd.to_numeric(df_all['longitude'], errors='coerce')
        
        # Filter for reasonable coordinate uncertainty (≤1km)
        df_filtered = df_all[
            (df_all['coordinate_uncertainty'].isna()) | 
            (df_all['coordinate_uncertainty'] <= 1000)
        ].copy()
        
        # Save to files
        df_all.to_csv('gbif_mushroom_data_raw.csv', index=False)
        df_filtered.to_csv('gbif_mushroom_data_filtered.csv', index=False)
        
        print(f"\n{'='*70}")
        print("OVERALL SUMMARY")
        print(f"{'='*70}")
        print(f"Total records found: {len(df_all):,}")
        print(f"Records with coordinate uncertainty ≤1km: {len(df_filtered):,}")
        
        if len(df_filtered) > 0:
            print(f"Date range: {df_filtered['year'].min()}-{df_filtered['year'].max()}")
            print(f"Elevation range: {df_filtered['elevation'].min():.0f}m to {df_filtered['elevation'].max():.0f}m")
            print(f"States: {df_filtered['state_province'].value_counts().index.tolist()}")
            
            # Monthly distribution
            print(f"\nFruiting seasonality (month):")
            monthly = df_filtered['month'].value_counts().sort_index()
            for month, count in monthly.items():
                if pd.notna(month):
                    month_name = datetime(2000, int(month), 1).strftime('%B')
                    print(f"  {month_name}: {count} observations")
            
            # Species summary
            print(f"\nSpecies summary:")
            species_counts = df_filtered['species_name'].value_counts()
            for species, count in species_counts.items():
                print(f"  {species}: {count} records")
        
        print(f"\nFiles saved:")
        print(f"  gbif_mushroom_data_raw.csv ({len(df_all):,} records)")
        print(f"  gbif_mushroom_data_filtered.csv ({len(df_filtered):,} records)")
        
    else:
        print("No records found for any target species.")
    
    # Save species summary
    pd.DataFrame(species_summary).to_csv('gbif_species_summary.csv', index=False)
    print(f"  gbif_species_summary.csv")

if __name__ == "__main__":
    main() 