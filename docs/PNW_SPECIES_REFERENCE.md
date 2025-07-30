# Pacific Northwest Forest Species Reference

This document provides species codes and names for common tree species found in the Pacific Northwest FIA databases (Washington, Oregon, Idaho).

## Key Species for Habitat Modeling

### Major Conifers
| Code | Common Name | Scientific Name | Notes |
|------|-------------|-----------------|--------|
| 202 | Douglas Fir | *Pseudotsuga menziesii* | Most abundant species in PNW |
| 122 | Ponderosa Pine | *Pinus ponderosa* | Second most abundant |
| 64 | Western Juniper | *Juniperus occidentalis* | ~14k trees in OR, good for arid habitat studies |
| 73 | Western Larch | *Larix occidentalis* | Deciduous conifer |
| 108 | Lodgepole Pine | *Pinus contorta* | High elevation specialist |
| 15 | White Fir | *Abies concolor* | Mixed conifer forests |
| 17 | Grand Fir | *Abies grandis* | Moist forest sites |
| 263 | Western Hemlock | *Tsuga heterophylla* | Coastal/moist forests |
| 242 | Pacific Silver Fir | *Abies amabilis* | High elevation |
| 264 | Mountain Hemlock | *Tsuga mertensiana* | Subalpine |

### Major Hardwoods  
| Code | Common Name | Scientific Name | Notes |
|------|-------------|-----------------|--------|
| 351 | Red Alder | *Alnus rubra* | Riparian/early succession |
| 312 | Bigleaf Maple | *Acer macrophyllum* | Pacific coast deciduous |
| 631 | Pacific Madrone | *Arbutus menziesii* | Mediterranean climate indicator |
| 361 | Oregon White Oak | *Quercus garryana* | Oak woodland/savanna |
| 11 | Balsam Poplar | *Populus balsamifera* | Riparian specialist |
| 81 | Quaking Aspen | *Populus tremuloides* | Montane/boreal |

### Other Notable Species
| Code | Common Name | Scientific Name | Notes |
|------|-------------|-----------------|--------|
| 69 | Oneseed Juniper | *Juniperus monosperma* | ⚠️ NOT Western Juniper |
| 242 | Noble Fir | *Abies procera* | High elevation, limited range |
| 93 | Engelmann Spruce | *Engelmannia* | Subalpine specialist |
| 374 | Water Birch | *Betula occidentalis* | Riparian/montane |

## Database Query Examples

### Find all species with 1000+ trees in Oregon:
```sql
SELECT DISTINCT t.SPCD, r.COMMON_NAME, r.SCIENTIFIC_NAME, COUNT(*) as tree_count 
FROM TREE t 
JOIN REF_SPECIES r ON t.SPCD = r.SPCD 
WHERE t.STATUSCD = 1 
GROUP BY t.SPCD, r.COMMON_NAME, r.SCIENTIFIC_NAME 
HAVING tree_count >= 1000 
ORDER BY tree_count DESC;
```

### Search for species by name:
```sql
SELECT SPCD, COMMON_NAME, SCIENTIFIC_NAME 
FROM REF_SPECIES 
WHERE COMMON_NAME LIKE '%juniper%' 
   OR SCIENTIFIC_NAME LIKE '%Juniperus%';
```

## Common Mistakes to Avoid

1. **Western Juniper**: Use code **64**, not 69 (which is Oneseed Juniper)
2. **Douglas Fir**: Use code **202** (most abundant PNW species)
3. **Ponderosa Pine**: Use code **122** (second most abundant)

## Species Distribution Script Usage

### For Western Juniper:
```bash
export CUSTOM_SPECIES_CODE=64
export CUSTOM_SPECIES_NAME="Western Juniper"
python scripts/visualization/create_species_distribution_maps_cached.py --workflow staged
```

### For Custom Species:
```bash
export CUSTOM_SPECIES_CODE=<code>
export CUSTOM_SPECIES_NAME="<name>"
python scripts/visualization/create_species_distribution_maps_cached.py
```

## Data Sources
- FIA Database: `data/raw/trees_SQLite_FIADB_{WA,OR,ID}.db`
- Species Reference Table: `REF_SPECIES`
- Tree Data Table: `TREE` (use `STATUSCD = 1` for live trees)

Last Updated: $(date) 