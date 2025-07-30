# Tasks

## Improve Tree Canopy Filtering
90m (factor3) can be quite large in some areas, and may miss more fine details 
30m will be too detailed to reasonably calculate across the entire region

could we do a layered approach?
- add areas to mask with no canopy using low detail filtering
- ignore areas for masking that have high tree canopy at low detail factor 3
- create 1 level of masking for entire region
- use more specific masking when processing a specific species focusing only on that species geographic area bounds
- consider using files like outputs/species_distribution/western_juniper_merged_habitat_regional_habitat_2km_radius_edited.geojson as an inverse mask so that we only focus on creating more detailed masking across these shaded areas
-- in these cases we would use *edited.geojson as the habitat region but we want to our final output to only include areas from that habitat region that have tree canopy. we would analyze the tree canopy data only within the *edited region to limit unnecessary processing