#!/usr/bin/env python3
"""
Bioregion Geometry Utilities

Shared functions for simplifying geometries and removing small holes
to optimize bioregion files for web mapping.
"""

from shapely.geometry import Polygon, MultiPolygon


def simplify_and_remove_holes(geom, tolerance=0.001, hole_threshold=0.001):
    """
    Simplify geometry and remove small holes
    
    Args:
        geom: Shapely geometry (Polygon or MultiPolygon)
        tolerance: Simplification tolerance in degrees (~100m at 0.001)
        hole_threshold: Minimum hole area in degrees² to keep (~0.03 km² at 0.001)
    
    Returns:
        Simplified geometry with small holes removed
    """
    if geom is None or geom.is_empty:
        return geom
    
    # First simplify
    simplified = geom.simplify(tolerance)
    
    # Remove small holes from polygons
    if hasattr(simplified, 'exterior'):  # Single Polygon
        exterior = simplified.exterior
        holes = [hole for hole in simplified.interiors 
                if Polygon(hole).area >= hole_threshold]
        return Polygon(exterior, holes)
    elif hasattr(simplified, 'geoms'):  # MultiPolygon
        cleaned_polys = []
        for poly in simplified.geoms:
            if hasattr(poly, 'exterior'):
                exterior = poly.exterior
                holes = [hole for hole in poly.interiors 
                        if Polygon(hole).area >= hole_threshold]
                cleaned_polys.append(Polygon(exterior, holes))
        return MultiPolygon(cleaned_polys) if cleaned_polys else simplified
    else:
        return simplified


def apply_geometry_optimization(gdf, tolerance=0.0015, hole_threshold=0.0015, coordinate_precision=3):
    """
    Apply complete geometry optimization to a GeoDataFrame
    
    Args:
        gdf: GeoDataFrame to optimize
        tolerance: Simplification tolerance in degrees
        hole_threshold: Minimum hole area to keep in degrees²
        coordinate_precision: Decimal places for coordinate rounding
    
    Returns:
        Optimized GeoDataFrame
    """
    if gdf.empty:
        return gdf
    
    print(f"Optimizing geometries (tolerance={tolerance}, hole_threshold={hole_threshold}, precision={coordinate_precision})")
    
    # Apply simplification and hole removal
    gdf = gdf.copy()
    gdf['geometry'] = gdf.geometry.apply(
        lambda g: simplify_and_remove_holes(g, tolerance, hole_threshold)
    )
    
    # Round coordinates
    from shapely.ops import transform
    
    def round_coordinates(geom):
        """Round coordinates to specified precision"""
        def round_coords(x, y, z=None):
            rounded_x = round(x, coordinate_precision)
            rounded_y = round(y, coordinate_precision)
            return (rounded_x, rounded_y) if z is None else (rounded_x, rounded_y, round(z, coordinate_precision))
        return transform(round_coords, geom)
    
    gdf['geometry'] = gdf.geometry.apply(round_coordinates)
    
    return gdf