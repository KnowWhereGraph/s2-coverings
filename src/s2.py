from __future__ import annotations
from pathlib import Path
from typing import Generator
import os
import argparse
from functools import partial
from multiprocessing import Pool

from rdflib import Graph, Namespace, Literal, URIRef, RDF, RDFS, XSD
from rdflib.namespace._GEO import GEO

from shapely.geometry import Point, LinearRing, Polygon, MultiPolygon, LineString
from shapely.geometry.polygon import signed_area
from s2geometry import S2CellId, S2RegionCoverer, S2Point, S2LatLng, S2Loop, S2Polyline, S2Polygon, S2Cell

from lib.geometric_feature import GeometricFeature
from lib.kwg_ont import KWGOnt

TOLERANCE = 1e-2

_PREFIX = {
    "kwgr": KWGOnt.KWGR,
    "kwg-ont": KWGOnt._NS,
    "geo": Namespace("http://www.opengis.net/ont/geosparql#"),
    "sf": Namespace("http://www.opengis.net/ont/sf#"),
    "rdf": RDF,
    "rdfs": RDFS,
    "xsd": XSD,
}


def create_cell_iri(cell_id: S2CellId) -> URIRef:
    """
    Creates an IRI for an individual cell, with a KnowWhereGraph domain

    Args:
        cell_id: The ID of the s2 cell
    Returns:
         A URI of the s2 cell
    """
    level = cell_id.level()
    id_str = cell_id.id()
    return KWGOnt.KWGR[f"{'s2.level'}{level}.{id_str}"]


def s2_from_coords(
        geometry: tuple | Point | LinearRing | LineString | Polygon | MultiPolygon
) -> S2Point | S2Loop | S2Polyline | S2Polygon:
    """
    Returns a corresponding S2 object whose vertices
    are obtained from the coordinates of the given geometry.
    Note that the returned object lives on the sphere and any
    linear edges are replaced by geodesic curves. As such, the resulting
    S2 object is an APPROXIMATION of the given geometry living on
    the sphere, and the approximation is worse when vertices are farther apart
    (because geodesics will curve more in that case).
    Args:
        geometry (tuple | Point | LinearRing | LineString | Polygon | MultiPolygon): a geometry
    Returns:
        S2Point | S2Loop | S2Polyline | S2Polygon: an S2 geometric object
    """
    if isinstance(geometry, tuple):
        return S2LatLng.FromDegrees(*geometry[::-1]).ToPoint()
    elif isinstance(geometry, Point):
        return S2LatLng.FromDegrees(*geometry.coords[0][::-1]).ToPoint()
    elif isinstance(geometry, LinearRing):
        s2_loop = S2Loop()
        s2_loop.Init(list(map(s2_from_coords, list(geometry.coords)[:-1])))
        return s2_loop
    elif isinstance(geometry, LineString):
        polyline = S2Polyline()
        polyline.InitFromS2Points(list(map(s2_from_coords, geometry.coords)))
        return polyline
    elif isinstance(geometry, (Polygon, MultiPolygon)):
        loops = map(s2_from_coords, map(orient, boundaries(geometry)))
        s2_polygon = S2Polygon()
        s2_polygon.InitNested(list(loops))
        return s2_polygon


def s2_approximation(
        geometry: S2Point | LinearRing | LineString | Polygon | MultiPolygon,
        tolerance: float = TOLERANCE
) -> S2Point | S2Loop | S2Polyline | S2Polygon:
    """
    Returns a corresponding S2 object that approximates the given
    Shapely object to a certain value of tolerance. Precisely, a given
    Shapely geometry is first segmentized, meaning extra vertices are added
    so that the width between any two adjacent vertices never exceeds the
    tolerance, and then these new coordinates are used as the vertices
    of an S2 object (see s2_from_coords).
    Args:
        geometry (LinearRing | LineString | Polygon | MultiPolygon): a Shapely geometry
        tolerance (float, optional): a maximal segment width Defaults to 1e-2.
    Returns:
        S2Point | S2Loop | S2Polyline | S2Polygon: an S2 geometry object
    """
    return s2_from_coords(geometry.segmentize(tolerance))


def orient(
        geometry: LinearRing | Polygon | MultiPolygon,
        sign: float = 1.0
) -> LinearRing | Polygon | MultiPolygon:
    """
    Returns a copy of the geometry with the specified orientation
    Args:
        geometry (LinearRing | Polygon | MultiPolygon): a ring or (multi)polygon
        sign (float, optional): an orientation. Defaults to 1.0.
    Returns:
        LinearRing | Polygon | MultiPolygon: a copy of the geometry with specified orientation
    """
    sign = float(sign)
    if isinstance(geometry, LinearRing):
        return geometry if signed_area(geometry) / sign >= 0 else geometry.reverse()
    elif isinstance(geometry, Polygon):
        exterior = orient(geometry.exterior, sign)
        oppositely_orient = partial(orient, sign=-sign)
        interiors = list(map(oppositely_orient, geometry.interiors))
        return Polygon(exterior, interiors)
    elif isinstance(geometry, MultiPolygon):
        orient_with_sign = partial(orient, sign=sign)
        return MultiPolygon(list(map(orient_with_sign, geometry.geoms)))


def boundaries(
        geometry: Polygon | MultiPolygon,
) -> Generator[LinearRing, None, None]:
    """
    Yields the boundary rings of a geometry with boundaries
    Args:
        geometry (Polygon | MultiPolygon): a shapely geometry with boundaries
    Yields:
        Generator[LinearRing, None, None]: a generator through the boundary rings
    """
    if isinstance(geometry, Polygon):
        polygons = [geometry]
    elif isinstance(geometry, MultiPolygon):
        polygons = geometry.geoms
    for polygon in polygons:
        yield polygon.exterior
        for interior in polygon.interiors:
            yield interior


def covering(
        geometry: Polygon | MultiPolygon,
        coverer: S2RegionCoverer,
        tolerance: float = TOLERANCE
) -> list[S2CellId]:
    """Returns a list of s2 cell IDs appearing in a homogeneous
    covering of a 2-dimensional geometry by level 13 cells to a
    certain value of tolerance
    Args:
        polygon (Polygon): a 2 dimensional geometry
        tolerance (float, optional): maximal segment width. Defaults to 1e-2.
    Returns:
        list[S2CellId]: a list of cell IDs in a level 13 covering
    """
    s2_obj = s2_approximation(geometry, tolerance)
    covering = coverer.GetCovering(s2_obj)
    return covering


def yield_file_paths(input_dir: str) -> Generator[str, None, None]:
    """yields those file_paths in input_dir except for .DS_Store,
    which is a file created by walking

    Args:
        input_dir (str): the name of the directory hosting graphical data

    Yields:
        Generator[str, None, None]: a generator of file paths
    """
    for (path, _, files) in os.walk(input_dir):
        for file in files:
            if not file == ".DS_Store":
                file_path = os.path.join(path, file)
                yield file_path


def yield_geometric_features(path: str) -> Generator[GeometricFeature, None, None]:
    if os.path.isfile(path):
        graph = Graph()
        with open(path, 'r') as read_stream:
            graph.parse(read_stream)
        result = graph.query("""
        PREFIX geo: <http://www.opengis.net/ont/geosparql#>
        SELECT ?feature_iri ?wkt 
        WHERE {
            ?feature_iri geo:hasGeometry ?geometry .
            ?geometry geo:asWKT ?wkt .
        }
        """)
        for query_solution in result:
            geometric_feature = GeometricFeature(query_solution)
            yield geometric_feature
    elif os.path.isdir(path):
        for file_path in yield_file_paths(path):
            for feature in yield_geometric_features(file_path):
                yield feature


def get_parents(cell_ids: list[S2CellId]) -> list[S2CellId]:
    if cell_ids[0].level() == 0:
        return
    parents = set()
    for cell_id in cell_ids:
        parents.add(cell_id.parent())
    return list(parents)


def write_to_rdf(cell_id_int: int, output_path: str) -> None:
    cell_id = S2CellId(cell_id_int)
    graph = graphify(cell_id)
    if graph:
        for prefix in _PREFIX:
            graph.bind(prefix, _PREFIX[prefix])
        file_name = str(cell_id.id()) + ".ttl"
        destination = os.path.join(output_path, file_name)
        graph.serialize(destination=destination, format="ttl")


def graphify(cell_id: S2CellId) -> Graph:
    graph = Graph()

    level = cell_id.level()
    id_int = cell_id.id()

    cell_iri = create_cell_iri(cell_id)
    p = RDF.type
    o = KWGOnt[f"S2Cell_Level{level}"]
    graph.add((cell_iri, p, o))

    label = f"S2 Cell at level {level} with ID {id_int}"
    p = RDFS.label
    o = Literal(label, datatype=XSD.string)
    graph.add((cell_iri, p, o))

    p = KWGOnt.cellID
    o = Literal(id_int, datatype=XSD.integer)
    graph.add((cell_iri, p, o))

    # p = KWG_ONT.cellLevel
    # o = Literal(level, datatype=XSD.integer)
    # graph.add((cell_iri, p, o))

    cell = S2Cell(cell_id)
    area_on_sphere = cell.ApproxArea()
    area_on_earth = area_on_sphere * (6.3781e6) * (6.3781e6)

    p = _PREFIX["geo"]["hasMetricArea"]
    o = Literal(area_on_earth, datatype=XSD.float)
    graph.add((cell_iri, p, o))

    geometry = get_vertex_polygon(cell=cell)
    geometry_iri = KWGOnt.KWGR[f"geometry.polygon.s2.level{level}.{id_int}"]
    p = GEO.hasGeometry
    graph.add((cell_iri, p, geometry_iri))

    p = _PREFIX["geo"]["hasDefaultGeometry"]
    graph.add((cell_iri, p, geometry_iri))

    p = RDF.type
    o = GEO.Geometry
    graph.add((geometry_iri, p, o))

    o = _PREFIX["sf"]["Polygon"]
    graph.add((geometry_iri, p, o))

    label = f"Geometry of the polygon formed from the vertices of the S2 Cell at level {level} with ID {id_int}"
    p = RDFS.label
    o = Literal(label, datatype=XSD.string)
    graph.add((geometry_iri, p, o))

    wkt = geometry.wkt
    p = GEO.asWKT
    o = Literal(wkt, datatype=GEO.wktLiteral)
    graph.add((geometry_iri, p, o))

    neighbors = cell_id.GetAllNeighbors(level)
    for neighbor in neighbors:
        p = KWGOnt.sfTouches
        neighbor_iri = create_cell_iri(neighbor)
        graph.add((cell_iri, p, neighbor_iri))
        graph.add((neighbor_iri, p, cell_iri))

    if level > 0:
        parent = cell_id.parent()
        parent_iri = create_cell_iri(parent)
        p = KWGOnt.sfWithin
        graph.add((cell_iri, p, parent_iri))

        p = KWGOnt.sfContains
        graph.add((parent_iri, p, cell_iri))

    return graph


def get_vertex_polygon(cell: S2Cell) -> Polygon:
    vertices = map(cell.GetVertex, range(4))
    lat_lngs = [S2LatLng(vertex) for vertex in vertices]
    coords = [[lat_lng.lng().degrees(), lat_lng.lat().degrees()] for lat_lng in lat_lngs]
    vertex_polygon = Polygon(coords)
    return vertex_polygon


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("level", type=int, help="Level at which the s2 cells are generated for")
    args = parser.parse_args()

    level = args.level
    output_folder = f"output"
    level_path = Path(f"level_{level}")
    output_path = os.path.join(output_folder, level_path)
    os.makedirs(output_path, exist_ok=True)

    beginning_id = S2CellId.Begin(level=level)
    ending_id = S2CellId.End(level=level)

    current_id = beginning_id
    cell_id_integers = []
    while True:
        id_as_integer = current_id.id()  # note that current_id is an instance of S2CellId
        print(id_as_integer)
        cell_id_integers.append(id_as_integer)
        current_id = current_id.next()
        if current_id == ending_id:
            break
    parents = set()
    print(f"Writing data for cells at level {level}...")
    write = partial(write_to_rdf, output_path=output_path)
    with Pool() as pool:
        pool.map(write, cell_id_integers)
    if level > 0:
        parents = set(S2CellId(cell_id_integer).parent() for cell_id_integer in cell_id_integers)
        cell_id_integers = [parent.id() for parent in parents]
        level -= 1