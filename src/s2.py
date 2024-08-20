from __future__ import annotations

import argparse
import os
from functools import partial
from multiprocessing import Pool
from pathlib import Path

from rdflib import RDF, RDFS, XSD, Graph, Literal
from rdflib.namespace._GEO import GEO
from s2geometry import (S2Cell, S2CellId, S2LatLng, S2Loop, S2Point, S2Polygon,
                        S2Polyline, S2RegionCoverer)
from shapely.geometry import Polygon

from lib.integrator import Integrator
from lib.kwg_ont import KWGOnt, generate_cell_iri, namespace_prefix


def write_to_rdf(cell_id_int: int, out_path: str, rdf_format: str) -> None:
    """
    Writes a s2 cell to disk

    Args
        cell_id_int: ID of the cell being written
        out_path: The location where the cell will be written
        rdf_format: Format of the RDF. Depends on the formats rdflib supports
    returns
        None
    """
    cell_id = S2CellId(cell_id_int)
    graph = graphify(cell_id)
    if graph:
        for pfx in namespace_prefix:
            graph.bind(pfx, namespace_prefix[pfx])
        file_extensions = {
            "ttl": ".ttl",
            "turtle": ".ttl",
            "xml": ".xml",
            "nq": ".nq",
            "n3": "n3",
            "nt": ".nt",
            "trix": ".trix",
            "trig": ".trig",
            "nquads": ".nq",
            "json-ld": ".jsonld",
        }
        file_name = str(cell_id.id()) + file_extensions[rdf_format]
        destination = os.path.join(out_path, file_name)
        graph.serialize(destination=destination, format=rdf_format)


def graphify(cell_id: S2CellId) -> Graph:
    graph = Graph()

    cell_level = cell_id.level()
    id_int = cell_id.id()

    cell_iri = generate_cell_iri(cell_id)
    p = RDF.type
    o = KWGOnt[f"S2Cell_Level{cell_level}"]
    graph.add((cell_iri, p, o))

    label = f"S2 Cell at level {cell_level} with ID {id_int}"
    p = RDFS.label
    o = Literal(label, datatype=XSD.string)
    graph.add((cell_iri, p, o))

    p = KWGOnt.cellID
    o = Literal(id_int, datatype=XSD.integer)
    graph.add((cell_iri, p, o))

    cell = S2Cell(cell_id)
    area_on_sphere = cell.ApproxArea()
    area_on_earth = area_on_sphere * 6.3781e6 * 6.3781e6

    p = namespace_prefix["geo"]["hasMetricArea"]
    o = Literal(area_on_earth, datatype=XSD.float)
    graph.add((cell_iri, p, o))

    geometry = get_vertex_polygon(cell=cell)
    geometry_iri = KWGOnt.KWGR[f"geometry.polygon.s2.level{cell_level}.{id_int}"]
    p = GEO.hasGeometry
    graph.add((cell_iri, p, geometry_iri))

    p = namespace_prefix["geo"]["hasDefaultGeometry"]
    graph.add((cell_iri, p, geometry_iri))

    p = RDF.type
    o = GEO.Geometry
    graph.add((geometry_iri, p, o))

    o = namespace_prefix["sf"]["Polygon"]
    graph.add((geometry_iri, p, o))

    label = f"Geometry of the polygon formed from the vertices of the S2 Cell at level {cell_level} with ID {id_int}"
    p = RDFS.label
    o = Literal(label, datatype=XSD.string)
    graph.add((geometry_iri, p, o))

    wkt = geometry.wkt
    p = GEO.asWKT
    o = Literal(wkt, datatype=GEO.wktLiteral)
    graph.add((geometry_iri, p, o))

    neighbors = cell_id.GetAllNeighbors(cell_level)
    for neighbor in neighbors:
        p = KWGOnt.sfTouches
        neighbor_iri = generate_cell_iri(neighbor)
        graph.add((cell_iri, p, neighbor_iri))
        graph.add((neighbor_iri, p, cell_iri))

    if cell_level > 0:
        parent = cell_id.parent()
        parent_iri = generate_cell_iri(parent)
        p = KWGOnt.sfWithin
        graph.add((cell_iri, p, parent_iri))

        p = KWGOnt.sfContains
        graph.add((parent_iri, p, cell_iri))

    return graph


def get_vertex_polygon(cell: S2Cell) -> Polygon:
    vertices = map(cell.GetVertex, range(4))
    lat_lngs = [S2LatLng(vertex) for vertex in vertices]
    coords = [
        [lat_lng.lng().degrees(), lat_lng.lat().degrees()] for lat_lng in lat_lngs
    ]
    vertex_polygon = Polygon(coords)
    return vertex_polygon


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--level", type=int, help="Level at which the s2 cells are generated for"
    )
    parser.add_argument(
        "--format",
        help="The format to write the RDF in. Options are xml, n3, turtle, nt, pretty-xml, trix, trig, nquads, "
        "json-ld, hext",
        type=str,
        nargs="?",
        default="ttl",
    )
    parser.add_argument(
        "--ni",
        help="When used, s2 integration is disabled",
        nargs="?",
        const=1,
        type=int,
    )
    parser.add_argument(
        "--compressed",
        help="use the S2 hierarchy to write a compressed collection of relations at various levels",
        type=bool,
        nargs="?",
        default=True,
    )
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
        id_as_integer = (
            current_id.id()
        )  # note that current_id is an instance of S2CellId
        cell_id_integers.append(id_as_integer)
        current_id = current_id.next()
        if current_id == ending_id:
            break
    parents = set()
    print(f"Writing data for cells at level {level}...")
    write = partial(write_to_rdf, out_path=output_path, rdf_format=args.format)
    with Pool() as pool:
        pool.map(write, cell_id_integers)
    if level > 0:
        parents = set(
            S2CellId(cell_id_integer).parent() for cell_id_integer in cell_id_integers
        )
        cell_id_integers = [parent.id() for parent in parents]
        level -= 1

    # Handle integration
    if not args.ni:
        Integrator(args.compressed)
