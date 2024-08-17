from rdflib import RDF, RDFS, XSD, URIRef
from rdflib.namespace import DefinedNamespace, Namespace
from s2geometry import S2Cell, S2CellId


class KWGOnt(DefinedNamespace):
    kwg_endpoint = "http://stko-kwg.geog.ucsb.edu/"

    KWGR = Namespace(f"{kwg_endpoint}lod/resource/")

    S2Cell_Level0: URIRef
    S2Cell_Level1: URIRef
    S2Cell_Level2: URIRef
    S2Cell_Level3: URIRef
    S2Cell_Level4: URIRef
    S2Cell_Level5: URIRef
    S2Cell_Level6: URIRef
    S2Cell_Level7: URIRef
    S2Cell_Level8: URIRef
    S2Cell_Level9: URIRef
    S2Cell_Level10: URIRef
    S2Cell_Level11: URIRef
    S2Cell_Level12: URIRef
    S2Cell_Level13: URIRef

    cellID: URIRef

    sfEquals: URIRef
    sfContains: URIRef
    sfWithin: URIRef
    sfTouches: URIRef
    sfOverlaps: URIRef
    sfCrosses: URIRef
    vertexPolygon: URIRef

    _NS = Namespace(f"{kwg_endpoint}lod/ontology/")


def generate_cell_iri(cell_id: S2CellId) -> URIRef:
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


namespace_prefix = {
    "kwgr": KWGOnt.KWGR,
    "kwg-ont": KWGOnt._NS,
    "geo": Namespace("http://www.opengis.net/ont/geosparql#"),
    "sf": Namespace("http://www.opengis.net/ont/sf#"),
    "rdf": RDF,
    "rdfs": RDFS,
    "xsd": XSD,
}
