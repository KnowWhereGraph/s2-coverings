from rdflib import URIRef
from rdflib.namespace import DefinedNamespace, Namespace


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