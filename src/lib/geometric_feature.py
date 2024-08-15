from rdflib.query import ResultRow
from shapely.wkt import loads


class GeometricFeature:
    """
    Represents an abstract geometric feature with an IRI
    """
    def __init__(self, query_solution: ResultRow) -> None:
        self.iri = query_solution['feature_iri']
        self.geometry = loads(query_solution['wkt'])

    def geometry(self):
        return self.geometry
