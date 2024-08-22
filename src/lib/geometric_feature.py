import os
from functools import partial
from pathlib import Path
from typing import Generator

from rdflib import Graph, URIRef
from rdflib.query import ResultRow
from s2geometry import (S2Cell, S2CellId, S2LatLng, S2Loop, S2Point, S2Polygon,
                        S2Polyline, S2RegionCoverer)
from shapely import (LinearRing, LineString, MultiLineString, MultiPolygon,
                     Point, Polygon, buffer)
from shapely.geometry.polygon import signed_area
from shapely.wkt import loads

from .config import config
from .kwg_ont import KWGOnt, generate_cell_iri


class GeometricFeature:
    """
    Represents an abstract geometric feature with an IRI
    """

    def __init__(self, query_solution: ResultRow) -> None:
        self.iri = query_solution["feature_iri"]
        self.geometry = loads(query_solution["wkt"])

    def geometry(self):
        return self.geometry

    def yield_overlapping_ids(
        self, geometry: Polygon | MultiPolygon, tolerance: float = config.tolerance
    ) -> Generator[S2CellId, None, None]:
        """yields the cell IDs of those level 13 cells that overlap
        the given 2D geometry to a certain value of tolerance

        Args:
            geometry (Polygon | MultiPolygon): a 2-dimensional geometry
            tolerance (float, optional): maximal segment width. Defaults to 1e-2.

        Yields:
            Generator[S2CellId, None, None]: a generator through the overlapping IDs
        """
        homogeneous_coverer = S2RegionCoverer()
        homogeneous_coverer.set_min_level(config.min_level)
        homogeneous_coverer.set_max_level(config.max_level)
        for boundary in self.boundaries(geometry):
            segmented_boundary = boundary.segmentize(tolerance)
            buff = buffer(segmented_boundary, tolerance / 100, 2)
            for cell_id in self.covering(
                buff, coverer=homogeneous_coverer, tolerance=tolerance
            ):
                yield cell_id

    def yield_crossing_ids(
        self,
        line_obj: LineString | MultiLineString,
        tolerance: float = config.tolerance,
    ) -> Generator[S2CellId, None, None]:
        """Yields those Cell IDs in a level 13 covering of a
        small (multi)polygon buffer around the given (multi)line string
        to a certain degree of tolerance. The buffer is constructed
        so that its border is at a distance of tolerance/100 from the
        (multi)line string.

        Args:
            line_obj (LineString | MultiLineString): a 1D geometry
            tolerance (float, optional): maximal segment width for the buffer. Defaults to 1e-2.

        Yields:
            Generator[S2CellId, None, None]: a generator of crossing cell IDs
        """
        homogeneous_coverer = S2RegionCoverer()
        homogeneous_coverer.set_min_level(config.min)
        homogeneous_coverer.set_max_level(config.max_level)
        buff = buffer(line_obj, tolerance / 100, 2)
        for cell_id in self.covering(buff, homogeneous_coverer, tolerance=tolerance):
            yield cell_id

    def yield_s2_relations(
        self, coverer: S2RegionCoverer, tolerance: float = config.tolerance
    ) -> Generator[tuple[URIRef, URIRef, URIRef], None, None]:
        if isinstance(self.geometry, (Polygon, MultiPolygon)):
            predicate = KWGOnt.sfContains
            inverse = KWGOnt.sfWithin
            for cell_id in self.filling(self.geometry, coverer, tolerance):
                yield self.iri, predicate, generate_cell_iri(cell_id)
                yield generate_cell_iri(cell_id), inverse, self.iri

            predicate = KWGOnt.sfOverlaps
            for cell_id in self.yield_overlapping_ids(self.geometry, tolerance):
                yield self.iri, predicate, generate_cell_iri(cell_id)
                yield generate_cell_iri(cell_id), predicate, self.iri

        elif isinstance(self.geometry, (LineString, MultiLineString)):
            predicate = KWGOnt.sfCrosses
            for cell_id in self.yield_crossing_ids(self.geometry, tolerance):
                yield self.iri, predicate, generate_cell_iri(cell_id)
                yield generate_cell_iri(cell_id), predicate, self.iri

        elif isinstance(self.geometry, Point):
            s2_point = self.s2_from_coords(self.geometry)
            cell_id = S2CellId(s2_point).parent(config.max_level)
            yield self.iri, KWGOnt.sfWithin, generate_cell_iri(cell_id)
            yield generate_cell_iri(cell_id), KWGOnt.sfContains, self.iri

        else:
            geom_type = self.geometry.geom_type
            msg = f"Geometry of type {geom_type} not supported for s2 relations"
            raise ValueError(msg)

    def s2_graph(
        self, coverer: S2RegionCoverer, tolerance: float = config.tolerance
    ) -> Graph:
        graph = Graph()
        for triple in self.yield_s2_relations(coverer, tolerance):
            graph.add(triple)
        return graph

    def covering(
        self,
        geometry: Polygon | MultiPolygon,
        coverer: S2RegionCoverer,
        tolerance: float = config.tolerance,
    ) -> list[S2CellId]:
        """Returns a list of s2 cell IDs appearing in a homogeneous
        covering of a 2-dimensional geometry by level 13 cells to a
        certain value of tolerance
        Args:
            geometry:
            coverer:
            tolerance (float, optional): maximal segment width. Defaults to 1e-2.
        Returns:
            list[S2CellId]: a list of cell IDs in a level 13 covering
        """
        s2_obj = self.s2_approximation(geometry, tolerance)
        covering = coverer.GetCovering(s2_obj)
        return covering

    def s2_from_coords(
        self,
        geometry: Point | LinearRing | LineString | Polygon | MultiPolygon,
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
            s2_loop.Init(list(map(self.s2_from_coords, list(geometry.coords)[:-1])))
            return s2_loop
        elif isinstance(geometry, LineString):
            polyline = S2Polyline()
            polyline.InitFromS2Points(list(map(self.s2_from_coords, geometry.coords)))
            return polyline
        elif isinstance(geometry, (Polygon, MultiPolygon)):
            loops = map(
                self.s2_from_coords, map(self.orient, self.boundaries(geometry))
            )
            s2_polygon = S2Polygon()
            s2_polygon.InitNested(list(loops))
            return s2_polygon

    def s2_approximation(
        self,
        geometry: S2Point | LinearRing | LineString | Polygon | MultiPolygon,
        tolerance: float = config.tolerance,
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
        return self.s2_from_coords(geometry.segmentize(tolerance))

    def orient(
        self, geometry: LinearRing | Polygon | MultiPolygon, sign: float = 1.0
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
            exterior = self.orient(geometry.exterior, sign)
            oppositely_orient = partial(self.orient, sign=-sign)
            interiors = list(map(oppositely_orient, geometry.interiors))
            return Polygon(exterior, interiors)
        elif isinstance(geometry, MultiPolygon):
            orient_with_sign = partial(self.orient, sign=sign)
            return MultiPolygon(list(map(orient_with_sign, geometry.geoms)))

    @staticmethod
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
        polygons = []
        if isinstance(geometry, Polygon):
            polygons = [geometry]
        elif isinstance(geometry, MultiPolygon):
            polygons = geometry.geoms
        for polygon in polygons:
            yield polygon.exterior
            for interior in polygon.interiors:
                yield interior

    def filling(
        self,
        polygon: Polygon | MultiPolygon,
        coverer: S2RegionCoverer,
        tolerance: float = 1e-2,
    ) -> list[S2CellId]:
        """
        Returns a list of cell IDs that constitute a filling
        of a 2-dimensional geometry that is saturated to the 13th level
        within a certain value of tolerance

        Args:
            polygon (Polygon | MultiPolygon): a 2-dimensional geometry
            coverer (S2RegionCoverer):
            tolerance (float, optional): maximal segment width. Defaults to 1e-2.

        Returns:
            list[S2CellId]: a list of s2 cell IDs in a saturated fill
        """
        filling = []
        s2_obj = self.s2_approximation(polygon, tolerance)
        for exponent in range(4, 9):
            max_cells = 10**exponent
            coverer.set_max_cells(max_cells)
            filling = coverer.GetInteriorCovering(s2_obj)
            num_cells = len(filling)
            if num_cells < 10 ** (exponent - 1):
                break  # iterate until the filling is saturated
        return filling


def yield_geometric_features(path: Path) -> Generator[GeometricFeature, None, None]:
    if os.path.isfile(path):
        graph = Graph()
        with open(path, "r") as read_stream:
            graph.parse(read_stream)
        result = graph.query(
            """
            PREFIX geo: <http://www.opengis.net/ont/geosparql#>
            SELECT ?feature_iri ?wkt 
            WHERE {
                ?feature_iri geo:hasGeometry ?geometry .
                ?geometry geo:asWKT ?wkt .
            }
            """
        )
        for query_solution in result:
            geometric_feature = GeometricFeature(query_solution)
            yield geometric_feature
    elif os.path.isdir(path):
        for file_path in yield_file_paths(path):
            for feature in yield_geometric_features(file_path):
                yield feature


def yield_file_paths(input_dir: Path) -> Generator[Path, None, None]:
    """
    Yields file_paths in input_dir

    Args:
        input_dir (Path): the name of the directory hosting graphical data

    Yields:
        Generator[str, None, None]: a generator of file paths
    """
    for path, _, files in os.walk(input_dir):
        for file in files:
            if not file.startswith("."):
                file_path = Path(os.path.join(path, file))
                yield file_path
