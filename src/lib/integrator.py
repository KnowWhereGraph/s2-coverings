from __future__ import annotations

import os
from functools import partial
from multiprocessing import Pool
from pathlib import Path

from .config import config
from .kwg_ont import namespace_prefix
from .geometric_feature import GeometricFeature, yield_geometric_features
from s2geometry import (S2CellId, S2RegionCoverer)
from shapely.geometry import MultiPolygon, Polygon


class Integrator:
    """
    Abstraction over the process for integrating s2 cells together with spatial relations.
    """

    def __init__(self, compressed: bool):
        if compressed:
            print("Compression is on. Relations will be compressed using the S2 hierarchy...")
        data_path = Path("./output/")
        output_folder = f"./output/{data_path.stem}"

        if compressed:
            output_folder += "_compressed"

        os.makedirs(output_folder, exist_ok=True)

        write = partial(
            self.write_all_relations,
            output_folder=output_folder,
            is_compressed=compressed,
        )

        with Pool() as pool:
            pool.map(write, enumerate(yield_geometric_features(data_path)))

        print(f"Done! \nRelations written in path '{output_folder}'.")

    @staticmethod
    def homogeneous_covering(
            geometry: Polygon | MultiPolygon,
            level: int,
            tolerance: float = config.tolerance,
    ) -> list[S2CellId]:
        homogeneous_coverer = S2RegionCoverer()
        homogeneous_coverer.set_min_level(level)
        homogeneous_coverer.set_max_level(level)
        return GeometricFeature().covering(
            geometry=geometry, coverer=homogeneous_coverer, tolerance=tolerance
        )

    @staticmethod
    def write_all_relations(
            indexed_feature: tuple[int, GeometricFeature],
            output_folder: str,
            is_compressed: bool,
            tolerance: float = config.tolerance,
    ) -> None:
        idx, feature = indexed_feature
        coverer = S2RegionCoverer()
        coverer.set_max_level(config.max_level)
        if not is_compressed:
            coverer.set_min_level(config.min_level)
        else:
            coverer.set_min_level(0)
        graph = feature.s2_graph(coverer, tolerance=tolerance)
        if graph:
            for pfx in namespace_prefix:
                graph.bind(pfx, namespace_prefix[pfx])
            file_name = f"{idx}.ttl"
            destination = os.path.join(output_folder, file_name)
            graph.serialize(destination=destination, format="ttl")
