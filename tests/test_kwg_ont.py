from rdflib import URIRef
from s2geometry import S2CellId

from ..src.lib.kwg_ont import generate_cell_iri


def test_generate_cell_iri():
    """
    Tests to ensure the S2 cell IRI is the expected form

    :return: None
    """
    cell_id_int = 288230376151711744
    cell_id = S2CellId(cell_id_int)
    assert generate_cell_iri(cell_id) == URIRef("http://stko-kwg.geog.ucsb.edu/lod/resource/s2.level1"
                                                ".288230376151711744")
