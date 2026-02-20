"""Tests for view_py.system."""
import json
import numpy as np

from view_py.nucleotide import Nucleotide, reset_nucleotide_index
from view_py.strand import Strand, reset_strand_index
from view_py.system import System


def _make_strand(n=3, x_start=0):
    s = Strand()
    for i in range(n):
        s.add_nucleotide(
            Nucleotide([x_start + i * 0.5, 0, 0], [1, 0, 0], [0, 0, 1], i % 4))
    return s


class TestSystem:
    def setup_method(self):
        reset_nucleotide_index()
        reset_strand_index()

    def test_add_strand(self):
        sys = System([100, 100, 100])
        sys.add_strand(_make_strand(5))
        assert sys.N == 5
        assert sys.N_strands == 1

    def test_multiple_strands(self):
        sys = System([100, 100, 100])
        sys.add_strand(_make_strand(5))
        sys.add_strand(_make_strand(3, x_start=10))
        assert sys.N == 8
        assert sys.N_strands == 2

    def test_to_oxview_dict_structure(self):
        sys = System([50, 50, 50])
        sys.add_strand(_make_strand(3))
        result = sys.to_oxview_dict()

        assert 'box' in result
        assert 'systems' in result
        assert len(result['systems']) == 1
        assert 'strands' in result['systems'][0]

    def test_to_oxview_dict_box(self):
        sys = System([123.4, 200, 300])
        sys.add_strand(_make_strand(1))
        result = sys.to_oxview_dict()
        assert result['box'] == [123, 200, 300]

    def test_to_oxview_dict_monomer_fields(self):
        sys = System([50, 50, 50])
        sys.add_strand(_make_strand(2))
        result = sys.to_oxview_dict()
        m = result['systems'][0]['strands'][0]['monomers'][0]

        assert 'id' in m
        assert 'type' in m
        assert 'class' in m
        assert m['class'] == 'DNA'
        assert 'p' in m
        assert 'a1' in m
        assert 'a3' in m
        assert len(m['p']) == 3

    def test_to_oxview_dict_connectivity(self):
        sys = System([50, 50, 50])
        sys.add_strand(_make_strand(3))
        result = sys.to_oxview_dict()
        monomers = result['systems'][0]['strands'][0]['monomers']

        # First monomer: no n3, has n5
        assert 'n3' not in monomers[0]
        assert 'n5' in monomers[0]
        # Middle monomer: has both
        assert 'n3' in monomers[1]
        assert 'n5' in monomers[1]
        # Last monomer: has n3, no n5
        assert 'n3' in monomers[2]
        assert 'n5' not in monomers[2]

    def test_to_oxview_dict_circular_connectivity(self):
        sys = System([50, 50, 50])
        s = _make_strand(3)
        s.make_circular()
        sys.add_strand(s)
        result = sys.to_oxview_dict()
        monomers = result['systems'][0]['strands'][0]['monomers']

        # All monomers should have both n3 and n5
        for m in monomers:
            assert 'n3' in m
            assert 'n5' in m

    def test_to_oxview_dict_pair(self):
        sys = System([50, 50, 50])
        s = Strand()
        n1 = Nucleotide([0, 0, 0], [1, 0, 0], [0, 0, 1], 0)
        n2 = Nucleotide([1, 0, 0], [-1, 0, 0], [0, 0, -1], 3)
        n1.pair = n2
        n2.pair = n1
        s.add_nucleotide(n1)
        s2 = Strand()
        s2.add_nucleotide(n2)
        sys.add_strand(s)
        sys.add_strand(s2)
        result = sys.to_oxview_dict()
        m0 = result['systems'][0]['strands'][0]['monomers'][0]
        m1 = result['systems'][0]['strands'][1]['monomers'][0]
        assert 'pair' in m0
        assert 'pair' in m1
        assert m0['pair'] == m1['id']
        assert m1['pair'] == m0['id']

    def test_to_oxview_string_valid_json(self):
        sys = System([50, 50, 50])
        sys.add_strand(_make_strand(3))
        s = sys.to_oxview_string()
        parsed = json.loads(s)
        assert 'box' in parsed

    def test_calc_clusters_single_strand(self):
        sys = System([50, 50, 50])
        sys.add_strand(_make_strand(3))
        sys.calc_clusters()
        for s in sys._strands:
            for nuc in s._nucleotides:
                assert nuc.cluster == 1
