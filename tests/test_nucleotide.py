"""Tests for view_py.nucleotide."""
import numpy as np

from view_py.nucleotide import Nucleotide, reset_nucleotide_index


class TestNucleotide:
    def setup_method(self):
        reset_nucleotide_index()

    def test_index_increments(self):
        n1 = Nucleotide([0, 0, 0], [1, 0, 0], [0, 0, 1], 0)
        n2 = Nucleotide([1, 0, 0], [1, 0, 0], [0, 0, 1], 1)
        assert n1.index == 0
        assert n2.index == 1

    def test_reset_index(self):
        Nucleotide([0, 0, 0], [1, 0, 0], [0, 0, 1], 0)
        reset_nucleotide_index()
        n = Nucleotide([0, 0, 0], [1, 0, 0], [0, 0, 1], 0)
        assert n.index == 0

    def test_a1_normalized(self):
        n = Nucleotide([0, 0, 0], [3, 0, 0], [0, 0, 1], 0)
        assert abs(np.linalg.norm(n._a1) - 1.0) < 1e-10

    def test_a3_normalized(self):
        n = Nucleotide([0, 0, 0], [1, 0, 0], [0, 0, 5], 0)
        assert abs(np.linalg.norm(n._a3) - 1.0) < 1e-10

    def test_pos_base(self):
        n = Nucleotide([0, 0, 0], [1, 0, 0], [0, 0, 1], 0)
        np.testing.assert_allclose(n.pos_base, [0.4, 0, 0])

    def test_pos_back(self):
        n = Nucleotide([0, 0, 0], [1, 0, 0], [0, 0, 1], 0)
        np.testing.assert_allclose(n.pos_back, [-0.4, 0, 0])

    def test_a2_perpendicular(self):
        n = Nucleotide([0, 0, 0], [1, 0, 0], [0, 0, 1], 0)
        assert abs(np.dot(n.a2, n._a1)) < 1e-10
        assert abs(np.dot(n.a2, n._a3)) < 1e-10

    def test_get_base_string(self):
        for base_int, name in [(0, 'A'), (1, 'G'), (2, 'C'), (3, 'T')]:
            n = Nucleotide([0, 0, 0], [1, 0, 0], [0, 0, 1], base_int)
            assert n.get_base() == name

    def test_base_from_string(self):
        n = Nucleotide([0, 0, 0], [1, 0, 0], [0, 0, 1], 'A')
        assert n._base == 0

    def test_distance(self):
        n1 = Nucleotide([0, 0, 0], [1, 0, 0], [0, 0, 1], 0)
        n2 = Nucleotide([1, 0, 0], [1, 0, 0], [0, 0, 1], 0)
        diff = n1.distance(n2)
        # Backbone positions differ by cm_pos difference (a1 is same)
        assert abs(np.linalg.norm(diff) - 1.0) < 1e-10

    def test_copy(self):
        n = Nucleotide([1, 2, 3], [1, 0, 0], [0, 0, 1], 2)
        n.cluster = 5
        n.color = 123
        c = n.copy()
        np.testing.assert_allclose(c.cm_pos, n.cm_pos)
        assert c._base == n._base
        assert c.cluster == n.cluster
        assert c.color == n.color
        # Ensure it's a deep copy
        c.cm_pos[0] = 999
        assert n.cm_pos[0] != 999
