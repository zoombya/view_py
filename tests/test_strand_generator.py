"""Tests for view_py.strand_generator."""
import math
import numpy as np

from view_py.nucleotide import reset_nucleotide_index
from view_py.strand import reset_strand_index
from view_py.strand_generator import StrandGenerator


class TestStrandGenerator:
    def setup_method(self):
        reset_nucleotide_index()
        reset_strand_index()

    def test_double_helix_returns_two_strands(self):
        gen = StrandGenerator()
        result = gen.generate_or_sq(10, double=True)
        assert isinstance(result, list)
        assert len(result) == 2

    def test_single_strand(self):
        gen = StrandGenerator()
        result = gen.generate_or_sq(10, double=False)
        assert hasattr(result, 'N')
        assert result.N == 10

    def test_correct_nucleotide_count(self):
        gen = StrandGenerator()
        s1, s2 = gen.generate_or_sq(20)
        assert s1.N == 20
        assert s2.N == 20

    def test_pair_assignment(self):
        gen = StrandGenerator()
        s1, s2 = gen.generate_or_sq(10)
        paired = sum(1 for n in s1._nucleotides if n.pair is not None)
        assert paired == 10

    def test_pair_reciprocal(self):
        gen = StrandGenerator()
        s1, s2 = gen.generate_or_sq(5)
        for n in s1._nucleotides:
            assert n.pair is not None
            assert n.pair.pair is n

    def test_custom_angle_array(self):
        gen = StrandGenerator()
        angles = [math.pi / 180 * 34] * 9
        s1, s2 = gen.generate_or_sq(10, angle=angles)
        assert s1.N == 10

    def test_scalar_angle(self):
        gen = StrandGenerator()
        s1, s2 = gen.generate_or_sq(10, angle=0.5)
        assert s1.N == 10

    def test_with_position_and_direction(self):
        gen = StrandGenerator()
        pos = np.array([10.0, 20.0, 30.0])
        direction = np.array([1.0, 0, 0])
        perp = np.array([0, 1.0, 0])
        s1, s2 = gen.generate_or_sq(5, pos=pos, direction=direction, perp=perp)
        assert s1.N == 5
        # First nucleotide should be near the starting position
        dist = np.linalg.norm(s1._nucleotides[0].cm_pos - pos)
        assert dist < 2.0  # within a reasonable distance

    def test_insertion_lengths(self):
        gen = StrandGenerator()
        s1, s2 = gen.generate_or_sq(
            10, lengths=[1], begin=[2], end=[5])
        assert s1.N == 10
