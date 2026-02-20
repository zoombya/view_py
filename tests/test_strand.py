"""Tests for view_py.strand."""
import numpy as np

from view_py.nucleotide import Nucleotide, reset_nucleotide_index
from view_py.strand import Strand, reset_strand_index


def _make_nuc(x=0):
    return Nucleotide([x, 0, 0], [1, 0, 0], [0, 0, 1], 0)


class TestStrand:
    def setup_method(self):
        reset_nucleotide_index()
        reset_strand_index()

    def test_add_nucleotide(self):
        s = Strand()
        s.add_nucleotide(_make_nuc(0))
        s.add_nucleotide(_make_nuc(1))
        assert s.N == 2

    def test_sequence_tracks_bases(self):
        s = Strand()
        n1 = Nucleotide([0, 0, 0], [1, 0, 0], [0, 0, 1], 0)  # A
        n2 = Nucleotide([1, 0, 0], [1, 0, 0], [0, 0, 1], 3)  # T
        s.add_nucleotide(n1)
        s.add_nucleotide(n2)
        assert s.sequence == [0, 3]

    def test_set_sequence(self):
        s = Strand()
        s.add_nucleotide(_make_nuc(0))
        s.add_nucleotide(_make_nuc(1))
        s.sequence = "GC"
        assert s._nucleotides[0]._base == 1  # G
        assert s._nucleotides[1]._base == 2  # C

    def test_append(self):
        s1 = Strand()
        s1.add_nucleotide(_make_nuc(0))
        s2 = Strand()
        s2.add_nucleotide(_make_nuc(1))
        s3 = s1.append(s2)
        assert s3.N == 2

    def test_get_slice(self):
        s = Strand()
        for i in range(5):
            s.add_nucleotide(_make_nuc(i))
        sliced = s.get_slice(1, 3)
        assert sliced.N == 2

    def test_make_circular(self):
        s = Strand()
        s.add_nucleotide(_make_nuc(0))
        s.add_nucleotide(_make_nuc(0.5))
        s.make_circular()
        assert s.is_circular()

    def test_cut_in_two(self):
        s = Strand()
        for i in range(10):
            s.add_nucleotide(_make_nuc(i))
        s1, s2 = s.cut_in_two()
        assert s1.N == 5
        assert s2.N == 5

    def test_prepare(self):
        s = Strand()
        s.add_nucleotide(_make_nuc(0))
        s.add_nucleotide(_make_nuc(1))
        s._prepare(3, 100)
        assert s.index == 3
        assert s._first == 100
        assert s._last == 101
        assert s._nucleotides[0].index == 100
        assert s._nucleotides[1].index == 101
