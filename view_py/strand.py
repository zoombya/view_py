import numpy as np
import math

_global_strand_index = 0


def reset_strand_index():
    global _global_strand_index
    _global_strand_index = 0


class Strand:
    def __init__(self):
        global _global_strand_index
        self.index = _global_strand_index
        _global_strand_index += 1
        self._first = -1
        self._last = -1
        self._nucleotides = []
        self._sequence = []
        self._circular = False

    @property
    def N(self):
        return len(self._nucleotides)

    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self, seq):
        from .utils import BASE_MAP
        if isinstance(seq, str):
            seq = [BASE_MAP[c] for c in seq]
        if len(seq) == len(self._nucleotides):
            for i, nuc in enumerate(self._nucleotides):
                nuc._base = seq[i]
            self._sequence = seq
        else:
            print(f"WARNING: Cannot change sequence: lengths don't match "
                  f"({len(seq)} vs {len(self._nucleotides)})")

    def _prepare(self, strand_idx, nuc_start):
        self.index = strand_idx
        self._first = nuc_start
        for i, nuc in enumerate(self._nucleotides):
            nuc.index = nuc_start + i
        self._last = nuc_start + len(self._nucleotides) - 1
        return nuc_start + len(self._nucleotides)

    def add_nucleotide(self, nuc):
        if len(self._nucleotides) == 0:
            self._first = nuc.index
        nuc.strand = self.index
        self._nucleotides.append(nuc)
        self._last = nuc.index
        self._sequence.append(nuc._base)

    def append(self, other):
        """Concatenate two strands, return a new strand."""
        s = Strand()
        for nuc in self._nucleotides:
            s.add_nucleotide(nuc)
        for nuc in other._nucleotides:
            s.add_nucleotide(nuc)
        return s

    def get_slice(self, start=0, end=None):
        if end is None:
            end = self.N
        s = Strand()
        for i in range(start, end):
            s.add_nucleotide(self._nucleotides[i].copy())
        return s

    def make_circular(self, check=False):
        if check:
            diff = self._nucleotides[-1].distance(self._nucleotides[0])
            dist = math.sqrt(np.dot(diff, diff))
            if dist > 1.0025:
                print("WARNING: Strand.make_circular(): ends seem too far apart.")
        self._circular = True

    def is_circular(self):
        return self._circular

    def cut_in_two(self, do_copy=True):
        s1 = Strand()
        s2 = Strand()
        half = len(self._nucleotides) // 2
        for i, nuc in enumerate(self._nucleotides):
            target = s1 if i < half else s2
            target.add_nucleotide(nuc.copy() if do_copy else nuc)
        return s1, s2
