import numpy as np
from .utils import BASE_NAMES, POS_BACK

_global_index = 0


def reset_nucleotide_index():
    global _global_index
    _global_index = 0


class Nucleotide:
    def __init__(self, cm_pos, a1, a3, base=None, v=None, L=None,
                 n3=-1, pair=None, cluster=None, color=None):
        global _global_index
        self.index = _global_index
        _global_index += 1

        self.cm_pos = np.array(cm_pos, dtype=float)
        self._a1 = np.array(a1, dtype=float)
        norm = np.linalg.norm(self._a1)
        if norm > 1e-10:
            self._a1 /= norm
        self._a3 = np.array(a3, dtype=float)
        norm = np.linalg.norm(self._a3)
        if norm > 1e-10:
            self._a3 /= norm

        if base is None:
            import random
            base = random.randint(0, 3)
        if isinstance(base, str):
            from .utils import BASE_MAP
            base = BASE_MAP.get(base, 0)
        self._base = base

        self._v = np.array(v if v is not None else [0, 0, 0], dtype=float)
        self._L = np.array(L if L is not None else [0, 0, 0], dtype=float)
        self.n3 = n3
        self.next = -1
        self.pair = pair
        self.cluster = cluster
        self.color = color
        self.strand = None

    @property
    def pos_base(self):
        return self.cm_pos + self._a1 * 0.4

    @property
    def pos_back(self):
        return self.cm_pos + self._a1 * (-0.4)

    @property
    def a2(self):
        return np.cross(self._a3, self._a1)

    def get_base(self):
        if self._base in BASE_NAMES:
            return BASE_NAMES[self._base]
        return str(self._base)

    def distance(self, other, pbc=False, box=None):
        """Backbone-backbone distance vector."""
        diff = other.pos_back - self.pos_back
        if pbc and box is not None:
            diff -= box * np.round(diff / box)
        return diff

    def copy(self):
        nuc = Nucleotide(
            self.cm_pos.copy(), self._a1.copy(), self._a3.copy(),
            self._base, self._v.copy(), self._L.copy(),
            self.n3, self.pair, self.cluster, self.color
        )
        return nuc
