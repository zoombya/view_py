import numpy as np
import math
from .nucleotide import Nucleotide
from .strand import Strand
from .utils import (POS_BACK, BASE_BASE, normalize,
                    quaternion_from_axis_angle, apply_quaternion,
                    quaternion_conjugate, rotate_vector_around_axis)


class StrandGenerator:
    def generate_or_sq(self, n_bp, pos=None, direction=None, perp=None,
                       double=True, rot=0.0, angle=None,
                       lengths=None, begin=None, end=None):
        """Generate a helix with per-base-pair twist angles.

        Args:
            n_bp: number of base pairs
            pos: starting position (np.array)
            direction: helix direction vector
            perp: perpendicular vector (backbone orientation)
            double: if True, generate both strands
            rot: initial rotation angle
            angle: list of twist angles (length n_bp-1), or a single float
            lengths: list of insertion lengths (for skips)
            begin: list of begin indices for insertion segments
            end: list of end indices for insertion segments
        """
        if pos is None:
            pos = np.array([0, 0, 0], dtype=float)
        if direction is None:
            direction = np.array([0, 0, 1], dtype=float)
        if lengths is None:
            lengths = []
        if begin is None:
            begin = []
        if end is None:
            end = []

        # Validate/fix length arrays
        if lengths and len(begin) != len(end):
            if len(end) + 1 == len(begin):
                print(f"WARNING: begin ({len(begin)}) and end ({len(end)}) "
                      f"array lengths mismatched; appending n_bp+1 to end")
                end.append(n_bp + 1)
            else:
                raise ValueError(
                    f"begin ({len(begin)}) and end ({len(end)}) "
                    f"array lengths unrecoverably mismatched")

        # Handle angle as list or scalar
        if angle is None:
            angle = [33.75 * math.pi / 180] * (n_bp - 1)
        elif isinstance(angle, (int, float)):
            angle = [angle] * n_bp
        elif len(angle) != n_bp - 1:
            if len(angle) >= n_bp:
                angle = angle[:n_bp - 1]
            else:
                raise ValueError(
                    f"incorrect angle array length ({len(angle)}), "
                    f"should be {n_bp - 1}")

        # Normalize direction
        norm = np.linalg.norm(direction)
        if norm < 1e-10:
            direction = np.array([0, 0, 1], dtype=float)
        else:
            direction = direction / norm

        # Generate perpendicular if not given
        if perp is None:
            perp = np.array([np.random.random(), np.random.random(),
                             np.random.random()])
            perp -= direction * np.dot(direction, perp)
            perp = normalize(perp)
        else:
            perp = perp.copy()

        strand1 = Strand()
        a1 = rotate_vector_around_axis(perp, direction, rot)
        cur_pos = pos.copy()
        a3 = direction.copy()

        # Store quaternions for complementary strand
        quats = []

        for i in range(n_bp):
            strand1.add_nucleotide(Nucleotide(
                cur_pos - a1 * POS_BACK,
                a1.copy(), a3.copy(), None
            ))

            if i != n_bp - 1:
                q = quaternion_from_axis_angle(direction, angle[i])
                quats.append(q)
                a1 = normalize(apply_quaternion(a1, q))
                cur_pos = cur_pos + normalize(a3) * BASE_BASE

                # Handle insertions/deletions (loops/skips)
                if lengths:
                    for k in range(len(lengths)):
                        if i >= begin[k] and i < end[k] and lengths[k]:
                            cur_pos += (normalize(a3) * BASE_BASE *
                                        (-lengths[k] / (end[k] - begin[k])))

        if double:
            # Build complementary strand (reverse direction)
            a1 = -normalize(a1)
            a3 = -direction.copy()
            strand2 = Strand()

            for i in range(n_bp):
                nuc1 = strand1._nucleotides[n_bp - i - 1]
                nuc2 = Nucleotide(
                    cur_pos - a1 * POS_BACK,
                    a1.copy(), a3.copy(), None,
                    pair=nuc1
                )
                nuc1.pair = nuc2
                strand2.add_nucleotide(nuc2)

                if i != n_bp - 1:
                    q = quaternion_conjugate(quats.pop())
                    a1 = normalize(apply_quaternion(a1, q))
                    cur_pos = cur_pos + normalize(a3) * BASE_BASE

                    if lengths:
                        for k in range(len(lengths)):
                            idx = n_bp - 2 - i
                            if idx >= begin[k] and idx < end[k] and lengths[k]:
                                cur_pos += (normalize(a3) * BASE_BASE *
                                            (-lengths[k] /
                                             (end[k] - begin[k])))

            return [strand1, strand2]

        return strand1
