"""
Convert cadnano JSON files to oxView format.

This is a port of the cadnano conversion logic from tacoxdna.js.
"""
import json
import math
import numpy as np
import random

from .nucleotide import Nucleotide, reset_nucleotide_index
from .strand import Strand, reset_strand_index
from .system import System
from .strand_generator import StrandGenerator
from .utils import (BASE_BASE, POS_BACK, SQ_LATTICE_SPACING,
                    HE_LATTICE_SPACING, BASE_MAP, expand_iupac_sequence,
                    normalize, rotate_vector_around_axis)


class CadnanoSquare:
    """Represents a [V_0, b_0, V_1, b_1] entry in cadnano."""

    def __init__(self, V_0=-1, b_0=-1, V_1=-1, b_1=-1):
        self.V_0 = V_0
        self.b_0 = b_0
        self.V_1 = V_1
        self.b_1 = b_1

    def type(self, vhelix, pos):
        """Determine the type of this square: 'empty', 'begin', 'end', 'continue'."""
        if self.V_0 == -1 and self.b_0 == -1:
            if self.V_1 == -1 and self.b_1 == -1:
                return 'empty'
            if self.V_1 == vhelix.num and abs(self.b_1 - pos) == 1:
                return 'begin'
        else:
            if self.V_0 == vhelix.num and abs(self.b_0 - pos) == 1:
                if self.V_1 == -1:
                    return 'end'
                if self.V_1 == vhelix.num and abs(self.b_1 - pos) == 1:
                    return 'continue'
                return 'end'
            if self.V_1 == vhelix.num and abs(self.b_1 - pos) == 1:
                return 'begin'
        print("WARNING: unexpected square array")
        return None


class VirtualHelix:
    """Represents a virtual helix from cadnano."""

    def __init__(self):
        self.stap_loop = []
        self.scaf_loop = []
        self.skip = []
        self.loop = []
        self.stap_colors = []
        self.row = 0
        self.col = 0
        self.num = 0
        self.stap = []
        self.scaf = []
        self.cad_index = -1
        self.skiploop_bases = 0

    @property
    def len(self):
        return max(len(self.scaf), len(self.stap))


class Segments:
    """Track begin/end segment positions for a virtual helix."""

    def __init__(self):
        self.begin = []
        self.end = []

    def add_begin(self, val):
        if val not in self.begin:
            self.begin.append(val)

    def add_end(self, val):
        if val not in self.end:
            self.end.append(val)


class TupleMap:
    """Map with tuple keys, mimicking the JS x class."""

    def __init__(self):
        self._data = {}

    def set(self, key, value):
        self._data[tuple(key)] = value

    def get(self, key, default=None):
        return self._data.get(tuple(key), default)

    def has(self, key):
        return tuple(key) in self._data

    def keys(self):
        return self._data.keys()

    def items(self):
        return self._data.items()

    def values(self):
        return self._data.values()

    @property
    def size(self):
        return len(self._data)


class NucMap(TupleMap):
    """Extended TupleMap for tracking nucleotide-to-strand mappings."""

    def __init__(self):
        super().__init__()
        self._scaf = TupleMap()
        self._stap = TupleMap()
        self.nuc_count = 0
        self.strand_count = 0

    def add_scaf(self, k1, k2, strand_id, nuc_ids):
        self._scaf.set([k1, k2], [strand_id, nuc_ids])

    def add_stap(self, k1, k2, strand_id, nuc_ids):
        self._stap.set([k1, k2], [strand_id, nuc_ids])

    def add_scaf_strand(self, strand_id, other, same_strand=False):
        count = 0
        old_size = self._scaf.size
        for (k1, k2), (sid, nuc_ids) in other._scaf.items():
            if sid == strand_id:
                new_ids = [n + self.nuc_count for n in nuc_ids]
                self.add_scaf(k1, k2, self.strand_count, new_ids)
                count += len(nuc_ids)
        self.nuc_count += count
        if self._scaf.size == old_size:
            return 1
        if not same_strand:
            self.strand_count += 1
        return 0

    def add_stap_strand(self, strand_id, other, same_strand=False):
        count = 0
        old_size = self._stap.size
        for (k1, k2), (sid, nuc_ids) in other._stap.items():
            if sid == strand_id:
                new_ids = [n + self.nuc_count for n in nuc_ids]
                self.add_stap(k1, k2, self.strand_count, new_ids)
                count += len(nuc_ids)
        self.nuc_count += count
        if self._stap.size == old_size:
            return 1
        if not same_strand:
            self.strand_count += 1
        return 0

    def add_strand(self, strand_id, other, same_strand=False):
        r1 = self.add_scaf_strand(strand_id, other, same_strand)
        r2 = self.add_stap_strand(strand_id, other, same_strand)
        return r1 and r2


class CadnanoDesign:
    """Collection of virtual helices."""

    def __init__(self):
        self.vhelices = []

    def add_vhelix(self, vh):
        self.vhelices.append(vh)

    def bbox(self):
        rows = [vh.row for vh in self.vhelices]
        cols = [vh.col for vh in self.vhelices]
        lengths = [vh.len for vh in self.vhelices]
        n = SQ_LATTICE_SPACING * (max(rows) - min(rows) + 2)
        o = SQ_LATTICE_SPACING * (max(cols) - min(cols) + 2)
        i = 0.34 * (max(lengths) + 2)
        return 2 * max(n, o, i) * 2


def _parse_cadnano_json(json_str):
    """Parse cadnano JSON string into a CadnanoDesign."""
    data = json.loads(json_str)
    design = CadnanoDesign()

    for vs in data['vstrands']:
        vh = VirtualHelix()
        for key, val in vs.items():
            if key == 'skip':
                vh.skip = [abs(x) for x in val]
            elif key == 'stap':
                vh.stap = [CadnanoSquare(*sq) for sq in val]
            elif key == 'scaf':
                vh.scaf = [CadnanoSquare(*sq) for sq in val]
            elif key == 'stap_colors':
                vh.stap_colors = val
            elif hasattr(vh, key):
                setattr(vh, key, val)

        vh.skiploop_bases = (len(vh.skip) + sum(vh.loop) - sum(vh.skip))
        design.add_vhelix(vh)

    return design


def _get_sq_angles(vhelix):
    """Generate square lattice twist angles (32-element periodic pattern)."""
    n = vhelix.len - 1
    angles = [0.0] * n

    for i in range(n):
        e = i % 32
        if e < 2:
            angles[i] = 28 * math.pi / 180
        elif e == 2:
            angles[i] = 36 * math.pi / 180
        elif e == 3:
            angles[i] = 54.375 * math.pi / 180
        elif e == 4:
            angles[i] = 37 * math.pi / 180
        elif e in (5, 6):
            angles[i] = 27.6666666666666 * math.pi / 180
        elif e == 7:
            angles[i] = 30.6666666666666 * math.pi / 180
        elif e in (8, 9):
            angles[i] = 29.3333333333 * math.pi / 180
        elif e == 10:
            angles[i] = 34.3333333333 * math.pi / 180
        elif e == 11:
            angles[i] = 54.5 * math.pi / 180
        elif e in (12, 13):
            angles[i] = 28.91666666666 * math.pi / 180
        elif e in (14, 15, 16, 17):
            angles[i] = 31.16666666666 * math.pi / 180
        elif e == 18:
            angles[i] = 35.5 * math.pi / 180
        elif e == 19:
            angles[i] = 52 * math.pi / 180
        elif e == 20:
            angles[i] = 35.5 * math.pi / 180
        elif e in (21, 22):
            angles[i] = 27.5 * math.pi / 180
        elif e == 23:
            angles[i] = 35.5 * math.pi / 180
        elif 24 <= e < 27:
            angles[i] = 30 * math.pi / 180
        elif e == 27:
            angles[i] = 52 * math.pi / 180
        elif e == 28:
            angles[i] = 35.5 * math.pi / 180
        else:
            angles[i] = 30.91666666666 * math.pi / 180

    # Fix last angle in each period of 32 to make total 1080 degrees
    total_31 = sum(angles[i] for i in range(min(31, n)))
    for i in range(n):
        if i % 32 == 31:
            angles[i] = 1080 * math.pi / 180 - total_31

    return angles


def _get_he_angles(vhelix):
    """Generate hexagonal lattice twist angles (21-element periodic pattern)."""
    n = vhelix.len - 1
    angles = [0.0] * n

    for i in range(n):
        e = i % 21
        if e == 0:
            angles[i] = 32.571 * math.pi / 180
        elif e == 1:
            angles[i] = 36 * math.pi / 180
        elif e in (1, 2, 3):
            angles[i] = 42 * math.pi / 180
        elif e in (5, 6, 7):
            angles[i] = 29.143 * math.pi / 180
        elif e == 8:
            angles[i] = 32 * math.pi / 180
        elif e in (9, 10):
            angles[i] = 44 * math.pi / 180
        elif e in (12, 13, 14):
            angles[i] = 28.571 * math.pi / 180
        elif e in (16, 17):
            angles[i] = 41.5 * math.pi / 180
        elif e in (19, 20):
            angles[i] = 28.476 * math.pi / 180
        else:
            angles[i] = 720 / 21 * (math.pi / 180)

    # Fix last angle in each period of 21 to make total 720 degrees
    total_20 = sum(angles[i] for i in range(min(20, n)))
    for i in range(n):
        if i % 21 == 20:
            angles[i] = 720 * math.pi / 180 - total_20

    return angles


def _generate_sq_helix(direction, perp, vhelix):
    """Generate a double helix for a square-lattice virtual helix.
    Port of function S() in tacoxdna.js."""
    gen = StrandGenerator()
    angles = _get_sq_angles(vhelix)

    if vhelix.num % 2 == 0:
        r = np.array([vhelix.col * SQ_LATTICE_SPACING,
                       vhelix.row * SQ_LATTICE_SPACING, 0], dtype=float)
        d = direction.copy()
        p = perp.copy()
        rot = 0
        h = angles[:]
    else:
        r = np.array([vhelix.col * SQ_LATTICE_SPACING,
                       vhelix.row * SQ_LATTICE_SPACING,
                       (vhelix.len - 1) * BASE_BASE], dtype=float)
        d = -direction.copy()
        p = -perp.copy()
        rot = -sum(angles) % (2 * math.pi)
        h = angles[::-1]

    result = gen.generate_or_sq(
        vhelix.len, r, d, p, double=True, rot=rot, angle=h)
    return result, angles, r, rot, d, p


def _generate_he_helix(direction, perp, vhelix):
    """Generate a double helix for a hex-lattice virtual helix.
    Port of function q() in tacoxdna.js."""
    gen = StrandGenerator()
    angles = _get_he_angles(vhelix)
    sp = HE_LATTICE_SPACING

    if vhelix.num % 2 == 0:
        r = np.array([vhelix.col * math.sqrt(3) * sp / 2,
                       3 * vhelix.row * sp / 2, 0], dtype=float)
        d = direction.copy()
        p = perp.copy()
        rot = 0
        h = gen.generate_or_sq(
            vhelix.len, r, d, p, double=True, rot=rot, angle=angles)
    else:
        r = np.array([vhelix.col * math.sqrt(3) * sp / 2,
                       3 * vhelix.row * sp / 2 + 1.275,
                       (vhelix.len - 1) * BASE_BASE], dtype=float)
        d = -direction.copy()
        p = -perp.copy()
        rot = -sum(angles) % (2 * math.pi)
        h = gen.generate_or_sq(
            vhelix.len, r, d, p, double=True, rot=rot,
            angle=angles[::-1])

    return h, angles, r, rot, d, p


def _detect_segments(vhelix):
    """Detect segment boundaries for a virtual helix.
    Port of function $() in tacoxdna.js."""
    segments = Segments()
    step = -1 if vhelix.num % 2 == 0 else 1

    for o in range(len(vhelix.scaf)):
        c = o - step  # previous position
        d = o + step  # next position

        # Get types at neighboring positions
        if 0 <= c < len(vhelix.scaf):
            scaf_prev = vhelix.scaf[c].type(vhelix, c)
            stap_prev = vhelix.stap[c].type(vhelix, c)
        else:
            scaf_prev = False
            stap_prev = False

        if 0 <= d < len(vhelix.scaf):
            scaf_next = vhelix.scaf[d].type(vhelix, d)
            stap_next = vhelix.stap[d].type(vhelix, d)
        else:
            scaf_next = False
            stap_next = False

        # Skip positions where both neighbors have matching begin/end types
        if (scaf_prev == stap_prev and scaf_prev in ('begin', 'end')):
            continue
        if (scaf_next == stap_next and scaf_next in ('begin', 'end')):
            continue

        scaf_type = vhelix.scaf[o].type(vhelix, o)
        stap_type = vhelix.stap[o].type(vhelix, o)

        if scaf_type == 'empty':
            if stap_type == 'begin':
                segments.add_end(o)
            elif stap_type == 'end':
                segments.add_begin(o)
        elif scaf_type == 'begin':
            if stap_type == 'empty':
                segments.add_begin(o)
            elif stap_type == 'continue':
                segments.add_begin(o)
                segments.add_end(o - step)
            elif stap_type == 'begin':
                segments.add_begin(o + step)
                segments.add_end(o - step)
            elif stap_type == 'end':
                segments.add_begin(o)
        elif scaf_type == 'end':
            if stap_type == 'empty':
                segments.add_end(o)
            elif stap_type == 'continue':
                segments.add_begin(o + step)
                segments.add_end(o)
            elif stap_type == 'begin':
                segments.add_end(o)
            elif stap_type == 'end':
                segments.add_begin(o + step)
                segments.add_end(o - step)
        elif scaf_type == 'continue':
            if stap_type == 'begin':
                segments.add_begin(o + step)
                segments.add_end(o)
            elif stap_type == 'end':
                segments.add_begin(o)
                segments.add_end(o - step)

    return segments


def _add_slice_strand(system, vhelix, begin_idx, end_idx, segments,
                      strand_type, helix_pos, helix_dir, helix_perp,
                      helix_rot, helix_angles, u):
    """Add a slice of a helix strand to the system.
    Port of function I() in tacoxdna.js.

    u: 0 for scaffold, 1 for staple
    """
    start = 0
    stop = 0

    if (vhelix.num % 2 + u) % 2 == 0:
        # Even direction
        skip_start = 0
        for s in vhelix.skip[:begin_idx]:
            skip_start -= s
        skip_end = 0
        for s in vhelix.skip[:end_idx + 1]:
            skip_end -= s
        loop_start = 0
        for l in vhelix.loop[:begin_idx]:
            loop_start += l
        loop_end = 0
        for l in vhelix.loop[:end_idx + 1]:
            loop_end += l
        start = begin_idx + skip_start + loop_start
        stop = end_idx + 1 + skip_end + loop_end
    else:
        # Odd direction
        skip_end = 0
        for s in vhelix.skip[end_idx:]:
            skip_end -= s
        skip_start = 0
        for s in vhelix.skip[begin_idx + 1:]:
            skip_start -= s
        loop_end = 0
        for l in vhelix.loop[end_idx:]:
            loop_end += l
        loop_start = 0
        for l in vhelix.loop[begin_idx + 1:]:
            loop_start += l
        start = vhelix.len - begin_idx - 1 + skip_start + loop_start
        stop = vhelix.len - end_idx + skip_end + loop_end

    # Build the helix with segments info, accounting for skips/loops
    # Port of the inner function in I()
    m = _generate_helix_with_segments(
        helix_pos, helix_dir, helix_perp, helix_rot,
        helix_angles, vhelix, segments)

    system.add_strand(m[u].get_slice(start, stop))
    return system


def _generate_helix_with_segments(pos, direction, perp, rot,
                                   all_angles, vhelix, segments):
    """Generate a helix accounting for per-segment angle adjustments.
    Port of the inner function of I() in tacoxdna.js."""
    seg_begin = segments.begin[:]
    seg_end = segments.end[:]

    if vhelix.num % 2 == 1:
        seg_begin = segments.begin[::-1]
        seg_end = segments.end[::-1]

    # Build per-segment adjusted angles
    angles_adjusted = all_angles[:]
    seg_skip_counts = []
    seg_avg_angles = []

    for t in range(len(seg_begin)):
        if vhelix.num % 2 == 0:
            s = seg_begin[t]
            p = seg_end[t]
            n = seg_begin[t]
            o = seg_end[t]
        else:
            s = vhelix.len - seg_begin[t] - 1
            p = vhelix.len - seg_end[t] - 1
            # Swap for reversed direction
            s, p = p, s
            n = vhelix.len - seg_begin[t] - 1
            o = vhelix.len - seg_end[t] - 1
            n, o = o, n

        if o - n != 0:
            skip_count = 0
            for sk in vhelix.skip[s:p + 1]:
                skip_count -= sk
            seg_skip_counts.append(skip_count)
            for lp in vhelix.loop[s:p + 1]:
                skip_count += lp  # reuse var for total adjustment

            avg = sum(all_angles[n:o]) / (o - n + skip_count)
            seg_avg_angles.append(avg)

            for idx in range(n, o):
                angles_adjusted[idx] = avg
        else:
            seg_skip_counts.append(0)
            if len(all_angles) > 0:
                avg = sum(all_angles) / len(all_angles)
            else:
                avg = 33.75 * math.pi / 180
            seg_avg_angles.append(avg)

    # Apply skip removals and loop insertions
    skip_removed = 0
    loop_added = 0
    skip_r = 0
    loop_a = 0

    for t in range(len(seg_begin)):
        skip_removed += skip_r
        loop_added += loop_a
        skip_r = 0
        loop_a = 0

        if vhelix.num % 2 == 0:
            s = seg_begin[t]
            p = seg_end[t]
            n = seg_begin[t]
            o = seg_end[t]
        else:
            s = vhelix.len - seg_begin[t] - 1
            p = vhelix.len - seg_end[t] - 1
            s, p = p, s
            n = vhelix.len - seg_begin[t] - 1
            o = vhelix.len - seg_end[t] - 1
            n, o = o, n

        # Remove skipped positions
        for sk in vhelix.skip[s:p + 1]:
            if sk == 1:
                idx = n - skip_removed + loop_added
                if 0 <= idx < len(angles_adjusted):
                    angles_adjusted.pop(idx)
                skip_r += 1

        # Insert loop angles
        for lp in vhelix.loop[s:p + 1]:
            for _ in range(lp):
                idx = n - skip_removed + loop_added
                if 0 <= idx < len(angles_adjusted):
                    angles_adjusted.insert(idx, seg_avg_angles[t])
                loop_a += 1

    # Build the adjusted begin/end arrays
    adj_begin = []
    adj_end = []
    offset = 0
    for t in range(len(seg_begin)):
        if vhelix.num % 2 == 0:
            n = seg_begin[t]
            o = seg_end[t]
        else:
            n = vhelix.len - seg_begin[t] - 1
            o = vhelix.len - seg_end[t] - 1
            n, o = o, n

        adj_begin.append(n + offset)
        adj_end.append(o + offset + seg_skip_counts[t])
        offset += seg_skip_counts[t]

    gen = StrandGenerator()
    return gen.generate_or_sq(
        len(angles_adjusted) + 1, pos, direction, perp,
        double=True, rot=rot, angle=angles_adjusted,
        lengths=seg_skip_counts, begin=adj_begin, end=adj_end)


def _map_nuc_to_strand(vhelix, strand_id, begin_idx, end_idx,
                        nuc_map, map_type):
    """Map nucleotide positions to strand indices.
    Port of function V() in tacoxdna.js.

    map_type: 2 for scaffold, 3 for staple
    """
    total_adjust = 0
    if (vhelix.num % 2 + map_type) % 2 == 0:
        for s in vhelix.skip[begin_idx:end_idx + 1]:
            total_adjust -= s
        for l in vhelix.loop[begin_idx:end_idx + 1]:
            total_adjust += l
    else:
        for s in vhelix.skip[end_idx:begin_idx + 1]:
            total_adjust -= s
        for l in vhelix.loop[end_idx:begin_idx + 1]:
            total_adjust += l

    if (map_type + vhelix.num % 2) % 2 == 0:
        count = end_idx - begin_idx + 1 + total_adjust
    else:
        count = begin_idx + 1 - end_idx + total_adjust

    h = 0
    skip_count = 0
    loop_count = 0
    while h < count:
        if (map_type + vhelix.num % 2) % 2 == 0:
            o = h + begin_idx + skip_count - loop_count
        else:
            o = begin_idx - h - skip_count + loop_count

        if vhelix.skip[o] != 1:
            if (map_type + vhelix.num % 2) % 2 == 0:
                nuc_ids = list(range(h, h + vhelix.loop[o] + 1))
            else:
                nuc_ids = list(range(h + vhelix.loop[o], -1, -1))
                nuc_ids = [h + x for x in range(vhelix.loop[o] + 1)][::-1]

            if map_type == 2:
                nuc_map.add_scaf(vhelix.num, o, strand_id, nuc_ids)
            elif map_type == 3:
                nuc_map.add_stap(vhelix.num, o, strand_id, nuc_ids)

            h += 1 + vhelix.loop[o]
            loop_count += vhelix.loop[o]
        else:
            if map_type == 2:
                nuc_map.add_scaf(vhelix.num, o, strand_id, [])
            elif map_type == 3:
                nuc_map.add_stap(vhelix.num, o, strand_id, [])
            skip_count += 1

    return nuc_map


def _detect_lattice(design):
    """Auto-detect lattice type from virtual helix lengths.

    If the helix length is a multiple of 32 → square lattice.
    If the helix length is a multiple of 21 → hexagonal lattice.
    Uses the first virtual helix for detection.
    """
    if not design.vhelices:
        return 'sq'

    length = design.vhelices[0].len
    if length % 32 == 0:
        return 'sq'
    elif length % 21 == 0:
        return 'he'
    else:
        print(f"WARNING: Could not auto-detect lattice type from helix "
              f"length {length} (not a multiple of 32 or 21). "
              f"Defaulting to square lattice.")
        return 'sq'


def convert_cadnano(json_str, grid=None, sequence=None, box_side=None,
                    default_val='N'):
    """Convert a cadnano JSON string to an oxView System.

    Args:
        json_str: cadnano JSON as a string
        grid: 'sq' for square lattice, 'he' for hexagonal lattice,
              or None to auto-detect from helix length
        sequence: optional scaffold sequence string
        box_side: optional box side length (auto-calculated if None)
        default_val: default base for staples without sequence ('N' = random)

    Returns:
        System object that can produce oxView JSON via to_oxview_string()
    """
    nuc_map_p = NucMap()
    nuc_map_f = NucMap()

    design = _parse_cadnano_json(json_str)

    if grid is None:
        grid = _detect_lattice(design)
        print(f"INFO: Auto-detected {grid} lattice")

    is_sq = grid == 'sq'
    is_he = grid == 'he'
    if not is_sq and not is_he:
        raise ValueError("grid must be 'sq', 'he', or None (auto-detect)")

    if box_side is None:
        box_side = design.bbox()
        print(f"INFO: Using default box size ({box_side:.1f}), "
              f"a factor 2 larger than the cadnano system")

    direction = np.array([0, 0, 1], dtype=float)
    perp = np.array([1, 0, 0], dtype=float)

    if is_sq:
        perp = rotate_vector_around_axis(perp, direction, 15 * math.pi / 180)
    elif is_he:
        perp = rotate_vector_around_axis(perp, direction, 160 * math.pi / 180)

    # Build temporary system for helix generation
    system_R = System(np.array([box_side, box_side, box_side]))
    # Final reversed system
    system_z = System(np.array([box_side, box_side, box_side]))

    cad_index = 0
    strand_W = -1  # global strand counter

    # Crossover tracking for scaffold
    scaf_crossover_targets = []
    scaf_crossover_chains = []
    # Crossover tracking for staple
    stap_crossover_targets = []
    stap_crossover_chains = []

    found_crossover = False

    begin_F = -1
    end_L = -1

    for vh in design.vhelices:
        vh.cad_index = cad_index

        if is_sq:
            helix_result, h_angles, h_pos, h_rot, h_dir, h_perp = \
                _generate_sq_helix(direction, perp, vh)
        else:
            helix_result, h_angles, h_pos, h_rot, h_dir, h_perp = \
                _generate_he_helix(direction, perp, vh)

        segments = _detect_segments(vh)

        # Process scaffold strand
        c = 0
        for sq in vh.scaf:
            if sq.V_0 == -1 and sq.b_0 == -1:
                if not (sq.V_1 == -1 and sq.b_1 == -1):
                    if sq.V_1 == vh.num and abs(sq.b_1 - c) == 1:
                        # Begin of scaffold
                        if vh.num % 2 == 0:
                            strand_W += 1
                        begin_F = c
                        if vh.num % 2 == 1:
                            system_R = _add_slice_strand(
                                system_R, vh, begin_F, end_L, segments, 0,
                                h_pos, h_dir, h_perp, h_rot, h_angles, 0)
                            nuc_map_p = _map_nuc_to_strand(
                                vh, strand_W, begin_F, end_L, nuc_map_p, 2)
            elif sq.V_0 == vh.num and abs(sq.b_0 - c) == 1:
                if sq.V_1 == -1 and sq.b_1 == -1:
                    # End of scaffold
                    if vh.num % 2 == 1:
                        strand_W += 1
                    end_L = c
                    if vh.num % 2 == 0:
                        system_R = _add_slice_strand(
                            system_R, vh, begin_F, end_L, segments, 0,
                            h_pos, h_dir, h_perp, h_rot, h_angles, 0)
                        nuc_map_p = _map_nuc_to_strand(
                            vh, strand_W, begin_F, end_L, nuc_map_p, 2)
                elif sq.V_1 == vh.num and abs(sq.b_1 - c) == 1:
                    pass  # Continue
                else:
                    # Crossover out
                    if vh.num % 2 == 1:
                        strand_W += 1
                    end_L = c
                    if vh.num % 2 == 0:
                        system_R = _add_slice_strand(
                            system_R, vh, begin_F, end_L, segments, 0,
                            h_pos, h_dir, h_perp, h_rot, h_angles, 0)
                        nuc_map_p = _map_nuc_to_strand(
                            vh, strand_W, begin_F, end_L, nuc_map_p, 2)

                    target = [sq.V_1, sq.b_1]
                    d = c
                    found_crossover = False
                    for e in range(len(scaf_crossover_targets)):
                        if [vh.num, d] == scaf_crossover_targets[e]:
                            scaf_crossover_chains[e].insert(0, strand_W)
                            found_crossover = True
                    if not found_crossover:
                        scaf_crossover_chains.append([strand_W])
                        scaf_crossover_targets.append(target)
                    found_crossover = False
            elif sq.V_1 == vh.num and abs(sq.b_1 - c) == 1:
                # Crossover in (begin)
                if vh.num % 2 == 0:
                    strand_W += 1
                begin_F = c
                if vh.num % 2 == 1:
                    system_R = _add_slice_strand(
                        system_R, vh, begin_F, end_L, segments, 0,
                        h_pos, h_dir, h_perp, h_rot, h_angles, 0)
                    nuc_map_p = _map_nuc_to_strand(
                        vh, strand_W, begin_F, end_L, nuc_map_p, 2)

                source = [sq.V_0, sq.b_0]
                found_crossover = False
                for e in range(len(scaf_crossover_targets)):
                    if [vh.num, c] == scaf_crossover_targets[e]:
                        scaf_crossover_chains[e].append(strand_W)
                        found_crossover = True
                if not found_crossover:
                    scaf_crossover_chains.append([strand_W])
                    scaf_crossover_targets.append(source)
                found_crossover = False
            c += 1

        if system_R.N_strands == 0:
            print(f"WARNING: No scaffold strand found in virtual helix "
                  f"n. {vh.num}: staples-only virtual helices are not "
                  f"supported when no scaffold has been seen yet")
            continue

        # Process staple strand
        c = 0
        for sq in vh.stap:
            if sq.V_0 == -1 and sq.b_0 == -1:
                if not (sq.V_1 == -1 and sq.b_1 == -1):
                    if sq.V_1 == vh.num and abs(sq.b_1 - c) == 1:
                        if vh.num % 2 == 1:
                            strand_W += 1
                        begin_F = c
                        if vh.num % 2 == 0:
                            system_R = _add_slice_strand(
                                system_R, vh, begin_F, end_L, segments, 0,
                                h_pos, h_dir, h_perp, h_rot, h_angles, 1)
                            nuc_map_p = _map_nuc_to_strand(
                                vh, strand_W, begin_F, end_L, nuc_map_p, 3)
            elif sq.V_0 == vh.num and abs(sq.b_0 - c) == 1:
                if sq.V_1 == -1 and sq.b_1 == -1:
                    if vh.num % 2 == 0:
                        strand_W += 1
                    end_L = c
                    if vh.num % 2 == 1:
                        system_R = _add_slice_strand(
                            system_R, vh, begin_F, end_L, segments, 0,
                            h_pos, h_dir, h_perp, h_rot, h_angles, 1)
                        nuc_map_p = _map_nuc_to_strand(
                            vh, strand_W, begin_F, end_L, nuc_map_p, 3)
                elif sq.V_1 == vh.num and abs(sq.b_1 - c) == 1:
                    pass  # Continue
                else:
                    # Crossover
                    if vh.num % 2 == 0:
                        strand_W += 1
                    end_L = c
                    if vh.num % 2 == 1:
                        system_R = _add_slice_strand(
                            system_R, vh, begin_F, end_L, segments, 0,
                            h_pos, h_dir, h_perp, h_rot, h_angles, 1)
                        nuc_map_p = _map_nuc_to_strand(
                            vh, strand_W, begin_F, end_L, nuc_map_p, 3)

                    d = c
                    found_crossover = False
                    for e in range(len(stap_crossover_targets)):
                        if [vh.num, d] == stap_crossover_targets[e]:
                            stap_crossover_chains[e].insert(0, strand_W)
                            found_crossover = True
                    if not found_crossover:
                        stap_crossover_chains.append([strand_W])
                        stap_crossover_targets.append([sq.V_1, sq.b_1])
                    found_crossover = False
            elif sq.V_1 == vh.num and abs(sq.b_1 - c) == 1:
                if vh.num % 2 == 1:
                    strand_W += 1
                begin_F = c
                if vh.num % 2 == 0:
                    system_R = _add_slice_strand(
                        system_R, vh, begin_F, end_L, segments, 0,
                        h_pos, h_dir, h_perp, h_rot, h_angles, 1)
                    nuc_map_p = _map_nuc_to_strand(
                        vh, strand_W, begin_F, end_L, nuc_map_p, 3)

                found_crossover = False
                for e in range(len(stap_crossover_targets)):
                    if [vh.num, c] == stap_crossover_targets[e]:
                        stap_crossover_chains[e].append(strand_W)
                        found_crossover = True
                if not found_crossover:
                    stap_crossover_chains.append([strand_W])
                    stap_crossover_targets.append([sq.V_0, sq.b_0])
                found_crossover = False
            c += 1

        cad_index += 1

    # Merge crossover chains
    all_chains = [scaf_crossover_chains, stap_crossover_chains]
    in_chain = set()
    for chain_group in all_chains:
        for chain in chain_group:
            in_chain.update(chain)

    # Add non-crossover strands to final system
    for t in range(len(system_R._strands)):
        if t not in in_chain:
            system_z.add_strand(system_R._strands[t])
            nuc_map_f.add_strand(t, nuc_map_p)

    # Process crossover chains (join strands)
    for chain_type in range(2):
        chains = all_chains[chain_type]
        circular_chains = []

        # Merge chains that connect
        done = False
        while not done:
            merged = False
            for i in range(len(chains)):
                if merged:
                    break
                for j in range(len(chains)):
                    if merged:
                        break
                    if i != j and chains[i][0] == chains[j][-1]:
                        chains[j].extend(chains[i][1:])
                        chains.pop(i)
                        merged = True
                    elif i == j and chains[i][0] == chains[i][-1]:
                        if i not in circular_chains:
                            circular_chains.append(i)
            if not merged:
                done = True

        # Build final strands from chains
        for i, chain in enumerate(chains):
            strand = system_R._strands[chain[0]]
            is_circ = i in circular_chains
            if is_circ:
                for k in range(1, len(chain) - 1):
                    strand = strand.append(system_R._strands[chain[k]])
                strand.make_circular(True)
            else:
                for k in range(1, len(chain)):
                    strand = strand.append(system_R._strands[chain[k]])

            system_z.add_strand(strand)

            if strand._circular:
                idx_range = list(range(len(chain) - 2))
            else:
                idx_range = list(range(len(chain) - 1))

            for k in idx_range:
                nuc_map_f.add_strand(chain[k], nuc_map_p, same_strand=True)
            last_idx = idx_range[-1] + 1 if len(idx_range) > 0 else 0
            nuc_map_f.add_strand(chain[last_idx], nuc_map_p,
                                 same_strand=False)

    # Now reverse all strands to get proper 3'→5' ordering
    # Port of the "Q" system construction in the JS
    system_Q = System(system_z._box)
    for strand in system_z._strands:
        nucs_reversed = list(reversed(strand._nucleotides))
        new_strand = Strand()
        for nuc in nucs_reversed:
            new_nuc = Nucleotide(
                nuc.cm_pos, nuc._a1, -nuc._a3, nuc._base)
            new_strand.add_nucleotide(new_nuc)
        if strand._circular:
            new_strand.make_circular()
        system_Q.add_strand(new_strand)

    # Rebuild nucleotide index map for the reversed system
    system_Q._prepare()

    # Build nuc_map for reversed system
    nuc_map_J = NucMap()
    strand_starts = []
    k = 0
    for strand in system_Q._strands:
        strand_starts.append(k)
        k += strand.N

    # Map from old nuc indices to reversed indices
    nuc_to_pos = {}
    for key, (sid, nuc_ids) in nuc_map_f._scaf.items():
        new_ids = []
        for nid in nuc_ids:
            s = system_Q._strands[sid]
            new_idx = s.N - 1 - (nid - strand_starts[sid]) + strand_starts[sid]
            new_ids.append(new_idx)
            nuc_to_pos[new_idx] = key
        nuc_map_J.add_scaf(key[0], key[1], sid, new_ids)

    for key, (sid, nuc_ids) in nuc_map_f._stap.items():
        new_ids = []
        for nid in nuc_ids:
            s = system_Q._strands[sid]
            new_idx = s.N - 1 - (nid - strand_starts[sid]) + strand_starts[sid]
            new_ids.append(new_idx)
            nuc_to_pos[new_idx] = key
        nuc_map_J.add_stap(key[0], key[1], sid, new_ids)

    if system_Q.N == 0:
        print("CRITICAL: The generated configuration is empty.")
        return system_Q

    # Build stap color map
    stap_color_map = TupleMap()
    for vh in design.vhelices:
        for pos, color in vh.stap_colors:
            stap_color_map.set([vh.num, pos], color)

    # Map position -> nucleotide list
    pos_to_nucs = TupleMap()
    for nuc_idx, pos_key in nuc_to_pos.items():
        if not pos_to_nucs.has(pos_key):
            pos_to_nucs.set(pos_key, [])
        pos_to_nucs.get(pos_key).append(nuc_idx)

    # Build index -> nucleotide lookup
    nuc_by_index = {}
    for strand in system_Q._strands:
        for nuc in strand._nucleotides:
            nuc_by_index[nuc.index] = nuc

    # Find scaffold strand (most common strand among paired positions)
    strand_counts = {}
    for nuc_list in pos_to_nucs.values():
        for nuc_idx in nuc_list:
            nuc = nuc_by_index.get(nuc_idx)
            if nuc:
                sid = nuc.strand
                strand_counts[sid] = strand_counts.get(sid, 0) + 1

    scaffold_strand_id = max(strand_counts, key=strand_counts.get)

    # Set base pairs and colors
    for strand in system_Q._strands:
        for nuc in strand._nucleotides:
            if nuc.index in nuc_to_pos:
                pos_key = nuc_to_pos[nuc.index]
                nuc_list = pos_to_nucs.get(pos_key)
                if nuc_list and len(nuc_list) > 1:
                    partner_idx = [x for x in nuc_list if x != nuc.index]
                    if partner_idx:
                        partner = nuc_by_index.get(partner_idx[0])
                        if partner:
                            nuc.pair = partner

                # Assign colors to staple nucleotides
                staple_nucs = [x for x in nuc_list
                               if nuc_by_index.get(x) and
                               nuc_by_index[x].strand != scaffold_strand_id]
                if staple_nucs:
                    staple_nuc = nuc_by_index.get(staple_nucs[0])
                    if nuc is staple_nuc:
                        if stap_color_map.has(pos_key):
                            nuc.color = stap_color_map.get(pos_key)
                    else:
                        nuc.color = 3633362  # default scaffold color

    # Propagate colors to whole strands
    for strand in system_Q._strands:
        strand_color = None
        for nuc in strand._nucleotides:
            if nuc.color is not None:
                strand_color = nuc.color
                break
        if strand_color is not None:
            for nuc in strand._nucleotides:
                nuc.color = strand_color

    # Compute clusters
    system_Q.calc_clusters()

    # Apply scaffold sequence
    scaffold = system_Q._strands[scaffold_strand_id]
    seq_array = None

    if sequence:
        if len(sequence) < scaffold.N:
            print(f"WARNING: Provided scaffold sequence is {len(sequence)}nt "
                  f"but needs to be at least {scaffold.N}")
        else:
            print("Applying custom sequence")
            seq_array = [BASE_MAP[c] for c in sequence]

    if seq_array is None:
        print("Applying random sequence")
        seq_array = [random.randint(0, 3) for _ in range(scaffold.N)]

    # Apply sequence (reversed order, matching JS)
    for i in range(scaffold.N):
        nuc = scaffold._nucleotides[i]
        nuc._base = seq_array[scaffold.N - i - 1]
        if nuc.pair:
            nuc.pair._base = 3 - nuc._base  # complement

    # Apply random/default bases to unpaired staple nucleotides
    for i, strand in enumerate(system_Q._strands):
        if i != scaffold_strand_id:
            for nuc in strand._nucleotides:
                if nuc.pair is None:
                    expanded = expand_iupac_sequence(default_val, is_dna=True)
                    nuc._base = BASE_MAP[expanded]

    return system_Q
