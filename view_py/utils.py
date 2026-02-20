import numpy as np
import math
import random

# oxDNA constants
POS_BACK = 0.4 + 0.2  # 0.6, backbone offset
BASE_BASE = 0.3897628551303122  # rise per base pair

# Lattice spacings
SQ_LATTICE_SPACING = 2.6  # square lattice
HE_LATTICE_SPACING = 2.55  # hexagonal lattice

# Base encoding
BASE_MAP = {'A': 0, 'a': 0, 'G': 1, 'g': 1, 'C': 2, 'c': 2,
            'T': 3, 't': 3, 'U': 3, 'u': 3, 'D': 4}
BASE_NAMES = {0: 'A', 1: 'G', 2: 'C', 3: 'T'}

COMPLEMENT = {0: 3, 1: 2, 2: 1, 3: 0}  # A-T, G-C


def random_base():
    return random.randint(0, 3)


def random_bases(n):
    return [random.randint(0, 3) for _ in range(n)]


def expand_iupac_sequence(seq, is_dna=True):
    """Expand IUPAC ambiguity codes to random concrete bases."""
    s_char = 'T' if is_dna else 'U'
    iupac = {
        'R': ['A', 'G'], 'Y': ['C', s_char], 'M': ['A', 'C'],
        'K': ['G', s_char], 'S': ['C', 'G'], 'W': ['A', s_char],
        'H': ['A', 'C', s_char], 'B': ['C', 'G', s_char],
        'V': ['A', 'C', 'G'], 'D': ['A', 'G', s_char],
        'N': ['A', 'C', 'G', s_char],
    }
    result = []
    for c in seq:
        if c.upper() in ('A', 'T', 'U', 'C', 'G'):
            result.append(c)
        elif c.upper() in iupac:
            result.append(random.choice(iupac[c.upper()]))
        else:
            raise ValueError(f"Unknown base code {c}")
    return ''.join(result)


def gram_schmidt(v1, v2, v3):
    """Gram-Schmidt orthogonalization, returns orthonormal basis."""
    v1 = np.array(v1, dtype=float)
    v2 = np.array(v2, dtype=float)
    v3 = np.array(v3, dtype=float)

    n1 = np.dot(v1, v1)
    o = np.dot(v2, v1)
    v2 = v2 - v1 * (o / n1)

    i = np.dot(v3, v1)
    a = np.dot(v3, v2)
    r = np.dot(v2, v2)
    v3 = v3 - v1 * (i / n1) - v2 * (a / r)

    v1 = v1 / n1
    v2 = v2 / r
    v3 = v3 / np.sqrt(np.dot(v3, v3))
    return v1, v2, v3


def normalize(v):
    n = np.linalg.norm(v)
    if n < 1e-10:
        return v
    return v / n


def rotate_vector_around_axis(vec, axis, angle):
    """Rotate vec around axis by angle (radians) using Rodrigues' formula."""
    axis = normalize(axis)
    cos_a = math.cos(angle)
    sin_a = math.sin(angle)
    return (vec * cos_a +
            np.cross(axis, vec) * sin_a +
            axis * np.dot(axis, vec) * (1 - cos_a))


def apply_quaternion(vec, quat):
    """Apply a quaternion rotation [x,y,z,w] to a vector.
    Uses Hamilton product: q * v * q_conjugate."""
    x, y, z, w = quat
    # quaternion * vector (treating vector as pure quaternion)
    t = 2 * np.cross(np.array([x, y, z]), vec)
    return vec + w * t + np.cross(np.array([x, y, z]), t)


def quaternion_from_axis_angle(axis, angle):
    """Create quaternion [x,y,z,w] from axis-angle."""
    axis = normalize(axis)
    half = angle / 2
    s = math.sin(half)
    return np.array([axis[0]*s, axis[1]*s, axis[2]*s, math.cos(half)])


def quaternion_conjugate(q):
    return np.array([-q[0], -q[1], -q[2], q[3]])
