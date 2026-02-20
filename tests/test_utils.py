"""Tests for view_py.utils."""
import math
import numpy as np
import pytest

from view_py.utils import (
    normalize, gram_schmidt, expand_iupac_sequence,
    rotate_vector_around_axis, apply_quaternion,
    quaternion_from_axis_angle, quaternion_conjugate,
    BASE_MAP, COMPLEMENT,
)


class TestNormalize:
    def test_unit_vector(self):
        v = normalize(np.array([3.0, 0, 0]))
        np.testing.assert_allclose(v, [1, 0, 0])

    def test_arbitrary_vector(self):
        v = normalize(np.array([1.0, 1.0, 1.0]))
        assert abs(np.linalg.norm(v) - 1.0) < 1e-10

    def test_zero_vector_unchanged(self):
        v = normalize(np.array([0.0, 0.0, 0.0]))
        np.testing.assert_allclose(v, [0, 0, 0])


class TestGramSchmidt:
    def test_orthonormality(self):
        v1, v2, v3 = gram_schmidt([1, 1, 0], [0, 1, 1], [1, 0, 1])
        assert abs(np.dot(v1, v2)) < 1e-10
        assert abs(np.dot(v1, v3)) < 1e-10
        assert abs(np.dot(v2, v3)) < 1e-10


class TestExpandIupac:
    def test_concrete_bases_unchanged(self):
        assert expand_iupac_sequence("ATCG") == "ATCG"

    def test_n_expands_to_valid_base(self):
        result = expand_iupac_sequence("N")
        assert result in ("A", "C", "G", "T")

    def test_length_preserved(self):
        result = expand_iupac_sequence("NNNNN")
        assert len(result) == 5

    def test_invalid_raises(self):
        with pytest.raises(ValueError):
            expand_iupac_sequence("X")


class TestQuaternions:
    def test_identity_rotation(self):
        q = quaternion_from_axis_angle(np.array([0, 0, 1]), 0)
        v = np.array([1.0, 0, 0])
        result = apply_quaternion(v, q)
        np.testing.assert_allclose(result, v, atol=1e-10)

    def test_90_degree_rotation(self):
        q = quaternion_from_axis_angle(np.array([0, 0, 1]), math.pi / 2)
        v = np.array([1.0, 0, 0])
        result = apply_quaternion(v, q)
        np.testing.assert_allclose(result, [0, 1, 0], atol=1e-10)

    def test_conjugate_reverses(self):
        q = quaternion_from_axis_angle(np.array([0, 0, 1]), math.pi / 3)
        v = np.array([1.0, 2.0, 0])
        rotated = apply_quaternion(v, q)
        restored = apply_quaternion(rotated, quaternion_conjugate(q))
        np.testing.assert_allclose(restored, v, atol=1e-10)


class TestRotateVector:
    def test_90_around_z(self):
        result = rotate_vector_around_axis(
            np.array([1.0, 0, 0]), np.array([0, 0, 1.0]), math.pi / 2)
        np.testing.assert_allclose(result, [0, 1, 0], atol=1e-10)

    def test_360_is_identity(self):
        v = np.array([1.0, 2.0, 3.0])
        result = rotate_vector_around_axis(v, np.array([0, 0, 1.0]), 2 * math.pi)
        np.testing.assert_allclose(result, v, atol=1e-10)


class TestConstants:
    def test_base_map_completeness(self):
        for base in "AaGgCcTtUu":
            assert base in BASE_MAP

    def test_complement_symmetry(self):
        for k, v in COMPLEMENT.items():
            assert COMPLEMENT[v] == k
