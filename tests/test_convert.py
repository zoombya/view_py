"""Tests for cadnano conversion using the linear actuator example data."""
import json
from pathlib import Path

import numpy as np
import pytest

from view_py import convert_cadnano

EXAMPLES = Path(__file__).parent.parent / "examples" / "linear_actuator"


@pytest.fixture
def slider_json():
    return (EXAMPLES / "slider.json").read_text()


@pytest.fixture
def slider_sequence():
    return (EXAMPLES / "slider_scaffold.txt").read_text().strip()


@pytest.fixture
def rail_json():
    return (EXAMPLES / "rail.json").read_text()


@pytest.fixture
def rail_sequence():
    return (EXAMPLES / "rail_scaffold.txt").read_text().strip()


class TestSliderConversion:
    def test_strand_count(self, slider_json):
        system = convert_cadnano(slider_json)
        assert system.N_strands == 102

    def test_nucleotide_count(self, slider_json):
        system = convert_cadnano(slider_json)
        assert system.N == 6340

    def test_with_sequence(self, slider_json, slider_sequence):
        system = convert_cadnano(slider_json, sequence=slider_sequence)
        assert system.N == 6340
        assert system.N_strands == 102

    def test_explicit_he_grid(self, slider_json):
        system = convert_cadnano(slider_json, grid='he')
        assert system.N == 6340

    def test_oxview_dict_valid(self, slider_json):
        system = convert_cadnano(slider_json)
        result = system.to_oxview_dict()
        assert 'box' in result
        assert 'systems' in result
        strands = result['systems'][0]['strands']
        assert len(strands) == 102

    def test_oxview_nucleotide_count(self, slider_json):
        system = convert_cadnano(slider_json)
        result = system.to_oxview_dict()
        strands = result['systems'][0]['strands']
        total = sum(len(s['monomers']) for s in strands)
        assert total == 6340

    def test_pair_count(self, slider_json):
        system = convert_cadnano(slider_json)
        result = system.to_oxview_dict()
        strands = result['systems'][0]['strands']
        n_paired = sum(
            1 for s in strands for m in s['monomers'] if 'pair' in m)
        assert n_paired > 0
        assert n_paired % 2 == 0  # pairs are reciprocal

    def test_pair_reciprocity(self, slider_json):
        system = convert_cadnano(slider_json)
        result = system.to_oxview_dict()
        strands = result['systems'][0]['strands']
        id_to_monomer = {}
        for s in strands:
            for m in s['monomers']:
                id_to_monomer[m['id']] = m

        for m in id_to_monomer.values():
            if 'pair' in m:
                partner = id_to_monomer[m['pair']]
                assert 'pair' in partner
                assert partner['pair'] == m['id']

    def test_clusters_assigned(self, slider_json):
        system = convert_cadnano(slider_json)
        result = system.to_oxview_dict()
        strands = result['systems'][0]['strands']
        has_cluster = any(
            'cluster' in m for s in strands for m in s['monomers'])
        assert has_cluster

    def test_strand_classes(self, slider_json):
        system = convert_cadnano(slider_json)
        result = system.to_oxview_dict()
        strands = result['systems'][0]['strands']
        for s in strands:
            assert s['class'] == 'NucleicAcidStrand'
            for m in s['monomers']:
                assert m['class'] == 'DNA'

    def test_sequence_complementarity(self, slider_json, slider_sequence):
        system = convert_cadnano(slider_json, sequence=slider_sequence)
        result = system.to_oxview_dict()
        strands = result['systems'][0]['strands']
        id_to_monomer = {}
        for s in strands:
            for m in s['monomers']:
                id_to_monomer[m['id']] = m

        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        for m in id_to_monomer.values():
            if 'pair' in m:
                partner = id_to_monomer[m['pair']]
                assert m['type'] == complement[partner['type']], (
                    f"Nucleotide {m['id']} ({m['type']}) paired with "
                    f"{partner['id']} ({partner['type']})")

    def test_oxview_string_valid_json(self, slider_json):
        system = convert_cadnano(slider_json)
        parsed = json.loads(system.to_oxview_string())
        assert 'box' in parsed

    def test_custom_box(self, slider_json):
        system = convert_cadnano(slider_json, box_side=500.0)
        result = system.to_oxview_dict()
        assert result['box'] == [500, 500, 500]


class TestRailConversion:
    def test_strand_count(self, rail_json):
        system = convert_cadnano(rail_json)
        assert system.N_strands == 260

    def test_nucleotide_count(self, rail_json):
        system = convert_cadnano(rail_json)
        assert system.N == 17080

    def test_with_sequence(self, rail_json, rail_sequence):
        system = convert_cadnano(rail_json, sequence=rail_sequence)
        assert system.N == 17080
        assert system.N_strands == 260


class TestCombined:
    def test_total_nucleotides(self, slider_json, rail_json):
        slider = convert_cadnano(slider_json)
        rail = convert_cadnano(rail_json)
        assert slider.N + rail.N == 23420


class TestLatticeDetection:
    def test_slider_autodetects_hexagonal(self, slider_json):
        # Should not raise; hexagonal is auto-detected
        system = convert_cadnano(slider_json)
        assert system.N == 6340

    def test_rail_autodetects_hexagonal(self, rail_json):
        system = convert_cadnano(rail_json)
        assert system.N == 17080


class TestEdgeCases:
    def test_random_sequence_when_none(self, slider_json):
        s1 = convert_cadnano(slider_json)
        s2 = convert_cadnano(slider_json)
        # Both should produce valid systems (sequences may differ)
        assert s1.N == s2.N == 6340

    def test_is_dna(self, slider_json):
        system = convert_cadnano(slider_json)
        assert system.isDNA is True
