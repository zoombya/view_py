"""Tests for the oxview_py CLI."""
import json
from pathlib import Path

import pytest

from view_py.cli import main

EXAMPLES = Path(__file__).parent.parent / "examples" / "linear_actuator"


class TestCLIConvert:
    def test_basic_convert(self, tmp_path):
        out = tmp_path / "out.oxview"
        main(["convert", str(EXAMPLES / "slider.json"), "-o", str(out)])
        assert out.exists()
        data = json.loads(out.read_text())
        assert 'box' in data
        strands = data['systems'][0]['strands']
        assert len(strands) == 102
        total = sum(len(s['monomers']) for s in strands)
        assert total == 6340

    def test_default_output_name(self, tmp_path, monkeypatch):
        # Copy slider.json to tmp_path so .oxview lands there too
        src = EXAMPLES / "slider.json"
        dst = tmp_path / "slider.json"
        dst.write_text(src.read_text())
        main(["convert", str(dst)])
        out = tmp_path / "slider.oxview"
        assert out.exists()

    def test_with_sequence_file(self, tmp_path):
        out = tmp_path / "out.oxview"
        main([
            "convert", str(EXAMPLES / "slider.json"),
            "-s", str(EXAMPLES / "slider_scaffold.txt"),
            "-o", str(out),
        ])
        data = json.loads(out.read_text())
        strands = data['systems'][0]['strands']
        total = sum(len(s['monomers']) for s in strands)
        assert total == 6340

    def test_with_grid_option(self, tmp_path):
        out = tmp_path / "out.oxview"
        main([
            "convert", str(EXAMPLES / "slider.json"),
            "-g", "he",
            "-o", str(out),
        ])
        assert out.exists()

    def test_with_box_option(self, tmp_path):
        out = tmp_path / "out.oxview"
        main([
            "convert", str(EXAMPLES / "slider.json"),
            "-b", "999",
            "-o", str(out),
        ])
        data = json.loads(out.read_text())
        assert data['box'] == [999, 999, 999]

    def test_no_command_exits(self):
        with pytest.raises(SystemExit):
            main([])

    def test_rail_convert(self, tmp_path):
        out = tmp_path / "rail.oxview"
        main(["convert", str(EXAMPLES / "rail.json"), "-o", str(out)])
        data = json.loads(out.read_text())
        strands = data['systems'][0]['strands']
        assert len(strands) == 260
        total = sum(len(s['monomers']) for s in strands)
        assert total == 17080
