#!/usr/bin/env python3
"""
Convert both rail and slider cadnano designs and write them as a single oxView file.

This example demonstrates converting two cadnano JSON files and combining them
into one oxView JSON file containing both systems.  The structures are the two
components of the linear actuator from:
    Benson, E., Carrascosa, R., Bath, J., Turberfield, A. J.,
    Strategies for Constructing and Operating DNA Origami Linear Actuators.
    Small 2021, 17, 2007704. https://doi.org/10.1002/smll.202007704

Usage:
    python convert_both.py

Requirements:
    pip install -e .  (from the repository root)
"""
import json
import os
from pathlib import Path

from view_py import convert_cadnano

HERE = Path(__file__).parent


def convert_one(name):
    """Convert a single cadnano design and return the System."""
    cadnano_path = HERE / f"{name}.json"
    scaffold_path = HERE / f"{name}_scaffold.txt"

    cadnano_json = cadnano_path.read_text()
    scaffold_seq = scaffold_path.read_text().strip()

    print(f"Converting {cadnano_path.name}  (scaffold: {len(scaffold_seq)} nt)")

    system = convert_cadnano(cadnano_json, sequence=scaffold_seq)

    print(f"  -> {system.N_strands} strands, {system.N} nucleotides")
    return system


def summarise(oxview_dict):
    """Print a short summary of an oxView dict."""
    total_strands = 0
    total_nucs = 0
    total_bp = 0
    for sys in oxview_dict["systems"]:
        for s in sys["strands"]:
            total_strands += 1
            total_nucs += len(s["monomers"])
            total_bp += sum(1 for m in s["monomers"] if "bp" in m)

    print(f"\nCombined result:")
    print(f"  Systems:     {len(oxview_dict['systems'])}")
    print(f"  Strands:     {total_strands}")
    print(f"  Nucleotides: {total_nucs}")
    print(f"  Base pairs:  {total_bp}")
    print(f"  Box:         {oxview_dict['box']}")


def main():
    # ---- 1. Convert both designs ----
    rail = convert_one("rail")
    slider = convert_one("slider")

    # ---- 2. Merge into a single oxView dict ----
    rail_dict = rail.to_oxview_dict()
    slider_dict = slider.to_oxview_dict()

    # Use the larger box to fit both structures
    box = [max(a, b) for a, b in zip(rail_dict["box"], slider_dict["box"])]

    combined = {
        "box": box,
        "systems": [
            rail_dict["systems"][0],
            slider_dict["systems"][0],
        ],
    }
    # Give each system a distinct id
    for i, sys in enumerate(combined["systems"]):
        sys["id"] = i

    summarise(combined)

    # ---- 3. Write combined oxView JSON ----
    out_path = HERE / "linear_actuator.oxview"
    out_path.write_text(json.dumps(combined))
    print(f"\nWrote {out_path.name} ({os.path.getsize(out_path) / 1024:.0f} KB)")
    print(f"Open in oxView: https://sulcgroup.github.io/oxdna-viewer/")


if __name__ == "__main__":
    main()
