#!/usr/bin/env python3
"""
Convert the rail cadnano design to oxView format using view_py.

This example demonstrates programmatic conversion of a cadnano JSON file
(hexagonal lattice) with a custom scaffold sequence into an oxView-compatible
JSON file that can be loaded in oxView (https://sulcgroup.github.io/oxdna-viewer/).

The rail is one component of the linear actuator from:
    Benson, E., Carrascosa, R., Bath, J., Turberfield, A. J.,
    Strategies for Constructing and Operating DNA Origami Linear Actuators.
    Small 2021, 17, 2007704. https://doi.org/10.1002/smll.202007704

Usage:
    python convert_rail.py

Requirements:
    pip install -e .  (from the repository root)
"""
import json
import os
from pathlib import Path

from view_py import convert_cadnano

HERE = Path(__file__).parent


def main():
    # ---- 1. Load cadnano design and scaffold sequence ----
    cadnano_path = HERE / "rail.json"
    scaffold_path = HERE / "rail_scaffold.txt"

    cadnano_json = cadnano_path.read_text()
    scaffold_seq = scaffold_path.read_text().strip()

    print(f"Loaded cadnano design: {cadnano_path.name}")
    print(f"Scaffold sequence: {len(scaffold_seq)} nt")

    # ---- 2. Convert to oxView format ----
    # Lattice type is auto-detected from helix length (42 = 21*2 -> hexagonal).
    # Pass the scaffold sequence so staples get proper complements.
    system = convert_cadnano(
        cadnano_json,
        sequence=scaffold_seq,
    )

    # ---- 3. Inspect the result ----
    result = system.to_oxview_dict()
    strands = result["systems"][0]["strands"]
    n_strands = len(strands)
    n_nucs = sum(len(s["monomers"]) for s in strands)
    n_bp = sum(
        1 for s in strands for m in s["monomers"] if "bp" in m
    )
    n_clusters = len(set(
        m["cluster"] for s in strands for m in s["monomers"]
        if "cluster" in m
    ))

    print(f"\nConversion result:")
    print(f"  Strands:     {n_strands}")
    print(f"  Nucleotides: {n_nucs}")
    print(f"  Base pairs:  {n_bp}")
    print(f"  Clusters:    {n_clusters}")
    print(f"  Box:         {result['box']}")

    # ---- 4. Write oxView JSON ----
    out_path = HERE / "rail.oxview"
    out_path.write_text(json.dumps(result))
    print(f"\nWrote {out_path.name} ({os.path.getsize(out_path) / 1024:.0f} KB)")
    print(f"Open in oxView: https://sulcgroup.github.io/oxdna-viewer/")


if __name__ == "__main__":
    main()
