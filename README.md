# view_py

A Python package for converting [cadnano](https://cadnano.org/) DNA origami designs into [oxView](https://github.com/sulcgroup/oxdna-viewer) format.

This is a standalone port of the cadnano conversion logic from the [tacoxDNA](http://tacoxdna.sissa.it/) JavaScript library used by oxView, enabling programmatic conversion without a browser.

## Installation

From the repository root:

```bash
pip install .
```

Or in development/editable mode:

```bash
pip install -e .
```

The only dependency is **numpy**.

## Quick Start

```python
from view_py import convert_cadnano

# Read a cadnano JSON file
with open("design.json") as f:
    cadnano_json = f.read()

# Convert to oxView format (lattice auto-detected, with scaffold sequence)
system = convert_cadnano(cadnano_json, sequence="ATCGATCG...")

# Write oxView JSON
with open("design.oxview", "w") as f:
    f.write(system.to_oxview_string())
```

## API Reference

### `convert_cadnano(json_str, grid, sequence, box_side, default_val)`

Main entry point. Converts a cadnano JSON string into a `System` object.

| Parameter     | Type         | Default | Description                                                     |
|---------------|--------------|---------|-----------------------------------------------------------------|
| `json_str`    | `str`        | —       | cadnano JSON file contents                                      |
| `grid`        | `str \| None`| `None`  | Lattice type: `'sq'` (square), `'he'` (hexagonal), or `None` to auto-detect from helix length (multiple of 32 → square, multiple of 21 → hexagonal) |
| `sequence`    | `str \| None`| `None`  | Scaffold sequence. If `None`, a random sequence is assigned.    |
| `box_side`    | `float \| None`| `None`| Simulation box side length. Auto-calculated if `None`.          |
| `default_val` | `str`        | `'N'`   | Default base for unsequenced staples (IUPAC code, `'N'`=random) |

**Returns:** `System` object.

### `System`

Represents a complete oxDNA system (strands in a simulation box).

| Method / Property   | Description                                        |
|---------------------|----------------------------------------------------|
| `to_oxview_string()`| Returns oxView JSON as a string                    |
| `to_oxview_dict()`  | Returns oxView JSON as a Python dict               |
| `N`                 | Total number of nucleotides                        |
| `N_strands`         | Total number of strands                            |
| `isDNA`             | `True` for DNA (default), `False` for RNA          |

## Supported Lattice Types

- **Square lattice** (`grid='sq'`): 32-element periodic twist angle pattern, spacing 2.6 oxDNA units
- **Hexagonal lattice** (`grid='he'`): 21-element periodic twist angle pattern, spacing 2.55 oxDNA units

## Output Format

The output follows the [oxView file format](https://github.com/sulcgroup/oxdna-viewer/blob/master/file-format.md) specification:

```json
{
  "box": [100, 100, 100],
  "systems": [{
    "id": 0,
    "strands": [{
      "id": 0,
      "end3": 0,
      "end5": 63,
      "class": "NucleicAcidStrand",
      "monomers": [{
        "id": 0, "type": "A", "class": "DNA",
        "p": [0, 0, 0],
        "a1": [1, 0, 0],
        "a3": [0, 0, 1],
        "n5": 1, "bp": 95,
        "cluster": 1, "color": 3633362
      }]
    }]
  }]
}
```

Each monomer includes:
- **id**: unique index
- **type**: base (A, T, C, G)
- **class**: `"DNA"` or `"RNA"`
- **p**: center-of-mass position (oxDNA units, 1 unit = 0.8518 nm)
- **a1**: backbone-base unit vector
- **a3**: stacking unit vector
- **n3/n5**: 3' and 5' neighbor indices
- **bp**: base pair partner index (if paired)
- **cluster**: cluster group index
- **color**: base-10 color value (from cadnano staple colors)

## Conversion Pipeline

1. Parse cadnano JSON and extract `vstrands` array
2. Generate 3D helix coordinates using square or hexagonal lattice twist angles
3. Walk scaffold/staple arrays to identify strand segments and crossovers
4. Join segments across crossovers into continuous strands
5. Reverse strands to establish proper 3'→5' ordering
6. Apply scaffold sequence (user-provided or random)
7. Assign Watson-Crick complements and staple colors
8. Compute connectivity clusters based on backbone distance
9. Output oxView JSON

## Limitations

- Only cadnano format is supported (not rpoly, tiamat, or PDB)
- Staple-only virtual helices (no scaffold) are not supported
- Structures should be relaxed with molecular dynamics after conversion, as crossover backbone distances may be non-physical
