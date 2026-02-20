# Linear Actuator Example

Convert the two-component linear actuator cadnano designs to oxView format.

The structure is from:
> Benson, E., Carrascosa, R., Bath, J., Turberfield, A. J., Strategies for Constructing and Operating DNA Origami Linear Actuators. Small 2021, 17, 2007704. https://doi.org/10.1002/smll.202007704

## Files

| File | Description |
|------|-------------|
| `slider.json` | Cadnano design for the slider component |
| `rail.json` | Cadnano design for the rail component |
| `slider_scaffold.txt` | Scaffold sequence for the slider (2880 nt) |
| `rail_scaffold.txt` | Scaffold sequence for the rail (8064 nt) |

## Scripts

| Script | Description | Output |
|--------|-------------|--------|
| `convert_slider.py` | Convert the slider only | `slider.oxview` |
| `convert_rail.py` | Convert the rail only | `rail.oxview` |
| `convert_both.py` | Convert both and combine into one file | `linear_actuator.oxview` |

## Usage

First install the package from the repository root:

```bash
pip install -e .
```

Then run any of the conversion scripts:

```bash
cd examples/linear_actuator
python convert_slider.py   # slider only (102 strands, 6340 nt)
python convert_rail.py     # rail only (260 strands, 17080 nt)
python convert_both.py     # both structures in one file (362 strands, 23420 nt)
```

The output `.oxview` files can be loaded directly in [oxView](https://sulcgroup.github.io/oxdna-viewer/).
