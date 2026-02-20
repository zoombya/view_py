"""Command-line interface for view_py."""
import argparse
import json
import sys
from pathlib import Path

from .cadnano_reader import convert_cadnano


def cmd_convert(args):
    """Run the convert subcommand."""
    cadnano_json = args.input.read_text()

    sequence = None
    if args.sequence:
        seq_path = Path(args.sequence)
        if seq_path.is_file():
            sequence = seq_path.read_text().strip()
        else:
            sequence = args.sequence

    system = convert_cadnano(
        cadnano_json,
        grid=args.grid,
        sequence=sequence,
        box_side=args.box,
        default_val=args.default_base,
    )

    if args.output is None:
        args.output = args.input.with_suffix('.oxview')

    args.output.write_text(json.dumps(system.to_oxview_dict()))

    result = system.to_oxview_dict()
    strands = result["systems"][0]["strands"]
    n_strands = len(strands)
    n_nucs = sum(len(s["monomers"]) for s in strands)

    print(f"{args.input.name} -> {args.output.name}")
    print(f"  {n_strands} strands, {n_nucs} nucleotides")


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="oxview_py",
        description="Convert DNA origami designs to oxView format.",
    )
    subparsers = parser.add_subparsers(dest="command")

    # -- convert subcommand --
    convert_parser = subparsers.add_parser(
        "convert",
        help="Convert a cadnano JSON file to oxView format.",
    )
    convert_parser.add_argument(
        "input", type=Path,
        help="Input cadnano JSON file",
    )
    convert_parser.add_argument(
        "-o", "--output", type=Path, default=None,
        help="Output oxView file (default: <input>.oxview)",
    )
    convert_parser.add_argument(
        "-g", "--grid", choices=["sq", "he"], default=None,
        help="Lattice type: sq (square) or he (hexagonal). "
             "Auto-detected if omitted.",
    )
    convert_parser.add_argument(
        "-s", "--sequence", default=None,
        help="Scaffold sequence: a nucleotide string or path to a text file. "
             "Random sequence if omitted.",
    )
    convert_parser.add_argument(
        "-b", "--box", type=float, default=None,
        help="Simulation box side length. Auto-calculated if omitted.",
    )
    convert_parser.add_argument(
        "-d", "--default-base", default="N",
        help="Default base for unsequenced staples (IUPAC code, default: N).",
    )
    convert_parser.set_defaults(func=cmd_convert)

    args = parser.parse_args(argv)
    if args.command is None:
        parser.print_help()
        sys.exit(1)

    args.func(args)


if __name__ == "__main__":
    main()
