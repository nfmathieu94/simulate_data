"""Command-line interface with auto-discovered simulation modules."""

import argparse
import importlib
import pkgutil
import sys

from simulate_data import modules


def _discover_modules():
    """Return sorted list of module names found in the modules package."""
    return sorted(
        name
        for _, name, _ in pkgutil.iter_modules(modules.__path__)
        if not name.startswith("_")
    )


def build_parser():
    """Construct the top-level argument parser with subcommands per module."""
    parser = argparse.ArgumentParser(
        prog="simulate-data",
        description="Simulate various types of genomic data.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    for name in _discover_modules():
        mod = importlib.import_module(f"simulate_data.modules.{name}")
        if not hasattr(mod, "register_parser") or not hasattr(mod, "main"):
            continue
        # Use hyphenated name as the subcommand (e.g., te_insertion -> te-insertion)
        subcommand_name = name.replace("_", "-")
        sub = subparsers.add_parser(
            subcommand_name,
            help=mod.__doc__.strip().splitlines()[0] if mod.__doc__ else name,
            description=mod.__doc__,
        )
        mod.register_parser(sub)
        sub.set_defaults(_main=mod.main)

    return parser


def main(argv=None):
    """Entry point for the ``simulate-data`` console script."""
    parser = build_parser()
    args = parser.parse_args(argv)
    return args._main(args)


if __name__ == "__main__":
    sys.exit(main())
