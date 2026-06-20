# Simulating Data

A modular library for simulating genomic data, designed to support the testing and benchmarking of bioinformatic tools.

## Overview

Each simulation type is implemented as an independent module that can be invoked from the command line. The project is built as a [pixi](https://pixi.sh)-managed Python package, providing reproducible conda environments and a consistent CLI interface.

New simulation modules are auto-discovered: drop a Python file in `src/simulate_data/modules/`, implement `register_parser()` and `main()`, and it becomes a new `simulate-data <module>` subcommand.

## Project Structure

```
simulate_data/
├── pyproject.toml          # Package metadata + pixi project config
├── src/simulate_data/      # Python package source
│   ├── cli.py              # CLI entry point (auto-discovers modules)
│   └── modules/            # Individual simulation modules (e.g. fastq.py)
├── scripts/                # SLURM submission scripts
├── tests/                  # pytest test suite
├── data/                   # Input data (read-only)
├── results/                # Simulation outputs
└── logs/                   # Slurm log files
```

## Quick Start

```bash
# Install the pixi environment
pixi install

# List available simulation modules
simulate-data --help

# Generate 1000 simulated FASTQ reads of length 150
simulate-data fastq -n 1000 -l 150 --seed 42 -o results/simulated_reads.fastq
```

## Running on the Cluster

Submit a SLURM job using the example script:

```bash
sbatch scripts/run.sh
```

## Adding a New Module

1. Create a new file in `src/simulate_data/modules/` (e.g., `bam.py`).
2. Implement two functions:
   - `register_parser(subparser)` — define the module's CLI arguments.
   - `main(args)` — execute the simulation.
3. The module is automatically available as a subcommand: `simulate-data bam`.

## Testing

```bash
pixi run test
```

## Code Quality

```bash
pixi run lint     # ruff check
pixi run format   # black
```
