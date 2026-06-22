# AGENTS.md

Repository-specific guidance for agents working in `simulate_data`.
The global `~/.config/opencode/AGENTS.md` (HPC bioinformatician defaults) applies on top of this.

## Commands

This is a pixi-managed project. All developer commands run through pixi:

```bash
pixi run test      # pytest tests/ -v
pixi run lint      # ruff check src/ tests/
pixi run format    # black src/ tests/
```

Run a single test file or case:

```bash
pixi run pytest tests/test_cli.py -v
pixi run pytest tests/test_cli.py::test_fastq_reproducible_with_seed -v
```

The CLI is installed as the `simulate-data` console script (`simulate_data.cli:main`). After environment changes, `pixi install` refreshes the editable install.

## Architecture

- **src layout**: package source lives in `src/simulate_data/`. Hatchling builds `src/simulate_data` as the wheel package.
- **CLI entry**: `simulate-data` → `simulate_data.cli:main`. `cli.py` auto-discovers modules and registers each as a subcommand.
- **Module contract** (the one hard rule for adding modules): drop a file in `src/simulate_data/modules/`, implement `register_parser(parser)` and `main(args)`, and it becomes `simulate-data <module>`. The module's `__doc__` first line becomes the subcommand's `--help` text. Modules prefixed with `_` are skipped. A module missing `register_parser` or `main` is silently ignored by discovery. Subcommand names use hyphens (e.g., `te_insertion.py` → `simulate-data te-insertion`).
- **Test layout**: pytest suite in `tests/`. Tests construct `argparse.Namespace` objects and call `module.main(ns)` directly rather than shelling out to the CLI. Namespace attributes must use argparse destination names (e.g., `num_reads`, not `n`).
- **Simulation modules**: `te_insertion.py` (TEvarSim), `sv_placement.py` (SURVIVOR), `reads_illumina.py` (ART), `reads_ont.py` (PBSIM3), `reads_pacbio.py` (PBSIM3). Each wraps an external tool via `subprocess` (see `utils.run_command`).
- **Chromosome selection**: TE insertion and SV placement modules accept `--chroms` with specs like `Chr1`, `Chr2-5`, `Chr1,Chr3`, or `all`. Parsed by `utils.parse_chromosome_spec`.
- **Shared utilities**: `src/simulate_data/utils.py` provides FASTA validation, chromosome extraction/merging, subprocess wrappers, and tool installation checks.

## Conventions

- **Formatting**: `black`, line-length 88. **Linting**: `ruff`, line-length 88, rules `E, F, I, W` (see `pyproject.toml`). Run `pixi run format` then `pixi run lint` before committing.
- **Python**: `>=3.10,<3.13`. Runtime deps: `pandas`, `biopython` (`biopython>=1.83`), `pysam`.
- **HPC / SLURM**: cluster jobs are submitted from the repo root via `sbatch scripts/run.sh`. Scripts use `set -euo pipefail` and `cd` to `PROJECT_ROOT` before running `simulate-data`. Do not run heavy generation on the login node — wrap it in a SLURM script under `scripts/`.
- **Directory roles**: `data/` is read-only input; `results/` holds simulation outputs; `logs/` holds Slurm logs (`%x_%j.err`). Contents of `results/`, `logs/`, and `data/` are gitignored (only `.gitkeep` is tracked).

## Verification

Before declaring work done, run from the repo root:

```bash
pixi run format
pixi run lint
pixi run test
```

`ruff` and `black` configs live in `pyproject.toml`; there is no separate `ruff.toml` or `.black` file.
