# Simulating Data

A modular library for simulating genomic data for testing and benchmarking bioinformatic tools.

## Table of Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Available Subcommands](#available-subcommands)
- [Usage Examples](#usage-examples)
  - [TE Insertion Only](#te-insertion-only)
  - [SV Placement Only](#sv-placement-only)
  - [TE Insertion + SV Placement (Combined Pipeline)](#te-insertion--sv-placement-combined-pipeline)
  - [Read Simulation from a Modified Genome](#read-simulation-from-a-modified-genome)
- [Parameter Reference](#parameter-reference)
  - [fastq](#fastq)
  - [te-insertion](#te-insertion)
  - [sv-placement](#sv-placement)
  - [reads-illumina](#reads-illumina)
  - [reads-ont](#reads-ont)
  - [reads-pacbio](#reads-pacbio)
- [Chromosome Selection](#chromosome-selection)
- [Running on HPC (SLURM)](#running-on-hpc-slurm)
- [Adding a New Module](#adding-a-new-module)
- [Testing & Code Quality](#testing--code-quality)
- [License](#license)

## Introduction

`simulate_data` provides a unified command-line interface for generating simulated genomic datasets. Each simulation type is an independent module auto-discovered by the CLI. Modules wrap established external tools (TEvarSim, SURVIVOR, ART, PBSIM3) and handle input validation, chromosome selection, and output management.

The project is built as a [pixi](https://pixi.sh)-managed Python package, providing reproducible conda environments and a consistent CLI interface.

## Installation

```bash
# Clone the repository
git clone <repo-url>
cd simulate_data

# Install the pixi environment (includes all dependencies)
pixi install
```

Pixi installs the `simulate-data` console script into `.pixi/envs/default/bin/`, which is **not** on your PATH by default. You have two options:

```bash
# Option A: prefix every command with pixi run (recommended)
pixi run simulate-data --help

# Option B: activate the environment once, then use bare commands
pixi shell
simulate-data --help
```

All examples below use the `pixi run` prefix. If you've activated the shell with `pixi shell`, you can drop the `pixi run` prefix.

## Available Subcommands

| Subcommand | Description | Underlying Tool |
|---|---|---|
| `fastq` | Generate random FASTQ reads (no reference needed) | Built-in |
| `te-insertion` | Insert transposable elements into a reference genome | TEvarSim |
| `sv-placement` | Place structural variants in a reference genome | SURVIVOR |
| `reads-illumina` | Simulate Illumina short reads from a reference | ART |
| `reads-ont` | Simulate Oxford Nanopore long reads from a reference | PBSIM3 |
| `reads-pacbio` | Simulate PacBio CLR or HiFi reads from a reference | PBSIM3 |

## Usage Examples

All examples assume you are in the repository root and have run `pixi install`.

### TE Insertion Only

Insert 100 copies of the mPing transposon into chromosomes 1–5 of the MSU7 rice reference genome:

```bash
pixi run simulate-data te-insertion \
    --ref data/ref_genome/MSU_r7.fa \
    --te data/TE_lib/mping.fa \
    --known-del data/ref_genome/MSU_r7.fa.RepeatMasker.out \
    --num 100 \
    --chroms Chr1-5 \
    --snp-rate 0.01 \
    --indel-rate 0.002 \
    --sense-strand-ratio 0.7 \
    --seed 42 \
    --output results/te_insertion/
```

**Output:** `results/te_insertion/modified_genome.fa` (TE-integrated genome) and a VCF file with insertion annotations.

### SV Placement Only

Place 50 structural variants (deletions, duplications, inversions, translocations) across the entire reference:

```bash
pixi run simulate-data sv-placement \
    --ref data/ref_genome/MSU_r7.fa \
    --num-sv 50 \
    --sv-types DEL,DUP,INV,TRA \
    --seed 123 \
    --output results/sv_placement/
```

**Output:** `results/sv_placement/modified_genome.fa` (SV-integrated genome) and a VCF file with variant annotations.

### TE Insertion + SV Placement (Combined Pipeline)

First insert TEs, then place SVs in the TE-modified genome:

```bash
# Step 1: Insert 100 mPing TEs into Chr1–Chr5
pixi run simulate-data te-insertion \
    --ref data/ref_genome/MSU_r7.fa \
    --te data/TE_lib/mping.fa \
    --known-del data/ref_genome/MSU_r7.fa.RepeatMasker.out \
    --num 100 \
    --chroms Chr1-5 \
    --seed 42 \
    --output results/te_insertion/

# Step 2: Place 50 SVs in the TE-modified genome
pixi run simulate-data sv-placement \
    --ref results/te_insertion/modified_genome.fa \
    --num-sv 50 \
    --sv-types DEL,DUP,INV,TRA \
    --seed 123 \
    --output results/sv_placement/
```

**Output:** `results/sv_placement/modified_genome.fa` is the final genome with both TEs and SVs, along with two VCF files (one per step) serving as ground truth.

### Read Simulation from a Modified Genome

After creating a modified genome (via TE insertion, SV placement, or both), simulate reads from it:

```bash
# Illumina paired-end reads (150 bp, 20x coverage)
pixi run simulate-data reads-illumina \
    --ref results/sv_placement/modified_genome.fa \
    --read-length 150 \
    --coverage 20 \
    --seed 456 \
    --output results/reads_illumina/

# Oxford Nanopore long reads (20x coverage)
pixi run simulate-data reads-ont \
    --ref results/sv_placement/modified_genome.fa \
    --coverage 20 \
    --seed 789 \
    --output results/reads_ont/

# PacBio HiFi reads (20x coverage, 10 passes)
pixi run simulate-data reads-pacbio \
    --ref results/sv_placement/modified_genome.fa \
    --coverage 20 \
    --read-type HiFi \
    --pass-num 10 \
    --seed 101 \
    --output results/reads_pacbio/
```

## Parameter Reference

### fastq

Generate random FASTQ reads without a reference genome. Useful for quick testing of FASTQ parsers and pipelines.

| Flag | Type | Default | Description |
|---|---|---|---|
| `-n, --num-reads` | int | 10 | Number of reads to generate |
| `-l, --read-length` | int | 150 | Length of each read (bp) |
| `-o, --output` | file | stdout | Output FASTQ file |
| `--seed` | int | None | Random seed for reproducibility |

### te-insertion

Insert transposable elements (TEs) from a consensus FASTA into specified chromosomes of a reference genome using TEvarSim.

| Flag | Type | Default | Description |
|---|---|---|---|
| `--ref` | path | **required** | Reference genome FASTA |
| `--te` | path | **required** | TE consensus FASTA (TE pool) |
| `--known-del` | path | **required** | Known TE deletion annotation: TEvarSim-compatible RepeatMasker `.out`/UCSC `.txt`, or RepeatMasker GFF/GFF3 |
| `--num` | int | **required** | Number of TE events to simulate per selected chromosome |
| `--output` | path | **required** | Output directory for modified genome and VCF |
| `--chroms` | string | `all` | Chromosomes to insert TEs into (see [Chromosome Selection](#chromosome-selection)) |
| `--seed` | int | None | Random seed for reproducibility |
| `--bed` | path | None | BED file of pre-generated TE positions |
| `--num-genomes` | int | 1 | Number of genomes/haplotypes for TEvarSim `Simulate` |
| `--ins-ratio` | float | 0.6 | Proportion of TE events that are insertions |
| `--te-type` | string | None | TE family/superfamily filter; can be repeated or comma-separated |
| `--snp-rate` | float | 0.02 | SNP mutation rate per base when generating inserted TE copies |
| `--indel-rate` | float | 0.005 | Indel mutation rate per base when generating inserted TE copies |
| `--indel-ins` | float | 0.4 | Proportion of TE-copy indels that are insertions |
| `--indel-geom-p` | float | 0.7 | Geometric distribution parameter for TE-copy indel lengths |
| `--truncated-ratio` | float | 0.3 | Proportion of inserted TE copies to truncate |
| `--truncated-max-length` | float | 0.5 | Maximum proportion of each inserted TE copy that can be truncated |
| `--polyA-ratio` | float | 0.8 | Proportion of inserted TE copies to add a polyA tail |
| `--polyA-min` | int | 5 | Minimum polyA tail length |
| `--polyA-max` | int | 20 | Maximum polyA tail length |
| `--sense-strand-ratio` | float | 0.5 | Proportion of TE insertions simulated on the sense strand |

### sv-placement

Place structural variants (SVs) in a reference genome using SURVIVOR. Supports deletions (DEL), duplications (DUP), inversions (INV), and translocations (TRA).

| Flag | Type | Default | Description |
|---|---|---|---|
| `--ref` | path | **required** | Reference genome FASTA |
| `--num-sv` | int | **required** | Number of SVs to simulate |
| `--output` | path | **required** | Output directory for modified genome and VCF |
| `--sv-types` | string | `DEL,DUP,INV,TRA` | Comma-separated SV types |
| `--chroms` | string | `all` | Chromosomes to place SVs in (see [Chromosome Selection](#chromosome-selection)) |
| `--seed` | int | None | Random seed for reproducibility |
| `--config` | path | None | Pre-generated SURVIVOR config file |

### reads-illumina

Simulate Illumina paired-end or single-end reads from a reference genome using ART with empirical error profiles.

| Flag | Type | Default | Description |
|---|---|---|---|
| `--ref` | path | **required** | Reference genome FASTA |
| `--output` | path | **required** | Output directory for FASTQ files |
| `--read-length` | int | 150 | Read length in bp |
| `--coverage` | float | 20.0 | Coverage depth (fold) |
| `--paired` | flag | True | Generate paired-end reads |
| `--fragment-size` | int | 300 | Mean fragment size for paired-end (bp) |
| `--fragment-std` | int | 30 | Fragment size standard deviation (bp) |
| `--profile` | path | None | ART error profile (e.g., HiSeq, NovaSeq) |
| `--seed` | int | None | Random seed for reproducibility |

### reads-ont

Simulate Oxford Nanopore (ONT) long reads from a reference genome using PBSIM3 with quality score hidden Markov models (QSHMM).

| Flag | Type | Default | Description |
|---|---|---|---|
| `--ref` | path | **required** | Reference genome FASTA |
| `--output` | path | **required** | Output directory for FASTQ files |
| `--coverage` | float | 20.0 | Coverage depth (fold) |
| `--read-length` | int | 9000 | Mean read length in bp |
| `--read-std` | int | 7000 | Read length standard deviation in bp |
| `--error-model` | string | `QSHMM-ONT` | PBSIM3 error model |
| `--qscore-model` | string | None | PBSIM3 quality score model |
| `--seed` | int | None | Random seed for reproducibility |

### reads-pacbio

Simulate PacBio continuous long reads (CLR) or HiFi circular consensus sequencing (CCS) reads from a reference genome using PBSIM3.

| Flag | Type | Default | Description |
|---|---|---|---|
| `--ref` | path | **required** | Reference genome FASTA |
| `--output` | path | **required** | Output directory for FASTQ/BAM files |
| `--coverage` | float | 20.0 | Coverage depth (fold) |
| `--read-type` | enum | `CLR` | Read type: `CLR` or `HiFi` |
| `--read-length` | int | 15000 | Mean read length in bp |
| `--read-std` | int | 13000 | Read length standard deviation in bp |
| `--pass-num` | int | 10 | Number of passes for HiFi sequencing |
| `--error-model` | string | `QSHMM-RSII` | PBSIM3 error model |
| `--qscore-model` | string | None | PBSIM3 quality score model |
| `--seed` | int | None | Random seed for reproducibility |

## Chromosome Selection

The `te-insertion` and `sv-placement` modules accept a `--chroms` argument that controls which chromosomes of the reference genome are modified. Unmodified chromosomes are preserved in the final output.

| Format | Example | Description |
|---|---|---|
| Single chromosome | `Chr1` | Insert/place in only Chr1 |
| Comma-separated list | `Chr1,Chr3,Chr5` | Insert/place in the listed chromosomes |
| Range (numeric suffix) | `Chr2-5` | Expands to Chr2, Chr3, Chr4, Chr5 |
| Range (full names) | `Chr2-Chr5` | Same as above |
| All chromosomes | `all` | Use the entire reference genome (default) |

You can combine ranges and single chromosomes: `Chr1-3,Chr8` → Chr1, Chr2, Chr3, Chr8.

## Running on HPC (SLURM)

Submit a SLURM job using one of the example scripts:

```bash
sbatch scripts/run.sh          # Simple fastq generation
sbatch scripts/run_pipeline.sh # Full TE → SV → reads pipeline
```

All SLURM scripts use `set -euo pipefail` and `cd` to the project root before running `simulate-data`. Do not run heavy generation on the login node — wrap it in a SLURM script under `scripts/`.

## Adding a New Module

1. Create a new file in `src/simulate_data/modules/` (e.g., `bam.py`).
2. Implement two functions:
   - `register_parser(parser)` — define the module's CLI arguments.
   - `main(args)` — execute the simulation.
3. The module is automatically available as a subcommand: `simulate-data bam`.

Modules prefixed with `_` (e.g., `_helpers.py`) are skipped by auto-discovery. A module missing `register_parser` or `main` is silently ignored.

## Testing & Code Quality

```bash
pixi run test      # Run pytest suite
pixi run lint      # Run ruff linter
pixi run format    # Run black formatter
```

Run a single test file or case:

```bash
pixi run pytest tests/test_cli.py -v
pixi run pytest tests/test_te_insertion.py::TestMain::test_main_all_chroms -v
```

## License

MIT
