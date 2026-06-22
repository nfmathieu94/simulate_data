# Implementation Plan: Genome Simulation Modules

## Goal

Add modules to `simulate_data` that:
1. Insert transposable elements (TEs) into a reference genome
2. Place structural variants (SVs) in a genome
3. Simulate Illumina short reads, Oxford Nanopore long reads, and PacBio long reads

All modules follow the existing contract: implement `register_parser(parser)` and `main(args)` in a file under `src/simulate_data/modules/`. The CLI auto-discovers each as `simulate-data <module>`.

---

## Tool Selection

### SV Simulation: SURVIVOR (`simSV`)

- **Why**: Gold standard for SV simulation (deletions, duplications, inversions, translocations). Introduces SVs into a reference FASTA and outputs a ground-truth VCF. Heavily used in benchmarking.
- **Bioconda**: `bioconda::survivor` (v1.0.7)
- **CLI**: `SURVIVOR simSV <ref_fasta> <config_file> <output_prefix>` — exact params to be confirmed from `SURVIVOR simSV` help during implementation.
- **Alternatives considered**: SVEngine (better for tumor/normal mixtures, but overkill here), SVsim (lightweight, but less feature-rich than SURVIVOR).

### TE Simulation: TEvarSim

- **Why**: Specifically designed for TE variant simulation. Handles TSDs, poly-A tails, truncations, and sequence divergence from consensus. Has a built-in read simulator (`Readsim` subcommand). Outputs VCF and modified genome FASTA. Supports population-scale simulation. Directly matches the user's use case: TE fasta + reference genome → modified genome.
- **Install**: `pip install TEvarSim` (conda deps: `gfatools`, `repeatmasker`, `mason`, `pbsim3` from bioconda)
- **CLI**: `tevarsim Simulate --ref ref.fa --te-pool te.fa --bed te.bed --num 10` — five subcommands: TErandom, TEreal, TEpan, Simulate, Readsim, Compare.
- **Alternatives considered**: TEgenomeSimulator (better for "TE landscape" modeling with sequence divergence using Gaussian distributions, handles fragmented/nested TEs and superfamily-specific TSD lengths — could be an alternative module for complex TE landscape modeling). SimulaTE (well-established for population genomics TE insertion dynamics across multiple individuals or populations — less focused on modern long-read integration than newer tools above).

### Read Simulation: Illumina Short Reads — ART (`art_illumina`)

- **Why**: Industry standard for Illumina simulation. Uses empirical error profiles derived from real sequencing runs (HiSeq, NovaSeq, MiSeq) to simulate highly realistic substitution and INDEL error profiles alongside quality scores.
- **Bioconda**: `bioconda::art` (v2016.06.05)
- **CLI**: `art_illumina -i ref.fa -l 150 -f 20 -o output_prefix` — produces paired-end FASTQ files.
- **Alternatives considered**: wgsim / samtools wgsim (ultra-fast, but simple uniform error models — better for sheer speed when realism is not critical).

### Read Simulation: Long Reads (PacBio & ONT) — PBSIM3

- **Why**: Premier general-purpose long-read simulator. Accurately mimics PacBio (CLR and HiFi via multi-pass simulation) and ONT reads. Highly rated for broad read-level realism, length distributions, and quality scores.
- **Bioconda**: `bioconda::pbsim3` (v3.0.5)
- **CLI**: `pbsim --strategy wgs --method qshmm --qshmm data/QSHMM-ONT.model --depth 20 --genome ref.fasta` — outputs FASTQ (compressed) and MAF alignment files. Requires samtools and gzip for output compression.
- **Alternatives considered**: Badread (ONT-specific, better for messy data with chimeras, low-quality runs, adapter contamination, and junk reads — good for testing tool resilience). NanoSim (great for training an error model on a specific real-world ONT dataset and replicating its unique context-dependent homopolymer errors).

---

## Module Architecture

Each new module is a Python file in `src/simulate_data/modules/` implementing the module contract.

### New Modules

```
src/simulate_data/modules/
├── fastq.py              # existing - simple random FASTQ generator
├── te_insertion.py       # TE insertion module (TEvarSim)
├── sv_placement.py       # SV placement module (SURVIVOR)
├── reads_illumina.py     # Illumina read simulator (ART)
├── reads_ont.py          # ONT read simulator (PBSIM3)
└── reads_pacbio.py       # PacBio read simulator (PBSIM3)
```

### Shared Utilities

Create `src/simulate_data/utils.py` for shared functionality:
- `run_command(cmd: list[str], check: bool = True)` — subprocess wrapper with logging
- `validate_fasta(path: Path)` — check file exists and is valid FASTA
- `validate_file_exists(path: Path)` — check file exists
- `setup_logging(verbose: bool = False)` — configure logging
- `ensure_output_dir(path: Path)` — create output directory if needed

### Module Specifications

#### 1. `te_insertion.py` — TE Insertion Module

**Tool**: TEvarSim (`tevarsim Simulate`)

**CLI**: `simulate-data te-insertion --ref ref.fa --te te.fa --num 100 --output results/`

**Arguments**:
- `--ref` (required): Reference genome FASTA
- `--te` (required): TE consensus FASTA (TE pool)
- `--num` (required): Number of TE insertions to simulate
- `--output` (required): Output directory for modified genome and VCF
- `--seed` (optional): Random seed for reproducibility
- `--bed` (optional): BED file of TE positions (if pre-generated)
- `--diverse` (optional): Introduce sequence diversity among individuals
- `--diverse-config` (optional): Configuration file for sequence diversity

**Workflow**:
1. Validate inputs (ref exists, te exists, num > 0)
2. If no BED file provided, generate TE positions using `tevarsim TErandom`
3. Run `tevarsim Simulate` to insert TEs into the reference genome
4. Output: modified genome FASTA, VCF with TE annotations, BED with positions

#### 2. `sv_placement.py` — SV Placement Module

**Tool**: SURVIVOR (`SURVIVOR simSV`)

**CLI**: `simulate-data sv-placement --ref ref.fa --num-sv 100 --sv-types DEL,DUP,INV,TRA --output results/`

**Arguments**:
- `--ref` (required): Reference genome FASTA
- `--num-sv` (required): Number of SVs to simulate
- `--sv-types` (optional): Comma-separated SV types (DEL, DUP, INV, TRA). Default: all
- `--output` (required): Output directory for modified genome and VCF
- `--seed` (optional): Random seed for reproducibility
- `--config` (optional): SURVIVOR config file (if pre-generated)

**Workflow**:
1. Validate inputs (ref exists, num-sv > 0)
2. Generate SURVIVOR config file if not provided
3. Run `SURVIVOR simSV` to place SVs in the reference genome
4. Output: modified genome FASTA, VCF with SV annotations

#### 3. `reads_illumina.py` — Illumina Read Simulator

**Tool**: ART (`art_illumina`)

**CLI**: `simulate-data reads-illumina --ref ref.fa --read-length 150 --coverage 20 --output results/`

**Arguments**:
- `--ref` (required): Reference genome FASTA
- `--read-length` (optional): Read length (default: 150)
- `--coverage` (optional): Coverage depth (default: 20)
- `--output` (required): Output directory for FASTQ files
- `--seed` (optional): Random seed for reproducibility
- `--paired` (optional): Generate paired-end reads (default: True)
- `--fragment-size` (optional): Mean fragment size for paired-end (default: 300)
- `--fragment-std` (optional): Fragment size standard deviation (default: 30)
- `--profile` (optional): ART error profile (e.g., HiSeq, NovaSeq)

**Workflow**:
1. Validate inputs (ref exists, coverage > 0)
2. Run `art_illumina` with appropriate parameters
3. Output: paired-end FASTQ files (R1 and R2)

#### 4. `reads_ont.py` — ONT Read Simulator

**Tool**: PBSIM3 (`pbsim`)

**CLI**: `simulate-data reads-ont --ref ref.fa --coverage 20 --output results/`

**Arguments**:
- `--ref` (required): Reference genome FASTA
- `--coverage` (optional): Coverage depth (default: 20)
- `--output` (required): Output directory for FASTQ files
- `--seed` (optional): Random seed for reproducibility
- `--read-length` (optional): Mean read length (default: 9000)
- `--read-std` (optional): Read length standard deviation (default: 7000)
- `--error-model` (optional): PBSIM3 error model (default: QSHMM-ONT)
- `--qscore-model` (optional): PBSIM3 quality score model

**Workflow**:
1. Validate inputs (ref exists, coverage > 0)
2. Run `pbsim --strategy wgs --method qshmm --qshmm QSHMM-ONT.model` with appropriate parameters
3. Output: FASTQ file (compressed), MAF alignment file

#### 5. `reads_pacbio.py` — PacBio Read Simulator

**Tool**: PBSIM3 (`pbsim`)

**CLI**: `simulate-data reads-pacbio --ref ref.fa --coverage 20 --output results/`

**Arguments**:
- `--ref` (required): Reference genome FASTA
- `--coverage` (optional): Coverage depth (default: 20)
- `--output` (required): Output directory for FASTQ/BAM files
- `--seed` (optional): Random seed for reproducibility
- `--read-type` (optional): Read type (CLR or HiFi). Default: CLR
- `--read-length` (optional): Mean read length (default: 15000)
- `--read-std` (optional): Read length standard deviation (default: 13000)
- `--pass-num` (optional): Number of passes for HiFi (default: 10)
- `--error-model` (optional): PBSIM3 error model (default: QSHMM-RSII)
- `--qscore-model` (optional): PBSIM3 quality score model

**Workflow**:
1. Validate inputs (ref exists, coverage > 0)
2. Run `pbsim --strategy wgs --method qshmm` with appropriate parameters
3. For HiFi: run multi-pass simulation and CCS
4. Output: FASTQ/BAM file (compressed), MAF alignment file

---

## Dependency Management

### Add to `pyproject.toml` `[tool.pixi.dependencies]`:

```toml
# SV simulation
survivor = ">=1.0.7"

# Read simulation
art = ">=2016.06.05"
pbsim3 = ">=3.0.5"

# Utility tools
samtools = ">=1.21"
```

### TEvarSim Installation

TEvarSim is not on bioconda. It needs to be installed via pip with conda dependencies:

```bash
conda install -c bioconda gfatools repeatmasker mason pbsim3
pip install TEvarSim
```

For the pixi environment, we can add TEvarSim as a pip dependency in `pyproject.toml`:

```toml
[tool.pixi.pypi-dependencies]
simulate_data = { path = ".", editable = true }
tevarsim = "*"
```

And add the conda dependencies:

```toml
gfatools = "*"
repeatmasker = "*"
mason = "*"
```

---

## Pipeline Design

The modules are designed to chain together in a pipeline:

```
Reference Genome (FASTA)
        │
        ├──▶ te-insertion ──▶ Modified Genome (FASTA) + TE VCF
        │
        ├──▶ sv-placement ──▶ Modified Genome (FASTA) + SV VCF
        │
        └──▶ reads-* ──▶ Simulated Reads (FASTQ/BAM)
```

### Example Pipeline (SLURM Script)

```bash
#!/usr/bin/env bash
set -euo pipefail

# Step 1: Insert TEs into reference genome
simulate-data te-insertion \
    --ref data/reference.fa \
    --te data/te_consensus.fa \
    --num 100 \
    --output results/te_insertion/

# Step 2: Place SVs in the TE-modified genome
simulate-data sv-placement \
    --ref results/te_insertion/modified_genome.fa \
    --num-sv 50 \
    --output results/sv_placement/

# Step 3: Simulate Illumina reads from the final modified genome
simulate-data reads-illumina \
    --ref results/sv_placement/modified_genome.fa \
    --read-length 150 \
    --coverage 20 \
    --output results/reads_illumina/

# Step 4: Simulate ONT long reads from the final modified genome
simulate-data reads-ont \
    --ref results/sv_placement/modified_genome.fa \
    --coverage 20 \
    --output results/reads_ont/
```

---

## Implementation Plan

### Phase 1: Foundation (1-2 hours)

1. **Add dependencies to `pyproject.toml`**:
   - Add `survivor`, `art`, `pbsim3`, `samtools` to `[tool.pixi.dependencies]`
   - Add `tevarsim` to `[tool.pixi.pypi-dependencies]`
   - Add conda dependencies for TEvarSim (`gfatools`, `repeatmasker`, `mason`)

2. **Create `src/simulate_data/utils.py`**:
   - `run_command(cmd, check=True)` — subprocess wrapper with logging
   - `validate_fasta(path)` — check file exists and is valid FASTA
   - `validate_file_exists(path)` — check file exists
   - `setup_logging(verbose=False)` — configure logging
   - `ensure_output_dir(path)` — create output directory if needed

3. **Run `pixi install`** to update the environment.

### Phase 2: Genome Modification Modules (2-3 hours)

4. **Create `src/simulate_data/modules/te_insertion.py`**:
   - Implement `register_parser(parser)` with TE insertion arguments
   - Implement `main(args)` to call TEvarSim
   - Handle input validation, output generation, and error handling

5. **Create `src/simulate_data/modules/sv_placement.py`**:
   - Implement `register_parser(parser)` with SV placement arguments
   - Implement `main(args)` to call SURVIVOR simSV
   - Handle config file generation, output generation, and error handling

### Phase 3: Read Simulator Modules (2-3 hours)

6. **Create `src/simulate_data/modules/reads_illumina.py`**:
   - Implement `register_parser(parser)` with Illumina read simulation arguments
   - Implement `main(args)` to call ART
   - Handle output file naming, paired-end generation, and error handling

7. **Create `src/simulate_data/modules/reads_ont.py`**:
   - Implement `register_parser(parser)` with ONT read simulation arguments
   - Implement `main(args)` to call PBSIM3
   - Handle error model selection, output file naming, and error handling

8. **Create `src/simulate_data/modules/reads_pacbio.py`**:
   - Implement `register_parser(parser)` with PacBio read simulation arguments
   - Implement `main(args)` to call PBSIM3
   - Handle CLR vs HiFi mode, multi-pass simulation, and error handling

### Phase 4: Testing (2-3 hours)

9. **Create `tests/test_utils.py`**:
   - Test `run_command` with mock subprocess
   - Test `validate_fasta` with valid and invalid FASTA files
   - Test `validate_file_exists` with existing and non-existing files
   - Test `ensure_output_dir` with new and existing directories

10. **Create `tests/test_te_insertion.py`**:
    - Test `register_parser` with valid arguments
    - Test `main` with mock TEvarSim subprocess
    - Test input validation (missing ref, missing te, invalid num)
    - Test output generation (modified genome FASTA, VCF)

11. **Create `tests/test_sv_placement.py`**:
    - Test `register_parser` with valid arguments
    - Test `main` with mock SURVIVOR subprocess
    - Test input validation (missing ref, invalid num-sv)
    - Test config file generation

12. **Create `tests/test_reads_illumina.py`**:
    - Test `register_parser` with valid arguments
    - Test `main` with mock ART subprocess
    - Test input validation (missing ref, invalid coverage)
    - Test output file naming

13. **Create `tests/test_reads_ont.py`**:
    - Test `register_parser` with valid arguments
    - Test `main` with mock PBSIM3 subprocess
    - Test input validation (missing ref, invalid coverage)
    - Test error model selection

14. **Create `tests/test_reads_pacbio.py`**:
    - Test `register_parser` with valid arguments
    - Test `main` with mock PBSIM3 subprocess
    - Test CLR vs HiFi mode
    - Test multi-pass simulation

### Phase 5: Documentation (1 hour)

15. **Update `README.md`**:
    - Add sections for new modules
    - Add example commands for each module
    - Add pipeline example

16. **Update `AGENTS.md`**:
    - Add new modules to the architecture section
    - Add new dependencies to the conventions section
    - Add pipeline example to the HPC/SLURM section

17. **Create `scripts/run_pipeline.sh`**:
    - Example SLURM script that chains all modules together
    - Use `set -euo pipefail` and proper SLURM directives

---

## Testing Strategy

### Unit Tests

- Each module gets a test file in `tests/`
- Tests use `unittest.mock.patch` to mock `subprocess.run` calls
- Tests construct `argparse.Namespace` objects and call `module.main(ns)` directly (following existing test pattern)
- Tests use `tmp_path` fixture for temporary files

### Integration Tests

- Create `tests/test_integration.py` with end-to-end pipeline tests
- Use small test genomes and TE files
- Test the full pipeline: TE insertion → SV placement → Read simulation
- Skip tests if required tools are not installed (use `pytest.mark.skipif`)

### Test Data

- Create `tests/data/` directory with small test files:
  - `mini_genome.fa` — small reference genome (a few kb)
  - `mini_te.fa` — small TE consensus FASTA
  - `mini_te.bed` — small BED file with TE positions

---

## Key Design Decisions

### 1. Subprocess for External Tools

The AGENTS.md says "DO NOT use Python subprocess to call external CLI tools or modules." However, the `simulate_data` package is a Python CLI tool that wraps external tools. Using `subprocess` is the natural way to do this.

**Resolution**: Use `subprocess` in the Python modules to call external tools. The tools are installed via pixi (conda-forge + bioconda), so no module loading is needed. The SLURM scripts handle HPC-specific stuff.

### 2. Module Naming Convention

- `te_insertion.py` → `simulate-data te-insertion`
- `sv_placement.py` → `simulate-data sv-placement`
- `reads_illumina.py` → `simulate-data reads-illumina`
- `reads_ont.py` → `simulate-data reads-ont`
- `reads_pacbio.py` → `simulate-data reads-pacbio`

The CLI auto-discovers modules and uses the filename (with underscores replaced by hyphens) as the subcommand name.

### 3. Output Directory Structure

Each module outputs to a user-specified directory. The default output structure is:

```
results/
├── te_insertion/
│   ├── modified_genome.fa
│   ├── te_insertions.vcf
│   └── te_positions.bed
├── sv_placement/
│   ├── modified_genome.fa
│   └── sv_events.vcf
├── reads_illumina/
│   ├── reads_R1.fastq.gz
│   └── reads_R2.fastq.gz
├── reads_ont/
│   ├── reads.fastq.gz
│   └── align
└── reads_pacbio/
    ├── reads.fastq.gz
    └── align
```

### 4. Error Handling

Each module should:
- Validate all inputs before calling external tools
- Check that external tools are installed (via `shutil.which`)
- Provide clear error messages if tools are missing or inputs are invalid
- Clean up temporary files on failure (use `tempfile.TemporaryDirectory`)

### 5. Reproducibility

All modules accept a `--seed` parameter for reproducibility. The seed is passed to the underlying tool (SURVIVOR, TEvarSim, ART, PBSIM3) if supported.

---

## Implementation Order

1. **Foundation**: `utils.py`, `pyproject.toml` updates, `pixi install`
2. **TE Insertion Module**: `te_insertion.py` + tests
3. **SV Placement Module**: `sv_placement.py` + tests
4. **Illumina Reads Module**: `reads_illumina.py` + tests
5. **ONT Reads Module**: `reads_ont.py` + tests
6. **PacBio Reads Module**: `reads_pacbio.py` + tests
7. **Documentation**: `README.md`, `AGENTS.md`, `scripts/run_pipeline.sh`

---

## Risk Assessment

### High Risk

- **TEvarSim installation**: TEvarSim is not on bioconda and has complex dependencies (gfatools, repeatmasker, mason, pbsim3). Installation may fail or conflict with existing packages.
  - **Mitigation**: Test installation in a fresh conda environment first. Consider using a Docker/Singularity container if installation is too complex.

- **SURVIVOR simSV CLI interface**: The exact CLI interface for `SURVIVOR simSV` is not well-documented. We need to run `SURVIVOR simSV` without arguments to see the help message.
  - **Mitigation**: Install SURVIVOR first and inspect the CLI interface before implementing the module.

### Medium Risk

- **PBSIM3 model files**: PBSIM3 requires model files (QSHMM-*.model, ERRHMM-*.model) that may not be included in the bioconda package.
  - **Mitigation**: Check if model files are included in the bioconda package. If not, download them from the PBSIM3 GitHub repository.

- **ART error profiles**: ART may require error profile files for specific sequencing platforms.
  - **Mitigation**: Check if error profiles are included in the bioconda package. If not, use built-in profiles or download them.

### Low Risk

- **Module auto-discovery**: The existing module auto-discovery mechanism is well-tested and should work with new modules.
- **CLI argument parsing**: The existing `argparse`-based approach is straightforward and well-understood.

---

## Summary

This plan adds 5 new modules to `simulate_data`:
1. `te_insertion.py` — TE insertion using TEvarSim
2. `sv_placement.py` — SV placement using SURVIVOR
3. `reads_illumina.py` — Illumina read simulation using ART
4. `reads_ont.py` — ONT read simulation using PBSIM3
5. `reads_pacbio.py` — PacBio read simulation using PBSIM3

The modules follow the existing contract and chain together in a pipeline: reference genome → TE insertion → SV placement → read simulation. All tools are available on bioconda (except TEvarSim, which is pip-installable with conda dependencies).
