# Simulation run preparation to-do

Date: 2026-07-09 America/Los_Angeles

## Current outcome (verified 2026-09-02)

This is a historical preparation checklist, not an active to-do list. The
simulation modules, somatic mPing panels, riceTElib short-read panel, divergence
panel, and multi-genome/pangenome panel were implemented. The riceTElib
long-read panel is 15 of 18 tasks complete; only the three HiFi 30x tasks remain
incomplete after reaching their 48-hour limit.

## Purpose

Prepare reproducible simulated datasets for developing and benchmarking a TE
caller, starting with mPing insertions and later adding one additional TE type,
small simulated genomes, Illumina reads, and long reads.

## Current Recommendation

Use a controlled insertion-only simulation path for the first benchmark series.
Do not rely on RepeatMasker placement as the primary strategy yet.

Rationale:

- The first goal is caller development, so exact truth coordinates and simple
  controlled difficulty tiers are more useful than complex repeat-context
  realism.
- TEvarSim requires a `--known-del` RepeatMasker/UCSC file even for
  insertion-only runs because `TErandom` uses it internally for candidate
  repeat context and placement.
- The downloaded rice repeat GFF3 is useful, but it is not directly compatible
  with the current TEvarSim GFF3 converter because it has `ID`, `Name`, and
  `Note` attributes, not `Class` and `Target`.
- Divergence from the mPing consensus should be controlled by simulation
  parameters such as SNP rate, indel rate, truncation, and optional polyA tails,
  not by the RepeatMasker file.

## Inputs to Confirm

1. Reference FASTA:

   ```text
   data/ref_genome/MSU_r7.fa
   ```

2. mPing consensus FASTA:

   ```text
   data/TE_lib/mping.fa
   ```

   Confirm the FASTA header is TEvarSim-compatible if using TEvarSim:

   ```text
   >mPing#DNA/Harbinger
   ```

3. Optional second TE consensus:

   Decide the second TE family/type later. It should have a clear class/family
   label, for example:

   ```text
   >SecondTE#DNA/Tourist
   ```

4. Downloaded rice repeat annotation:

   ```text
   /rhome/nmath020/wessler_bigdata/rice/ref_data/repeats/osa1_r7.asm.repeat_masked.gff3.gz
   ```

   This file contains mPing-like records in `Note=...`, but no `Class=` or
   `Target=` attributes.

## Output Contract

Genome-modifying workflows should produce:

```text
results/<run_name>/final_genome.fa
results/<run_name>/truth_te.vcf
```

Read simulation workflows should consume `final_genome.fa` and write outputs
under separate read-specific directories:

```text
results/<run_name>/reads_illumina/
results/<run_name>/reads_ont/
results/<run_name>/reads_pacbio/
```

## Preferred Implementation To-Do

### 1. Add a direct insertion-only TE simulator module

Create a new module:

```text
src/simulate_data/modules/te_insert_direct.py
```

CLI subcommand should become:

```bash
pixi run simulate-data te-insert-direct --help
```

Required arguments:

```text
--ref              reference genome FASTA
--te               TE consensus FASTA
--output           output directory
--num              number of insertions
--chroms           chromosome selection, default all
--seed             random seed
```

Important simulation arguments:

```text
--snp-rate
--indel-rate
--indel-ins
--indel-geom-p
--truncated-ratio
--truncated-max-length
--min-distance
--avoid-bed
```

Expected behavior:

- Validate reference FASTA.
- Validate TE FASTA.
- Parse selected chromosomes using existing `parse_chromosome_spec`.
- Sample insertion positions deterministically from `--seed`.
- Avoid overlapping or very close insertions using `--min-distance`.
- Optionally avoid intervals from `--avoid-bed`.
- Mutate each inserted TE copy according to SNP/indel/truncation parameters.
- Insert TE copies into the reference genome.
- Write:

  ```text
  <output>/final_genome.fa
  <output>/truth_te.vcf
  <output>/truth_te.bed
  <output>/truth_te.tsv
  ```

VCF should include enough truth detail for caller evaluation:

```text
CHROM
POS
ID
REF
ALT=<INS:ME:mPing>
INFO fields:
  SVTYPE=INS
  TE=mPing
  END=<pos>
  SVLEN=<inserted_length>
  STRAND=+|-
  CONSENSUS=mPing
  SNP_COUNT=<n>
  INDEL_COUNT=<n>
  TRUNCATED_BASES=<n>
```

### 2. Add focused tests for the direct simulator

Use mini FASTA fixtures in `tests/data/`.

Tests should verify:

- Same seed gives identical `final_genome.fa` and truth files.
- Different seeds change insertion positions or inserted copies.
- `--num` controls the number of truth records.
- Insertions are present in the final FASTA.
- `--chroms Chr1` restricts insertions to Chr1.
- invalid chromosome names fail clearly.
- divergence parameters change inserted sequence content.
- `truth_te.vcf` and `truth_te.bed` have expected coordinates.

### 3. Add a small run script for mPing-only data

Create:

```text
scripts/run_mping_direct.sh
```

SLURM requirements:

```bash
#!/usr/bin/bash -l
set -euo pipefail

#SBATCH --job-name=sim_mping
#SBATCH -p batch
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
```

Path handling:

```bash
BASE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$BASE_DIR"
mkdir -p logs results
```

Example direct run:

```bash
pixi run simulate-data te-insert-direct \
    --ref data/ref_genome/MSU_r7.fa \
    --te data/TE_lib/mping.fa \
    --num 100 \
    --chroms Chr1-5 \
    --seed 42 \
    --snp-rate 0.01 \
    --indel-rate 0.002 \
    --truncated-ratio 0.1 \
    --output results/mping_direct_seed42/
```

### 4. Add read simulation steps

Use the stable final genome path:

```text
results/mping_direct_seed42/final_genome.fa
```

Illumina:

```bash
pixi run simulate-data reads-illumina \
    --ref results/mping_direct_seed42/final_genome.fa \
    --read-length 150 \
    --coverage 20 \
    --seed 456 \
    --output results/mping_direct_seed42/reads_illumina/
```

ONT:

```bash
pixi run simulate-data reads-ont \
    --ref results/mping_direct_seed42/final_genome.fa \
    --coverage 20 \
    --read-length 9000 \
    --read-std 7000 \
    --seed 789 \
    --output results/mping_direct_seed42/reads_ont/
```

PacBio HiFi:

```bash
pixi run simulate-data reads-pacbio \
    --ref results/mping_direct_seed42/final_genome.fa \
    --coverage 20 \
    --read-type HiFi \
    --pass-num 10 \
    --seed 101 \
    --output results/mping_direct_seed42/reads_pacbio/
```

Before relying on these wrappers, confirm the external tools are installed:

```bash
pixi run which art_illumina
pixi run which pbsim
```

If missing, add the required packages to the pixi environment or use cluster
modules in the SLURM scripts.

## TEvarSim Workaround Path

This was the fallback proposed before direct insertion simulation was
implemented. It is retained only as historical troubleshooting context.

### 1. Convert mPing-like GFF3 records into a TEvarSim-compatible `.out`

The downloaded GFF3 has mPing records in `Note=...`, not in
RepeatMasker-style `Class=` and `Target=` attributes. Convert the mPing-like
records into a minimal `.out` file:

```bash
mkdir -p data/repeats

gzip -cd /rhome/nmath020/wessler_bigdata/rice/ref_data/repeats/osa1_r7.asm.repeat_masked.gff3.gz \
  | awk '
BEGIN { OFS=" " }
!/^#/ && NF == 9 && tolower($0) ~ /mping/ {
    chrom = $1
    start = $4
    end = $5
    len = end - start + 1
    id = "mping_" NR
    print 0, 0.0, 0.0, 0.0, chrom, start, end, "(0)", "+", "mPing", "DNA/Harbinger", 1, len, "(0)", id
}
' > data/repeats/mping_known_del.out
```

Check the converted file:

```bash
awk '!/^#/ && NF >= 11 {n++} END {print n}' data/repeats/mping_known_del.out
awk '{print $5}' data/repeats/mping_known_del.out | sort | uniq -c
awk '{print $11}' data/repeats/mping_known_del.out | sort | uniq -c
```

### 2. Run TEvarSim insertion-only

```bash
pixi run simulate-data te-insertion \
    --ref data/ref_genome/MSU_r7.fa \
    --te data/TE_lib/mping.fa \
    --known-del data/repeats/mping_known_del.out \
    --te-type Harbinger \
    --ins-ratio 1.0 \
    --num 100 \
    --chroms Chr1-5 \
    --seed 42 \
    --snp-rate 0.01 \
    --indel-rate 0.002 \
    --truncated-ratio 0.1 \
    --output results/mping_tevarsim_seed42/
```

Expected outputs:

```text
results/mping_tevarsim_seed42/final_genome.fa
results/mping_tevarsim_seed42/truth_te.vcf
```

## Suggested Benchmark Tiers

Start simple and add difficulty gradually.

Tier 1:

```text
mPing only
100 insertions
Chr1 only or Chr1-3
no divergence
Illumina 20x
```

Tier 2:

```text
mPing only
100 insertions
Chr1-5
low divergence: --snp-rate 0.005 --indel-rate 0.001
Illumina 20x and ONT 20x
```

Tier 3:

```text
mPing only
200 insertions
Chr1-12
moderate divergence: --snp-rate 0.01 --indel-rate 0.002
some truncation
Illumina, ONT, and PacBio HiFi
```

Tier 4:

```text
mPing plus one second TE type
separate truth labels per TE family
moderate divergence
mixed short and long reads
```

## Validation Checklist

Before declaring a simulation run usable:

```bash
pixi run simulate-data --help
pixi run simulate-data te-insert-direct --help
pixi run simulate-data reads-illumina --help
pixi run simulate-data reads-ont --help
pixi run simulate-data reads-pacbio --help
```

Check final genome:

```bash
test -s results/<run_name>/final_genome.fa
grep -c "^>" results/<run_name>/final_genome.fa
```

Check truth records:

```bash
test -s results/<run_name>/truth_te.vcf
grep -vc "^#" results/<run_name>/truth_te.vcf
```

Check reads:

```bash
find results/<run_name> -maxdepth 2 -type f | sort
```

Run code checks before committing wrapper changes:

```bash
pixi run format
pixi run lint
pixi run test
```

If `pixi run format` hangs on the HPC filesystem, record the failure and run:

```bash
.pixi/envs/default/bin/ruff check src/ tests/
.pixi/envs/default/bin/pytest tests/ -v
git diff --check
```

## Open Decisions

- Choose the second TE family after mPing.
- Decide whether to simulate on full chromosomes or extract smaller test
  chromosomes/windows for faster early runs.
- Decide whether insertion positions should be uniform random, avoid genes,
  avoid repeats, or deliberately include repeat-rich hard cases.
- Decide whether direct simulator should write haploid-only genomes initially
  or support diploid/haplotype-aware truth sets later.
