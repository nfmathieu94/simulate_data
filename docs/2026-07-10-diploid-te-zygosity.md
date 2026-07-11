# Diploid TE zygosity simulation

Date: 2026-07-10 America/Los_Angeles

## Purpose

Support RelocaTE-style homozygous and heterozygous insertion benchmarks from
`simulate-data te-insertion`.

## Status

- Added `--af-min` and `--af-max` to `te-insertion`.
- `--num-genomes 2` now works with chromosome subsetting.
- For multi-haplotype output, `final_genome.fa` contains one FASTA record per
  chromosome per haplotype, named `<chrom>_Hap1`, `<chrom>_Hap2`, etc.
- Unmodified chromosomes are duplicated across haplotypes so read simulators
  sample both alleles genome-wide.
- Added `truth_te_zygosity.tsv`, derived from TEvarSim VCF haplotype GT columns.

## Truth Labels

For two haplotypes:

- `1,1` -> `homozygous`
- `1,0` or `0,1` -> `heterozygous`
- `0,0` -> `absent`

The `absent` rows are not truth positives for TE insertion recall but are useful
for confirming that randomized TEvarSim candidate events were not placed in
either haplotype.

## Commands

Focused verification:

```bash
.pixi/envs/default/bin/python -m pytest tests/test_te_insertion.py -v
.pixi/envs/default/bin/ruff check src/simulate_data/modules/te_insertion.py tests/test_te_insertion.py
.pixi/envs/default/bin/simulate-data te-insertion --help
```

## Recommended RelocaTE3 Benchmark Settings

For a mixed homozygous/heterozygous panel:

```bash
simulate-data te-insertion \
  --num-genomes 2 \
  --af-min 0.5 \
  --af-max 0.5
```

This samples each candidate event independently on each haplotype, giving an
expected mixture of homozygous, heterozygous, and absent events.
