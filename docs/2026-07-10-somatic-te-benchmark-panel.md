# Somatic TE benchmark panel

Date: 2026-07-10 America/Los_Angeles

## Purpose and status

The `te-benchmark-panel` module creates mixed germline and subclonal somatic TE
datasets for RelocaTE development. Somatic insertions are heterozygous in a
fraction of cells, so insertion VAF is half the configured cellular fraction.

The module uses separate `catalog` and `reads` stages. This permits one catalog
to be shared safely by SLURM coverage/replicate array jobs.

## Default biological model

| Class | Cellular fraction | Expected insertion VAF |
|---|---:|---:|
| homozygous | 100% | 100% |
| heterozygous | 100% | 50% |
| somatic | 10% | 5% |
| somatic | 20% | 10% |
| somatic | 40% | 20% |

The nested mixture is 60% baseline, 20% clone40, 10% clone20, and 10% clone10.

## Commands

```bash
pixi run simulate-data te-benchmark-panel \
  --config config/somatic_mping_panel.toml \
  --output results/somatic_mping_panel \
  --stage catalog

pixi run simulate-data te-benchmark-panel \
  --config config/somatic_mping_panel.toml \
  --output results/somatic_mping_panel \
  --stage reads \
  --coverage 15 \
  --replicate 1
```

The read stage writes pooled FASTQ files, a component manifest, matched
reference-only FASTQs, and `observed_support.tsv`. That QC table is based on
ART read origins and reports the insertion-junction/reference-spanner balance
used by RelocaTE's characterizer.

Coordinates in `truth_events.tsv` and VCF are the 1-based reference anchor base
immediately before the insertion. `component_coordinates.tsv` uses 0-based
coordinates in each generated component FASTA because it describes ART origins.

## Idempotence

Completed catalog and read stages have `.catalog.complete` and `.complete`
sentinels. Completed work is skipped. Non-empty incomplete output directories
cause a failure rather than silent overwrite.
