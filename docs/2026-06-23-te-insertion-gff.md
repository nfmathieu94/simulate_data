# TE insertion GFF compatibility

Date: 2026-06-23

## Purpose

Support the mPing simulation workflow using a RepeatMasker GFF3 annotation as
the TEvarSim known deletion source.

## Status

- Added RepeatMasker GFF3 conversion inside `te-insertion`.
- Added `--te-type` forwarding to `tevarsim TErandom`.
- Fixed per-chromosome TEvarSim output merging into `final_genome.fa`.
- Updated the external mPing simulation script to use:
  - `input/ref_genome/MSU_r7.fa.RepeatMasker.out.gff`
  - `input/TE_lib/mping_superfam_header.fa`
  - `--te-type Harbinger`
  - `--ins-ratio 1.0`

## Validation

- `ruff check src/ tests/`: passed.
- `pytest tests/ -v`: 107 passed.
- Smoke-tested one real Chr10 mPing insertion into `/tmp`; outputs included
  `Sim_Chr10.fa`, `Sim_Chr10.vcf`, and `final_genome.fa`.

## Notes

The requested chromosomes had fewer than 200 qualifying Harbinger deletion
records in the GFF, so the workflow uses insertion-only simulation
(`--ins-ratio 1.0`) for the mPing genome.
