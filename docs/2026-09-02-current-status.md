# Current project status

Date/time: 2026-09-02 21:21 PDT (America/Los_Angeles)

## Purpose

Record one current, evidence-based snapshot across the simulation toolkit, the
operational rice data workspace, RelocaTE3, and its benchmark runner. Older
dated documents remain useful design and failure records, but their original
status text should not be read as the present state.

## Repositories

- `simulate_data`: branch `feat/pangenome-panel-design`; reusable simulation
  modules and long-read helpers are implemented.
- RelocaTE3: branch `main`; the production pipeline is still centered on the
  short-read workflow and needs a dedicated long-read input/calling path.
- `relocate-benchmark`: branch `main`; it is the maintained benchmark runner.
  It already has uncommitted user work, which this documentation cleanup did
  not modify.
- `make_simulation_new`: operational data workspace, not a Git repository.

No pull-request identifier is recorded in these local checkouts.

## Toolkit status

The CLI exposes TE insertion, SV placement, Illumina, ONT, PacBio, and mixed
TE benchmark-panel modules. PBSIM3 model paths, platform-specific error ratios,
MAF origin parsing, gzip pooling, and PacBio CCS generation are implemented.
ART, PBSIM3, SURVIVOR, samtools, and TEvarSim are declared in the Pixi manifest
and their executables are present in the current environment.

## Operational data status

| Dataset | Status |
| --- | --- |
| Somatic mPing | Complete: 9 mixed samples plus 9 controls |
| Chr1-only somatic mPing | Complete: 9 mixed samples plus 9 controls |
| riceTElib short reads | Complete: 500-event truth, 9 mixed samples plus 9 controls |
| riceTElib divergence | Complete: 18 scenarios, 54 mixed-read tasks, 9 shared controls |
| riceTElib long reads | Incomplete: 15 of 18 tasks complete |
| Multi-genome/pangenome | Phases 1-3 complete, including reads, graph, long-read index, and alignments |
| Standalone `output/mping_only` | Not run |

The incomplete long-read tasks are HiFi 30x replicates 1-3. Their SLURM jobs
reached the 48-hour limit. Each has a partial main FASTQ, an empty control
FASTQ, and no completion sentinel. These files must not be treated as valid
benchmark inputs.

All retained panel reads are gzip-compressed. A scan found no plain `.fastq`
or `.fq` files in the operational result trees.

## Current priority

First finish and validate the three HiFi 30x tasks. Then add a long-read input
seam to RelocaTE3 so users can provide either FASTQ reads or an existing
BAM/CRAM alignment, while the downstream caller consumes one normalized
alignment/evidence representation. Extend `relocate-benchmark` to run and score
that path against the fixed riceTElib truth set.

The implementation plan is
`plans/2026-09-02-relocate3-long-read-support.md`.

The repository intentionally ignores `plans/` as agent working material, so
that plan is present locally but will not appear in ordinary Git status or a
commit unless the repository policy is changed deliberately.

## Commands and validation

The status was checked from completion sentinels, manifest row counts, file
sizes, failed SLURM logs, Git branch/status output, the Pixi manifest, and the
installed executable paths. Repository validation reported:

```text
pixi run format  -> 24 files unchanged
pixi run lint    -> All checks passed
pixi run test    -> 161 passed in 47.56 seconds
git diff --check -> passed
```

No large simulation or caller run was performed on the login node.

## Failures and interpretation

The only data-generation failure in the current priority is the three HiFi
30x SLURM tasks described above. During this documentation check, each Pixi
task printed a successful result, but the Pixi parent process remained open
afterward and was interrupted manually. The formatter, linter, and pytest
results themselves had already completed successfully; no source-test failure
was reported.

The status logic is deliberately sentinel-based: a large or apparently valid
file is not considered complete without the workflow's completion marker and
required companion outputs.

## Next step

Implement only Phase 1 of the linked plan first: preserve the failed HiFi 30x
artifacts, split or resize that work for a valid SLURM partition, rerun the
three tasks, and produce panel-wide realism and integrity QC.
