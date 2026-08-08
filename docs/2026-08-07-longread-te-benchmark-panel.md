# Design: Long-read TE benchmark panel (riceTElib-LR)

Date/time: 2026-08-07 America/Los_Angeles

## Purpose

Generate simulated **long-read** (ONT high-accuracy and PacBio HiFi) benchmark
data over the **existing riceTElib truth set**, so that adding long-read support
to RelocaTE3 can be developed and scored against the same 500 insertion events
already used for the short-read benchmark.

The panel is a datagen product. It lives on the `make_simulation_new` side of the
datagen / benchmark boundary; `relocate-benchmark` consumes it later by path.

## Decisions

| Decision | Choice | Rationale |
|---|---|---|
| Platforms | ONT-HQ (~99%) + PacBio HiFi (~99.9%) | Both are what people actually run for TE detection today. Read accuracy is the binding constraint on TSD-exact calling. |
| Truth set | Reuse `results/riceTElib_benchmark/` catalog, read-only | Makes short-read vs long-read a controlled comparison where read technology is the only variable. |
| Coverage | 5x, 15x, 30x; 3 replicates | Point-for-point comparability with the short-read panel. |
| Mixture | Unchanged: 4 components x 2 haplotypes | Component genomes already exist. Long reads should be markedly better at low-VAF somatic detection; that is worth measuring. |
| Genome scope | Whole genome (12 chroms), insertions Chr1-only | Inherited from the short-read panel. Preserves realistic mapping ambiguity at the cost of ~92% background reads. |
| Engine location | `simulate_data` toolkit | Mirrors how `build_multite_panel.py` already imports `base_panel._simulate_reads`. Avoids forking read-simulation logic into a non-git-tracked directory. |
| Compression | Pool directly into gzip; never write uncompressed FASTQ | pbsim3 emits `.fq.gz` and ccs emits `.fastq.gz` natively, so the short-read write-then-compress step is unnecessary. |

## Verified facts (2026-08-07)

Established by direct execution, not assumption:

- `pbsim3` is already pinned in the `simulate_data` pixi env.
- **`reads_ont.py` and `reads_pacbio.py` are broken and have never been run.**
  They pass a bare model name; PBSIM3 requires a path:
  `ERROR: Cannot open file: QSHMM-ONT`. Works with
  `--qshmm $CONDA_PREFIX/data/QSHMM-ONT.model`.
- Seven models ship in `$CONDA_PREFIX/data/`: `QSHMM-{ONT,ONT-HQ,RSII}.model`,
  `ERRHMM-{ONT,ONT-HQ,RSII,SEQUEL}.model`.
- PBSIM3 output contract, **per contig**: `<prefix>_0001.fq.gz`,
  `<prefix>_0001.maf.gz`, `<prefix>_0001.ref`. Single-end. **No SAM** — so the
  existing `_count_origin_support()` pysam path does not apply.
- `--difference-ratio` default `6:55:39` is PacBio RS II. PBSIM3 documents
  `39:24:36` for ONT and `22:45:33` for Sequel.
- HiFi requires `ccs`, absent from the env. `pbccs 6.4.0` exists on bioconda
  linux-64.
- **The full HiFi chain was tested end-to-end and works**: pbsim3
  `--pass-num 10` subread BAM (correct `@RG PL:PACBIO READTYPE=SUBREAD`,
  `movie/zmw/idx` naming) -> `ccs` -> HiFi FASTQ. 37 reads from 40 ZMWs,
  mean length 14.5 kb.
- **ccs yield is ~92.5%**, so requested pbsim depth overstates delivered HiFi
  depth. Actual yield must be recorded per component.
- Existing catalog is complete: 500 events, 8 component genomes + reference,
  ~378 MB each. Short-read panel totals 65 GB.

## Architecture

Two layers, matching the existing split.

### Layer 1 — `simulate_data` (toolkit, git-tracked)

- `src/simulate_data/longread.py` (new) — model-name resolution against
  `$CONDA_PREFIX/data/`, pbsim3 invocation, per-contig output collection, MAF
  parsing, ccs wrapper.
- `src/simulate_data/modules/reads_ont.py` — fix model path; add
  `--difference-ratio`, `--accuracy-mean`, `--length-min/max`, errhmm method.
- `src/simulate_data/modules/reads_pacbio.py` — same fix; real HiFi via ccs.
- `src/simulate_data/modules/te_benchmark_panel.py` — add
  `_simulate_long_reads(settings, output, coverage, replicate, platform)`
  mirroring `_simulate_reads`: same 8-component loop, same `_seed()` derivation,
  same pooling and manifests; single-end, MAF-based QC, gzip output.
- `pyproject.toml` — add `pbccs ==6.4.0` to pixi dependencies.

### Layer 2 — `make_simulation_new` (panel driver)

- `config/riceTElib_lr_benchmark.toml` — catalog path, one table per platform.
- `pipeline/make_riceTElib_lr_benchmark/`
  - `build_longread_panel.py` — validates and reuses the existing catalog; runs
    the read stage. Never rebuilds the catalog.
  - `01_link_catalog.sh` — validate `.catalog.complete`, link genomes read-only.
  - `02_simulate_reads.sh` — SLURM array, 18 tasks.
  - `submit.sh`, `test_build_longread_panel.py`, `README.md`
- `docs/2026-08-07-longread-benchmark-workflow.md`

## Data flow

```
results/riceTElib_benchmark/          (existing, READ-ONLY)
  genomes/{baseline,clone40,clone20,clone10}.hap{1,2}.fa
  genomes/reference.fa
  truth_events.tsv, component_coordinates.tsv
        |
        |  per platform x coverage x replicate -> 8 component runs
        |  at depth = coverage * weight / 2
        |
        +-- ONT-HQ : --method qshmm  --qshmm QSHMM-ONT-HQ.model
        |            --difference-ratio 39:24:36, 1 pass      -> _NNNN.fq.gz
        |
        +-- HiFi   : --method errhmm --errhmm ERRHMM-SEQUEL.model
                     --difference-ratio 22:45:33 --pass-num 10 -> _NNNN.bam
                     -> ccs                                     -> .fastq.gz
        |
        |  pool into a single gzip stream, read names prefixed
        |  <sample>:<component>:h<hap>:s<seed>
        v
results/riceTElib_lr_benchmark/reads/<platform>/cov{5,15,30}x_rep{1,2,3}/
  reads.fastq.gz
  component_manifest.tsv     (incl. actual ccs yield)
  observed_support.tsv, panel_summary.tsv
  .complete
results/riceTElib_lr_benchmark/reads/<platform>/cov..._reference_control/
  reads.fastq.gz, .complete
```

18 array tasks = 2 platforms x 3 coverages x 3 replicates.

## Origin-support QC

The short-read panel counts ART-origin junction support from a SAM. PBSIM3 emits
MAF instead, so `longread.py` provides a MAF reader yielding
`(contig, ref_start, ref_length)` per read. The counting logic is otherwise
unchanged: the same bisect over `component_coordinates.tsv`, the same
present/absent tally, the same `observed_support.tsv` and `panel_summary.tsv`
schema. Reusing the schema keeps long-read VAF QC directly comparable to
short-read.

## Compression contract

Explicit requirement. No uncompressed FASTQ is ever written:

- pbsim3 emits `.fq.gz`; ccs emits `.fastq.gz`. Both are read through `gzip`
  and appended into an output `gzip` stream.
- Pooled output is `reads.fastq.gz` at the same `GZIP_LEVEL` used by
  `pipeline/fastq_compression.py`, for consistency with the short-read panel.
- Post-stage validation refuses any residual `*.fastq` / `*.fq`, matching
  `02_simulate_reads.sh` in the short-read pipeline.
- `run_metadata.json` records compression and level.

## Error handling

Follows the conventions already in these pipelines:

- Validate the borrowed catalog's `.catalog.complete` before any work; the
  long-read panel never writes into the catalog directory.
- Sentinel-gated resume: skip on `.complete`; refuse a non-empty directory that
  lacks its sentinel rather than overwriting.
- Fail on empty or unreadable gzip output, and on residual uncompressed FASTQ.
- HiFi: fail if ccs yield falls below a configurable floor, since silent
  under-coverage would otherwise look like a caller recall problem.

## Testing

Mirrors `test_build_multite_panel.py`:

- Unit: model-name resolution, pbsim command construction per platform, MAF
  parsing, mixture/seed math, gzip pooling round-trip.
- End-to-end: tiny fixture genome, both platforms, asserting read counts,
  compression, sentinels, and manifest schema.

## Risks

1. **ccs yield ~92.5%** — delivered HiFi depth is below requested. Recorded per
   component; a yield floor fails loudly rather than silently under-covering.
2. **Storage** — ~65 GB per technology, ~130 GB total. `/bigdata` is at 86%
   (3.7 TB free). Adequate, but not free.
3. **HiFi 30x is the expensive corner** — 10 passes over a 380 Mb genome x 8
   components, plus CPU-heavy polishing. Gets its own thread/time allocation
   rather than uniform array resources.
4. **Whole-genome scope** means ~92% background reads. Restricting to Chr1 would
   cut cost ~12x at the price of realistic mapping ambiguity. Deliberately not
   taken for the headline panel.

## Next steps

Implementation plan follows this document.
