# Multi-genome / pangenome TE benchmark panel

Date/time: 2026-08-08 America/Los_Angeles

## Purpose

Ground truth for two upcoming RelocaTE3 capabilities:

1. **Multi-FASTA support.** A user supplies several reference genomes plus
   sample reads; RelocaTE3 classifies each TE as *non-reference* or
   *reference*, and for reference TEs reports **which** genome(s) carry it.
2. **Pangenome graph support**, later, against a graph built from the same
   genomes.

Phase 1 builds the TE-augmented reference genomes; Phase 2 builds the sample
and its reads. No pangenome graph is built yet.

## Design

- Output: `results/pangenome_panel/`
- Chromosome: Chr1
- Genomes: 4 MAGIC16 accessions spanning the cultivated diversity

  | Cultivar | Accession | Group |
  |---|---|---|
  | Nipponbare | GCA_001433935.1 | GJ-tmp (temperate japonica) |
  | Azucena | GCA_009830595.1 | GJ-trop1 (tropical japonica) |
  | IR 64 | GCA_009914875.1 | XI-1B1 (indica) |
  | N22 | GCA_001952365.3 | cA1 (aus) |

- Sharing patterns: all 15 non-empty subsets of the 4 genomes
- Events: 20 per pattern = 300; each genome carries 8 patterns = 160
- TE source: riceTElib (`rice7.0.0.liban`), reusing the 10 TE-group
  definitions from `config/riceTElib_benchmark.toml`

**IRGSP-1.0 is downloaded rather than reusing the local `MSU_r7.fa`.** They
are different Nipponbare assemblies; mixing one MSU genome with three MAGIC16
genomes would confound assembly differences with real biological divergence.

## Why orthology is a pipeline stage

In the riceTElib short- and long-read panels, every component genome was
Nipponbare, so a coordinate meant the same thing everywhere. Here the genomes
differ by megabases of indels. Placing a TE at Nipponbare 5,000,000 and IR64
5,000,000 creates **two unrelated loci**: a graph will not merge them,
RelocaTE3 will correctly report two insertions, and the truth table is wrong
with nothing raising an error.

So anchors are chosen in Nipponbare and lifted into every other genome:

```
minimap2 -x asm5 <target>.fa <Nipponbare>.fa > Nipponbare_to_<target>.paf
paftools.js liftover -l 50000 -q 5 <paf> candidate_anchors.bed
```

**Argument order is load-bearing.** minimap2 takes `<target> <query>`, and
`paftools.js liftover` maps *query* coordinates into *target* coordinates.
Anchors live in Nipponbare, so Nipponbare must be the **query**.

An anchor is accepted only if it lifts to **exactly one** position in **every**
genome. If fewer than 300 survive, reduce `events_per_pattern` — do not relax
the filter. A weaker filter produces a silently wrong truth set.

Background riceTElib-family TEs are annotated per genome, and sites within
`background_te_buffer` of a pre-existing element of the same family in *any*
genome are rejected, so a correct call is never scored as a false positive.

## Commands

Submit the whole Phase 1 chain:

```bash
bash pipeline/make_pangenome_panel/submit.sh
```

Stages individually (each is rerun-safe and skips on its sentinel):

```bash
sbatch pipeline/make_pangenome_panel/01_download_genomes.sh
sbatch pipeline/make_pangenome_panel/02_prepare_chromosomes.sh
sbatch pipeline/make_pangenome_panel/03_build_orthology_map.sh
sbatch --array=0-3 pipeline/make_pangenome_panel/04_annotate_background_tes.sh
sbatch pipeline/make_pangenome_panel/05_build_panel.sh
```

Check anchor yield without writing genomes:

```bash
module load minimap2/2.30
pixi run \
  --manifest-path /rhome/nmath020/bigdata/github/github_tools/data_sim/simulate_data/pyproject.toml \
  python pipeline/make_pangenome_panel/05_build_pangenome_panel.py \
  --config config/pangenome_panel.toml \
  --output results/pangenome_panel \
  --validate-only
```

Tests:

```bash
pixi run \
  --manifest-path /rhome/nmath020/bigdata/github/github_tools/data_sim/simulate_data/pyproject.toml \
  python -m unittest -v pipeline.make_pangenome_panel.test_build_pangenome_panel
```

## Output contract

```
results/pangenome_panel/
├── genomes/raw/<name>/           as downloaded, MD5-verified (read-only)
├── genomes/chr1/<name>.chr1.fa   PanSN-named: <name>#1#Chr1
├── genomes/augmented/            <name>.chr1.te.fa  <- the deliverable
├── orthology/
│   ├── Nipponbare_to_<name>.paf
│   ├── candidate_anchors.bed, lifted_<name>.bed
│   ├── alignment_stats.tsv
│   └── anchor_acceptance.tsv     <- headline diagnostic
├── annotation/<name>.background_te.bed
├── truth_events.tsv              one row per event, per-genome positions
├── sharing_matrix.tsv            long form: event_id, genome, present, position
├── run_metadata.json             SHA256 of every input and output
└── .panel.complete
```

`truth_events.tsv` keys each event by a 4-bit `pattern` over
(Nipponbare, Azucena, IR64, N22) plus a `genomes` column naming the carriers —
the direct answer key for the "which genome(s)" report. Absent genomes carry
`-` for position.

**PanSN naming** (`Nipponbare#1#Chr1`) is applied at extraction time, not
retrofitted, because renaming later would invalidate every coordinate.

## Rerun behavior

Every stage is sentinel-gated and safe to rerun. Downloads re-verify MD5 and
skip only when it still matches; a mismatch triggers a re-download rather than
silent acceptance.

## Caveats

- **IR64 (indica) gates anchor yield**, not N22: measured syntenic fractions
  are Azucena 0.902, N22 0.704, IR64 0.641. The deep indica/japonica split is
  the primary division in cultivated rice.
- **Anchor acceptance is 47%** (2779 of 5915), well above the 300 needed.
  `anchor_acceptance.tsv` is the headline diagnostic.
- **The accepted set is mildly optimistic, and it is measured** — see
  `orthology/anchor_bias.tsv`. Accepted anchors sit a median 123 bp from the
  nearest TE vs 20 bp for rejected, but they are not TE deserts: 37% lie
  inside an existing TE and 79% within 1 kb (vs 47% and 85% rejected). The
  substantive exclusion is structural (non-syntenic regions), not repeat
  context. Treat Phase 1 numbers as an upper bound.
- RepeatMasker against riceTElib is used instead of EDTA/panEDTA: Phase 1 only
  needs to know where riceTElib *families* already sit, and de-novo annotation
  of four genomes would dominate runtime.

## Phase 2: the sample

`results/pangenome_panel/sample/SampleA/` holds the genome whose reads a
caller actually sees. It is derived from the reference background and carries
three deliberately different event classes:

| Class | Count | Expected call | Answer key |
|---|---|---|---|
| `reference` | 150 | `reference` | the panel genomes carrying that event |
| `non_reference` | 100 | `non_reference` | — |
| `absent` | 150 | **no call** | over-calling control |

The `absent` class is the one that makes the benchmark honest: without it, a
caller that simply reports every panel event would score perfectly.

Both `reference` and `absent` are stratified across all 15 sharing patterns,
so attribution is tested for every pattern and no pattern can vanish into one
class.

A subtlety worth knowing: a sample event may be `reference` by virtue of a
genome *other* than the background it was derived from — e.g. an
Azucena-only TE placed at the orthologous Nipponbare locus. `truth_events.tsv`
records `-` for the non-carrying genome, so the sample builder takes those
coordinates from the anchor map instead (`reference_position_for`).

Build and simulate:

```bash
module load minimap2/2.30
pixi run --manifest-path /rhome/nmath020/bigdata/github/github_tools/data_sim/simulate_data/pyproject.toml   python pipeline/make_pangenome_panel/06_build_sample_genome.py   --config config/pangenome_panel.toml --output results/pangenome_panel

sbatch --array=0-2 pipeline/make_pangenome_panel/07_simulate_sample_reads.sh
# task 0 Illumina, 1 ONT-HQ, 2 PacBio HiFi
```

Sample outputs:

```
sample/SampleA/
├── SampleA.chr1.fa        PanSN: SampleA#1#Chr1
├── sample_truth.tsv       per-event expected_call + expected_genomes
├── run_metadata.json
├── .complete
└── reads/{illumina,ont-hq,hifi}/   all gzip compressed
```

## Next phase

- **Phase 3:** Minigraph-Cactus graph (`module load cactus/3.2.0` provides
  `cactus-pangenome`, `vg` 1.74.0, `minigraph` 0.21), vg Giraffe mapping, and
  GraffiTE as the comparison baseline.

## Toolkit goal

The generic parts of this pipeline — orthologous anchor selection via
liftover, sharing-pattern enumeration, multi-genome TE insertion, and SV
placement — are intended to migrate into `simulate_data` as a reusable
pangenome-simulation module. Species-specific choices (rice, Chr1, riceTElib)
are kept in `config/pangenome_panel.toml` rather than in code, so that move
stays cheap.
