# Design: multi-genome / pangenome TE benchmark panel

**Date:** 2026-08-08 America/Los_Angeles

**Current outcome (verified 2026-09-02):** This is the historical design
record. Phases 1-3 were completed in the linked operational workspace: the
four-genome panel, SampleA truth and Illumina/ONT-HQ/HiFi reads,
Minigraph-Cactus graph, long-read index, and all three Giraffe alignments have
completion sentinels.

## Purpose

Create simulated ground-truth data for two upcoming RelocaTE3 capabilities:

1. **Multi-FASTA support (near term).** A user supplies several reference
   genomes plus sample reads; RelocaTE3 reports each detected TE as
   *non-reference* (in no supplied genome) or *reference* (in at least one),
   and for reference TEs, **which** genome or genomes carry it.
2. **Pangenome graph support (later).** The same question answered against a
   graph built from those genomes.

This document originally scoped Phase 1 to the TE-augmented genomes. Later
work completed the sample/read and graph/alignment phases as recorded in
`docs/2026-08-08-pangenome-panel-workflow.md`.

## Decisions

| Decision | Choice | Rationale |
|---|---|---|
| Genomes | Nipponbare, Azucena, IR64, N22 | Spans temperate japonica, tropical japonica, indica, aus. Gives one closely related pair (Nip/Azucena) and diverged pairs, so detection can be measured against divergence. |
| Source | MAGIC16 / Rice Population Reference Panel | Near-T2T, >100x PacBio, uniformly assembled and annotated. Uniform provenance beats mixing sources. |
| Sharing design | All 15 non-empty subsets of 4 genomes | Exhaustively tests the "which genome(s)" report, including 3-of-4 patterns where partial-sharing bugs hide. |
| Events | 20 per pattern = 300 total | Each genome carries 8 of 15 patterns = 160 insertions. |
| TE source | riceTElib (`rice7.0.0.liban`, 2627 records) | Same library as the existing short/long-read panels, so results are comparable. |
| Scope | Chr1 only | Matches existing riceTElib panels; keeps graph construction and iteration fast during development. |
| Sample design | 5th Nipponbare-derived genome (Phase 2) | Only design giving ground truth for *both* output classes: reference (with attribution) and non-reference. |

## Verified facts (2026-08-08)

Established by execution, not assumption.

**Assemblies** — all four resolve on the NCBI FTP site (~120 MB gzipped each):

| Accession | Directory name | Cultivar | Group |
|---|---|---|---|
| GCA_001433935.1 | `GCA_001433935.1_IRGSP-1.0` | Nipponbare | GJ-tmp |
| GCA_009830595.1 | `GCA_009830595.1_AzucenaRS1` | Azucena | GJ-trop1 |
| GCA_009914875.1 | `GCA_009914875.1_OsIR64RS1` | IR 64 | XI-1B1 |
| GCA_001952365.3 | `GCA_001952365.3_OsN22RS2` | N22 | cA1 (aus) |

**Naming mismatch is real and must be handled.** MAGIC16 FASTA headers use
GenBank accessions:

```
>CM020633.1 Oryza sativa Japonica Group cultivar Azucena chromosome 1, ...
```

while the local `MSU_r7.fa` uses `Chr1..Chr12`. Chromosome identity should be
resolved from each assembly's `*_assembly_report.txt` (Sequence-Name /
GenBank-Accn / Assigned-Molecule), not by parsing the description text.

**Tooling is already installed on the cluster:**

| Tool | Source | Version |
|---|---|---|
| `cactus-pangenome` (Minigraph-Cactus) | `module load cactus/3.2.0` | 3.2.0 |
| `vg` | inside the cactus module | 1.74.0 |
| `minigraph` | inside the cactus module | 0.21-r606 |
| `minimap2`, `paftools.js`, `k8` | `module load minimap2/2.30` | 2.30 |
| EDTA | `module load edta/2.2.0` | 2.2.0 |
| RepeatMasker | `module load RepeatMasker/4.1.8` | 4.1.8 |
| seqkit, sratoolkit, samtools | modules | — |

No `vg`/`minigraph` install is needed; they ship inside the Cactus module.
`/bigdata` has 688 TB free.

**Liftover interface:**

```
paftools.js liftover [options] <aln.paf> <query.bed>
  -q INT    min mapping quality [5]
  -l INT    min alignment length [50000]
```

It lifts **query → target**, so the alignment must place Nipponbare as the
*query*: `minimap2 -x asm5 <other>.fa <nipponbare>.fa`. The defaults
(50 kb minimum block, MAPQ 5) usefully enforce the syntenic, unambiguous
regions we want.

## Tooling survey (literature)

**GraffiTE** (Nat Commun 2024) is the closest existing tool and the natural
comparison baseline. Its stack defines the field-standard approach:

| Stage | Tools |
|---|---|
| SV discovery, assembly→ref | minimap2 + SVIM-asm |
| SV discovery, long reads→ref | minimap2 + Sniffles2 |
| SV merging | SURVIVOR |
| TE annotation | RepeatMasker + OneCodeToFindThemAll, ≥80% coverage |
| Graph construction | vg |
| Short-read genotyping | PanGenie or vg Giraffe |
| Long-read genotyping | GraphAligner |

**GraffiTE degrades sharply on TE-rich plant genomes**: >80% recall/precision
on human, but 70–80% recall and 58–82% precision on maize (~80% TE content).
Rice is ~35–40% TE. That gap is RelocaTE3's opportunity, and the benchmark
should be built to expose it.

**Graph construction: use Minigraph-Cactus.** The 2025 GigaScience crop-plant
evaluation could not build mapping indexes for a PGGB graph of 6 sorghum
assemblies due to memory; Minigraph-Cactus was the balanced choice, with
Giraffe achieving >98% correct mapping. Minigraph alone misses <50 bp
variation. 28% of graph variants matched repeats, so graphs do represent TE
variation well.

**vg Giraffe now maps long reads** (2025) comparably to Minimap2+Sniffles2,
with the highest recall of methods tested — one graph mapper can serve both
the short-read and long-read panels.

**panEDTA** (2024) is the pangenome-aware EDTA, designed for consistent TE
annotation across multiple genomes.

## The orthology problem

This is what makes the multi-genome panel structurally different from every
previous panel.

In the riceTElib short- and long-read panels every component genome was
Nipponbare, so a coordinate meant the same thing in all of them. Here the
genomes differ by megabases of indels. Placing a TE at Nipponbare 5,000,000
and IR64 5,000,000 creates **two different loci**. The graph will not collapse
them into one bubble, RelocaTE3 will correctly report two separate insertions,
and the truth table will be wrong — silently, with no error raised.

**Solution.** Anchor sites are chosen in Nipponbare coordinates and lifted to
every other genome:

1. `minimap2 -x asm5 <other>Chr1.fa <nipponbare>Chr1.fa` → PAF (query = Nipponbare)
2. `paftools.js liftover` a BED of candidate Nipponbare anchors
3. Accept an anchor only if it lifts **1:1 and unambiguously to all four
   genomes**; reject multi-mapping, unlifted, or near-indel sites
4. Insert each event at each carrying genome's own lifted coordinate

Sites surviving this filter sit in syntenic, alignable regions, which is
exactly where a pangenome graph can represent them. The filter is expected to
reject a substantial fraction of candidates; that is correct behavior, not a
bug, and the acceptance rate should be reported.

## Background TE interference

The genomes already contain thousands of real riceTElib-family TEs. A
synthetic insertion designated "unique to Azucena" may land where IR64 already
carries a real element of the same family, which would make a correct
RelocaTE3 call look like a false positive.

Mitigation: annotate background TEs in all four genomes (panEDTA, or
RepeatMasker with riceTElib for speed), and reject anchor sites within a
configurable distance of any existing TE of the same family **in any genome**.
Report how many candidates this removes.

## Architecture

```
make_simulation_new/
├── config/
│   └── pangenome_panel.toml
├── pipeline/make_pangenome_panel/
│   ├── 01_download_genomes.sh       # download + checksum verify
│   ├── 02_prepare_chromosomes.sh    # Chr1 extraction, PanSN rename
│   ├── 03_build_orthology_map.sh    # minimap2 asm5 + paftools liftover
│   ├── 04_annotate_background_tes.sh
│   ├── 05_build_pangenome_panel.py  # anchors, sharing matrix, insertion
│   ├── submit.sh
│   ├── test_build_pangenome_panel.py
│   └── README.md
└── results/pangenome_panel/
    ├── genomes/raw/                 # as downloaded (read-only)
    ├── genomes/chr1/                # PanSN-named Chr1 per genome
    ├── genomes/augmented/           # TE-augmented output
    ├── orthology/                   # PAF + liftover maps + acceptance stats
    ├── annotation/                  # background TE annotation
    ├── truth_events.tsv
    ├── sharing_matrix.tsv
    ├── run_metadata.json
    └── .panel.complete
```

**PanSN naming is applied at extraction time**, not retrofitted:
`Nipponbare#1#Chr1`, `Azucena#1#Chr1`, `IR64#1#Chr1`, `N22#1#Chr1`. Graph
tools expect this, and renaming later would invalidate every coordinate.

## Truth schema

One row per event, keyed by a 4-bit presence mask over
(Nipponbare, Azucena, IR64, N22):

```
event_id  te_id  te_family  te_group  pattern  n_genomes
          nipponbare_pos  azucena_pos  ir64_pos  n22_pos
          strand  tsd  tsd_length  inserted_length
```

Absent genomes carry `-` for position. `sharing_matrix.tsv` gives the
pattern → genome-set expansion in long form, which is the direct answer key
for the "which genome(s)" report.

TE groups are assigned round-robin across sharing patterns so TE family and
sharing pattern stay uncorrelated; otherwise a family-specific bug would look
like an attribution bug.

## Phases

- **Phase 1 — complete.** Directories, download, Chr1 + PanSN, orthology map,
  background annotation, TE insertion → 4 augmented genomes + truth.
- **Phase 2 — complete.** 5th Nipponbare-derived sample genome carrying a mix of
  panel-matching insertions (→ reference, with known attribution) and novel
  insertions (→ non-reference), plus short/long read simulation reusing the
  existing panel machinery. Also deliberately *omits* some panel TEs, to test
  over-calling.
- **Phase 3 — complete.** Minigraph-Cactus graph and vg Giraffe mapping, with
  GraffiTE retained as the
  comparison baseline.

## Risks

1. **Anchor attrition.** Requiring a clean 1:1 lift to all four genomes may
   reject most candidates, especially near the Chr1 centromere. Mitigation:
   report acceptance rate; relax `-l` or fall back to 3-of-4 patterns only if
   300 events prove unreachable.
2. **N22 divergence.** Aus is the most diverged from japonica; its liftover
   yield will be the lowest and will gate the whole panel.
3. **Assembly chromosome naming** varies per accession; resolve from
   `*_assembly_report.txt`, never from header text.
4. **Background TE collisions** as described above.
5. **Phase 3 compute.** Minigraph-Cactus on 4 rice genomes is substantial even
   for Chr1; not a Phase 1 concern but sizes the later work.

## Next steps

The original implementation plan is
`docs/2026-08-08-pangenome-panel-implementation-plan.md`. The next work is to
connect RelocaTE3 and the standalone benchmark runner to the completed panel.
