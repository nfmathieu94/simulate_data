# Pangenome TE panel — Phase 1 workflow record

**Date:** 2026-08-08 America/Los_Angeles

## Purpose

Build four TE-augmented rice Chr1 genomes as ground truth for RelocaTE3
multi-FASTA support, and later pangenome-graph support.

## Current status: Phase 1 COMPLETE

300 events across all 15 sharing patterns; 160 insertions per genome. All
verification checks pass. **No reads simulated and no graph built**, per the
agreed scope.

| Stage | Job | Elapsed |
|---|---|---|
| Download (4 assemblies, MD5-verified) | 27293618 | 5:45 |
| Chr1 extraction + PanSN rename | 27293652 | 3:21 |
| Orthology (minimap2 asm5 + liftover) | 27293809 | 4:00 |
| Background TE annotation (RepeatMasker array) | 27293793 | 7:00–9:49 |
| Panel build | 27293818 | 0:25 |

## Results

**Chr1 lengths — the reason orthology is a pipeline stage:**

| Genome | Chr1 length | vs Nipponbare |
|---|---|---|
| Nipponbare | 43,270,923 | — |
| Azucena | 44,011,168 | +740 kb |
| IR64 | 44,350,042 | +1.08 Mb |
| N22 | 42,787,722 | −483 kb |

**Syntenic coverage (blocks ≥50 kb, fraction of Nipponbare Chr1):**

| Pair | Syntenic fraction |
|---|---|
| Nipponbare → Azucena | 0.9018 |
| Nipponbare → N22 | 0.7044 |
| Nipponbare → IR64 | 0.6411 |

**IR64 gates anchor yield, not N22** — the prediction in the design doc was
wrong. This is consistent with the deep indica/japonica (XI/GJ) split being
the primary division in cultivated rice; aus sits closer to japonica than
indica does by syntenic block coverage.

**Anchor acceptance — the headline diagnostic:**

```
candidates_sampled      5915
rejected_unlifted       3136
rejected_multimapping      0
rejected_background_te     0
accepted                2779   (47.0%)
needed                   300
```

47% acceptance yields 2779 usable anchors against 300 needed, so no
relaxation of the orthology filter was necessary. Zero multi-mapping
rejections indicates `paftools.js` defaults (`-l 50000 -q 5`) already exclude
ambiguous regions before our own check.

Zero background-TE rejections is a real result, not a dead filter: family
names were confirmed to match (262 of 300 truth families appear in the
annotation), and the filter was verified to reject inside a real element and
within its 1 kb buffer while passing when clear. With ~20 elements per family
over 43 Mb, expected collisions across 300 events are ≈1.

**Panel:**

| Genome | Insertions | Length before → after | Delta |
|---|---|---|---|
| Nipponbare | 160 | 43,270,923 → 43,545,687 | +274,764 |
| Azucena | 160 | 44,011,168 → 44,271,831 | +260,663 |
| IR64 | 160 | 44,350,042 → 44,601,971 | +251,929 |
| N22 | 160 | 42,787,722 → 43,069,227 | +281,505 |

Every delta equals the sum of inserted TE plus TSD lengths exactly.

Pattern distribution: 80 private, 120 pairwise, 80 triple, 20 core. All 10 TE
groups appear in every one of the 15 patterns, so TE family and sharing
pattern are fully decorrelated.

## Anchor bias measurement

The orthology filter keeps 47% of candidates, so the panel's realism depends
on how unrepresentative that 47% is. Measured against the Nipponbare
background annotation (`orthology/anchor_bias.tsv`):

| Metric (median) | Accepted | Rejected | Ratio |
|---|---|---|---|
| Distance to nearest TE (bp) | 123 | 20 | 6.15 |
| TE count in 20 kb window | 19 | 21 | 0.90 |
| TE-masked fraction of window | 0.318 | 0.400 | 0.79 |
| Position on Chr1 (Mb) | 24.9 | 19.7 | 1.26 |

| TE adjacency | Accepted | Rejected |
|---|---|---|
| Inside a TE | 37.0% | 47.3% |
| Within 100 bp of a TE | 48.1% | 57.2% |
| Within 1 kb of a TE | 79.3% | 85.0% |

**Interpretation.** The bias is real but modest. Accepted anchors are not TE
deserts: 37% sit inside an existing TE and 79% are within 1 kb of one, versus
47% and 85% for rejected sites. The 6.15x median-distance ratio comes from the
tail, not from wholesale exclusion of repeat-adjacent sites.

The substantive exclusion is **structural, not repeat context**: the rejected
53% are largely sites in non-syntenic regions, which also skew toward the Chr1
centromere (median 19.7 Mb vs 24.9 Mb).

Phase 1 therefore measures an optimistic bound, but a mildly optimistic one.
Report it as an upper bound; the real-TE dataset proposed below is what
brackets it from the realistic side.

## Phase 2 results (sample + reads)

`SampleA#1#Chr1`, derived from the Nipponbare background:
43,270,923 -> 43,715,474 bp. Length delta 444,551 bp, exactly the sum of
inserted TE plus TSD lengths.

| Class | Count | Expected call |
|---|---|---|
| `reference` | 150 | `reference`, attributed to the carrying panel genomes |
| `non_reference` | 100 | `non_reference` |
| `absent` | 150 | no call (over-calling control) |

Both `reference` and `absent` span all 15 sharing patterns.

**Reads (30x requested, all gzip):**

| Technology | Job elapsed | Reads | Delivered coverage |
|---|---|---|---|
| Illumina 150 bp PE | 12:02 | 4,370,384 pairs | 30.0x |
| ONT-HQ | 8:38 | 111,052 | 30.0x |
| PacBio HiFi | 6:10:11 | 79,459 | 27.2x |

HiFi delivers 27.2x rather than 30x because **ccs yield is 93.0%**
(79,459 of 85,441 ZMWs) — the same shortfall the long-read panel's yield
accounting was built for, reproducing at full scale. HiFi also costs ~43x the
ONT wall time (10 passes plus consensus polishing).

### Known waste: HiFi MAF intermediate

`pacbio_reads_0001.maf.gz` is **8.8 GB**, two thirds of the 13 GB sample read
total, because PBSIM3 records every one of the 10 subread passes. Nothing in
Phase 2 consumes it — sample truth is known by construction from
`sample_truth.tsv`, not derived from alignments. `07_simulate_sample_reads.sh`
should drop MAF/ref intermediates by default with a flag to retain them. The
existing files were left in place rather than deleted, since regenerating them
costs a 6-hour job.

## Failures / issues found

1. **`short` partition caps at 2 h**, not 4. Download script adjusted.
2. **RepeatMasker defaulted to the HMMER engine**, which cannot `hmmpress` a
   nucleotide library: `Error invoking hmmpress on rice7.0.0.liban`. Fixed
   with `-engine rmblast` (the engine GraffiTE also uses). RepeatMasker also
   drops `RM_*` working directories into the CWD, so it now runs from scratch.
3. **`paftools.js liftover` requires the `cg` (CIGAR) tag**, which
   `minimap2 -x asm5 --cs` does not emit. Failed with
   `unable to find the 'cg' tag`. Fixed by adding `-c`.
4. **`paftools.js liftover` discards the input BED name**, rewriting column 4
   as `<query_chrom>_<start>_<end>`. Anchors keyed by name matched nothing, so
   the first successful run reported 0/5915 accepted. Fixed by keying on the
   original reference coordinate. The unit test had encoded the assumed
   format, so it passed while the code was wrong — the test now uses real
   paftools output.
5. **The background-TE filter was dead code** on first write: it was tested
   but never called, because the TE family is not known until assignment. It
   now runs in the assignment loop.
6. **Sample build crashed on `int('-')`.** A panel event whose pattern excludes
   the reference (e.g. Azucena-only) has no reference position in
   `truth_events.tsv`. The anchor still has a coordinate in every genome, so
   `reference_position_for` now sources it from the anchor map. This is the
   case where a TE is "reference" by virtue of a genome other than the
   sample's own background — arguably the most interesting class in the
   benchmark, and it would have been silently dropped.

## Decisions / logic

- **IRGSP-1.0 downloaded rather than reusing local MSU_r7.fa.** Different
  Nipponbare assemblies; mixing sources would confound assembly differences
  with biological divergence.
- **Anchors chosen in Nipponbare and lifted**, accepted only on a 1:1 lift to
  all four genomes. Verified working: core event PGTE000281 sits at
  32,025,584 / 32,783,537 / 33,412,435 / 31,855,575 across the four genomes —
  naive same-coordinate placement would have been wrong by ~1.4 Mb.
- **RepeatMasker + riceTElib instead of EDTA/panEDTA.** Phase 1 only needs to
  know where riceTElib families already sit; de-novo annotation of four
  genomes would dominate runtime.
- **PanSN naming applied at extraction**, verified preserved through
  augmentation (`>Nipponbare#1#Chr1` etc.).

## Commands

```bash
# Whole chain
bash pipeline/make_pangenome_panel/submit.sh

# Anchor yield only, no genomes written
module load minimap2/2.30
pixi run --manifest-path /rhome/nmath020/bigdata/github/github_tools/data_sim/simulate_data/pyproject.toml \
  python pipeline/make_pangenome_panel/05_build_pangenome_panel.py \
  --config config/pangenome_panel.toml \
  --output results/pangenome_panel --validate-only

# Tests (25)
pixi run --manifest-path /rhome/nmath020/bigdata/github/github_tools/data_sim/simulate_data/pyproject.toml \
  python -m unittest -v pipeline.make_pangenome_panel.test_build_pangenome_panel
```

## Next steps

1. **Phase 2:** 5th Nipponbare-derived sample genome carrying a mix of
   panel-matching insertions (→ reference, with known attribution) and novel
   insertions (→ non-reference), deliberately omitting some panel TEs to test
   over-calling; then read simulation reusing the existing panel machinery.
2. **Wire RelocaTE3 multi-FASTA support** and score against
   `truth_events.tsv` / `sharing_matrix.tsv`.
3. **Phase 3:** Minigraph-Cactus graph (`module load cactus/3.2.0`), vg
   Giraffe mapping, GraffiTE as comparison baseline.
4. Consider scaling beyond Chr1 once the approach is validated; 2779 accepted
   anchors on Chr1 alone suggests headroom for a much larger event count. Note
   the exhaustive sharing design does not scale in genome count: patterns go as
   2^n - 1, so 8 genomes would need 255 patterns and 16 would need 65,535.
   Beyond 4 genomes, switch to sampled or phylogeny-structured patterns.
5. **Real-TE complement.** Hold out a real MAGIC16 genome, derive truth from
   assembly-to-assembly comparison (minimap2 + SVIM-asm + RepeatMasker, i.e.
   GraffiTE's assembly mode, which its own benchmarking found most reliable),
   and simulate reads from it. Real positions, nesting, truncation and family
   spectrum; high-confidence rather than exact truth. Brackets the synthetic
   panel from the realistic side.
6. **Toolkit goal.** The generic parts of this pipeline -- orthologous anchor
   selection, sharing-pattern enumeration, multi-genome TE insertion, SV
   placement -- are to migrate into `simulate_data` as a reusable pangenome
   simulation module, so this is not a rice-only script.
