# Pangenome TE Benchmark Panel — Phase 1 Implementation Plan

**Current outcome (verified 2026-09-02):** This historical Phase 1 plan was
implemented. Phases 2 and 3 were also completed later; see
`docs/2026-08-08-pangenome-panel-workflow.md` for realized results. The steps
below describe the original implementation sequence and are not an active
to-do list.

**Goal:** Produce four TE-augmented rice Chr1 genomes (Nipponbare, Azucena, IR64, N22) carrying 300 synthetic riceTElib insertions distributed across all 15 non-empty sharing patterns, with a truth table recording which genome(s) carry each event.

**Architecture:** Anchor sites are chosen in Nipponbare coordinates and lifted to every other genome with `minimap2 -x asm5` + `paftools.js liftover`, so a "shared" insertion sits at genuinely orthologous loci. Sites are rejected unless they lift 1:1 to all four genomes and avoid pre-existing TEs of the same family. Scripts live in `make_simulation_new/pipeline/make_pangenome_panel/`, reusing `simulate_data` helpers, matching the existing panel-builder pattern.

**Tech Stack:** Python 3.12 (via the simulate_data pixi env), Biopython, minimap2 2.30 + paftools.js, RepeatMasker 4.1.8 / EDTA 2.2.0, SLURM.

**Design doc:** `docs/2026-08-08-pangenome-panel-design.md`

---

## Repository paths

- **PANEL** = `/bigdata/wesslerlab/shared/Rice/Nathan/rice/make_simulated_genome/make_simulation_new` (NOT a git repo)
- **TOOLKIT** = `/bigdata/stajichlab/nmath020/github/github_tools/data_sim/simulate_data` (git repo)

Run Python through the toolkit env:

```bash
pixi run --manifest-path $TOOLKIT/pyproject.toml python <script>
```

Run the panel tests with:

```bash
cd $PANEL && pixi run --manifest-path $TOOLKIT/pyproject.toml \
  python -m unittest -v pipeline.make_pangenome_panel.test_build_pangenome_panel
```

---

## Critical domain facts

Verified 2026-08-08. Do not re-derive; do not guess past these.

1. **Assemblies** (NCBI FTP, ~120 MB gzipped each):

   | Accession | Directory | Cultivar |
   |---|---|---|
   | GCA_001433935.1 | `GCA_001433935.1_IRGSP-1.0` | Nipponbare |
   | GCA_009830595.1 | `GCA_009830595.1_AzucenaRS1` | Azucena |
   | GCA_009914875.1 | `GCA_009914875.1_OsIR64RS1` | IR 64 |
   | GCA_001952365.3 | `GCA_001952365.3_OsN22RS2` | N22 |

   URL shape:
   `https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/<3>/<3>/<3>/<dir>/<dir>_genomic.fna.gz`
   e.g. `.../GCA/009/830/595/GCA_009830595.1_AzucenaRS1/GCA_009830595.1_AzucenaRS1_genomic.fna.gz`

   Each directory also has `md5checksums.txt` and `<dir>_assembly_report.txt`.

2. **Chromosome naming differs per assembly.** Headers are GenBank accessions
   (`>CM020633.1 ... chromosome 1, ...`). Resolve chromosome identity from
   `*_assembly_report.txt` columns (Sequence-Name, GenBank-Accn,
   Assigned-Molecule) — **never** by parsing the description text.

3. **`paftools.js liftover <aln.paf> <query.bed>` lifts query → target.** To
   map Nipponbare anchors into Azucena, Nipponbare must be the *query*:
   `minimap2 -x asm5 azucena_chr1.fa nipponbare_chr1.fa > nip_to_azu.paf`.
   Defaults `-l 50000` (min block) and `-q 5` (MAPQ) enforce syntenic,
   unambiguous regions — keep them.

4. **PanSN naming** (`sample#haplotype#contig`, e.g. `Nipponbare#1#Chr1`) is
   applied at Chr1 extraction time. Graph tools expect it and retrofitting
   later invalidates every coordinate.

5. **Modules:** `cactus/3.2.0` (provides `cactus-pangenome`, `vg` 1.74.0,
   `minigraph` 0.21), `minimap2/2.30` (provides `paftools.js`, `k8`),
   `edta/2.2.0`, `RepeatMasker/4.1.8`, `seqkit/2.4.0`, `samtools/1.22.1`.

6. **TE library:** `input/TE_lib/other_TEs/riceTElib/rice7.0.0.liban`,
   2627 records. The 10 TE groups and their TSD rules are already defined in
   `config/riceTElib_benchmark.toml` — reuse those definitions verbatim.

7. **Local Nipponbare** `input/ref_genome/MSU_r7.fa` is MSU r7 with `Chr1`
   naming and Chr1 = 43,270,923 bp. This is a *different assembly* from
   IRGSP-1.0 (GCA_001433935.1). **Download IRGSP-1.0 rather than reusing
   MSU_r7**, so all four genomes come from one uniformly assembled panel.
   Note the discrepancy in the README.

---

## Task 1: Directory scaffold and configuration

**Files:**
- Create: `$PANEL/config/pangenome_panel.toml`
- Create: `$PANEL/pipeline/make_pangenome_panel/` (directory)

**Step 1: Write the config**

```toml
# Multi-genome TE panel for RelocaTE3 multi-FASTA and pangenome support.
# Phase 1: build TE-augmented genomes. Reads come in Phase 2.

[inputs]
te_library = "input/TE_lib/other_TEs/riceTElib/rice7.0.0.liban"
# TE group definitions (names, TSD rules, classifications) are reused from:
te_group_config = "config/riceTElib_benchmark.toml"

[panel]
chrom = "Chr1"
events_per_pattern = 20        # x 15 patterns = 300 events
seed = 917                     # distinct from the riceTElib panels (916)

[[genomes]]
name = "Nipponbare"
accession = "GCA_001433935.1"
directory = "GCA_001433935.1_IRGSP-1.0"
group = "GJ-tmp"
reference = true               # anchor coordinate system

[[genomes]]
name = "Azucena"
accession = "GCA_009830595.1"
directory = "GCA_009830595.1_AzucenaRS1"
group = "GJ-trop1"

[[genomes]]
name = "IR64"
accession = "GCA_009914875.1"
directory = "GCA_009914875.1_OsIR64RS1"
group = "XI-1B1"

[[genomes]]
name = "N22"
accession = "GCA_001952365.3"
directory = "GCA_001952365.3_OsN22RS2"
group = "cA1"

[orthology]
minimap2_preset = "asm5"
min_alignment_length = 50000   # paftools.js liftover -l
min_mapping_quality = 5        # paftools.js liftover -q
# Reject an anchor if the lifted position shifts by more than this relative
# to its neighbours, which indicates a nearby indel.
max_anchor_drift = 50

[selection]
min_distance = 5000            # between synthetic insertions
flank_acgt = 200               # unambiguous bases required around a site
min_te_length = 80
max_te_length = 20000
sense_strand_ratio = 0.5
# Reject a site within this distance of a pre-existing TE of the same family
# in ANY genome.
background_te_buffer = 1000

[download]
# NCBI FTP; each directory also carries md5checksums.txt
base_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA"
```

**Step 2: Verify it parses**

```bash
cd $PANEL && python3.12 -c "
import tomllib
c = tomllib.load(open('config/pangenome_panel.toml','rb'))
print([g['name'] for g in c['genomes']])
print('events:', c['panel']['events_per_pattern'] * 15)"
```
Expected: `['Nipponbare', 'Azucena', 'IR64', 'N22']` and `events: 300`

---

## Task 2: Genome download script

This is the "first script" — downloads four related rice genomes, verifies
checksums, and is safe to rerun.

**Files:**
- Create: `$PANEL/pipeline/make_pangenome_panel/01_download_genomes.sh`

**Requirements:**

- `#!/usr/bin/bash -l`, `set -euo pipefail`, date-stamped logging.
- SLURM header: `-p short`, `--mem=4gb`, `--cpus-per-task=2`,
  `--time=04:00:00`, `-o logs/pangenome_panel/download.%j.log`.
- For each genome in the config, download into
  `results/pangenome_panel/genomes/raw/<name>/`:
  - `<dir>_genomic.fna.gz`
  - `<dir>_assembly_report.txt`
  - `md5checksums.txt`
- **Verify the MD5** of the FASTA against `md5checksums.txt`. Fail loudly on
  mismatch; do not proceed with a corrupt genome.
- **Idempotent:** skip a genome whose `.download.complete` sentinel exists and
  whose MD5 still matches. Never re-download silently.
- Use `curl --fail --location --continue-at -` so partial downloads resume
  rather than truncating.
- Write `download_manifest.tsv`: name, accession, url, bytes, md5, timestamp.

**Verification:**

```bash
bash -n pipeline/make_pangenome_panel/01_download_genomes.sh   # syntax
sbatch pipeline/make_pangenome_panel/01_download_genomes.sh    # ~480 MB total
```

Expect four directories, four verified MD5s, and 12 chromosomes per assembly.

---

## Task 3: Chromosome extraction and PanSN renaming

**Files:**
- Create: `$PANEL/pipeline/make_pangenome_panel/02_prepare_chromosomes.sh`
- Create: `$PANEL/pipeline/make_pangenome_panel/prepare_chromosomes.py`
- Test: `$PANEL/pipeline/make_pangenome_panel/test_build_pangenome_panel.py`

**Step 1: Write the failing test**

```python
class TestChromosomeResolution(unittest.TestCase):
    """Chromosome identity comes from the assembly report, not header text."""

    REPORT = (
        "# Assembly name: AzucenaRS1\n"
        "# Sequence-Name\tSequence-Role\tAssigned-Molecule\t"
        "Assigned-Molecule-Location/Type\tGenBank-Accn\n"
        "1\tassembled-molecule\t1\tChromosome\tCM020633.1\n"
        "2\tassembled-molecule\t2\tChromosome\tCM020634.1\n"
    )

    def test_resolves_accession_for_chromosome(self):
        with tempfile.TemporaryDirectory() as tmp:
            report = Path(tmp) / "report.txt"
            report.write_text(self.REPORT)
            self.assertEqual(
                builder.accession_for_chromosome(report, "1"), "CM020633.1"
            )

    def test_missing_chromosome_raises(self):
        with tempfile.TemporaryDirectory() as tmp:
            report = Path(tmp) / "report.txt"
            report.write_text(self.REPORT)
            with self.assertRaises(ValueError):
                builder.accession_for_chromosome(report, "13")


class TestPanSN(unittest.TestCase):
    def test_pansn_name(self):
        self.assertEqual(builder.pansn_name("Azucena", "Chr1"), "Azucena#1#Chr1")
```

**Step 2: Implement**

`prepare_chromosomes.py` must:
- Parse `*_assembly_report.txt` into {Assigned-Molecule → GenBank-Accn}.
- Extract the requested chromosome from the gzipped FASTA by accession.
- Write `results/pangenome_panel/genomes/chr1/<name>.chr1.fa` with a single
  record named `<name>#1#Chr1`.
- Also write a `.fai` (via `samtools faidx`) and record the length.
- Refuse to overwrite; skip on a `.complete` sentinel.

**Step 3: Verify**

```bash
grep '^>' results/pangenome_panel/genomes/chr1/*.chr1.fa
cut -f1,2 results/pangenome_panel/genomes/chr1/*.fai
```
Expect four records named `<Name>#1#Chr1`. Nipponbare Chr1 ≈ 43.3 Mb; the
others will differ by megabases — that difference is the whole reason the
orthology stage exists.

---

## Task 4: Orthology map

**Files:**
- Create: `$PANEL/pipeline/make_pangenome_panel/03_build_orthology_map.sh`

**Requirements:**

- SLURM: `-p epyc`, `--mem=32gb`, `--cpus-per-task=16`, `--time=08:00:00`.
- For each non-reference genome, with **Nipponbare as the query**:

  ```bash
  minimap2 -x asm5 -t "$THREADS" --cs \
      "results/pangenome_panel/genomes/chr1/${TARGET}.chr1.fa" \
      "results/pangenome_panel/genomes/chr1/Nipponbare.chr1.fa" \
      > "results/pangenome_panel/orthology/Nipponbare_to_${TARGET}.paf"
  ```

  The query/target order is load-bearing: `paftools.js liftover` maps query
  coordinates into target coordinates, and anchors are chosen in Nipponbare.

- Report per pair: number of alignment blocks, total aligned bp, and the
  fraction of Nipponbare Chr1 covered by blocks ≥50 kb. Write to
  `orthology/alignment_stats.tsv`.
- Skip on a `.complete` sentinel.

**Sanity check before proceeding:** Azucena (japonica) should show the highest
syntenic coverage and N22 (aus) the lowest. If N22 coverage is very low, the
300-event target may be unreachable — see the fallback in Task 6.

---

## Task 5: Background TE annotation

**Files:**
- Create: `$PANEL/pipeline/make_pangenome_panel/04_annotate_background_tes.sh`

**Requirements:**

- SLURM: `-p epyc`, `--mem=32gb`, `--cpus-per-task=16`, `--time=24:00:00`.
- `module load RepeatMasker/4.1.8`.
- For each genome's Chr1, annotate with the **riceTElib library** (not the
  RepeatMasker default library) so families match the synthetic insertions:

  ```bash
  RepeatMasker -pa "$((THREADS / 4))" -lib "$TE_LIBRARY" -dir "$OUTDIR" \
      -nolow -no_is -gff "$CHR1_FASTA"
  ```

  `-pa` in RepeatMasker counts *jobs of 4 threads*, so divide.
- Convert each `.out` to BED with family names, into
  `annotation/<name>.background_te.bed`.
- Skip on a `.complete` sentinel.

Note: EDTA/panEDTA gives a better de-novo annotation but takes far longer.
RepeatMasker against riceTElib is the right Phase 1 choice because we only
need to know where *riceTElib families* already sit. Record this choice in
the README.

---

## Task 6: Anchor selection

**Files:**
- Create: `$PANEL/pipeline/make_pangenome_panel/05_build_pangenome_panel.py`
- Test: `test_build_pangenome_panel.py`

This is the correctness-critical stage.

**Step 1: Write the failing tests**

```python
class TestAnchorFiltering(unittest.TestCase):
    def test_anchor_must_lift_to_every_genome(self):
        lifted = {"Azucena": 100, "IR64": 200, "N22": None}
        self.assertFalse(builder.anchor_is_usable(lifted, ["Azucena", "IR64", "N22"]))

    def test_anchor_accepted_when_all_lift(self):
        lifted = {"Azucena": 100, "IR64": 200, "N22": 300}
        self.assertTrue(builder.anchor_is_usable(lifted, ["Azucena", "IR64", "N22"]))

    def test_multimapping_anchor_rejected(self):
        # paftools.js emits more than one interval for an ambiguous anchor
        self.assertFalse(builder.anchor_is_usable({"Azucena": "MULTI"}, ["Azucena"]))


class TestSharingPatterns(unittest.TestCase):
    def test_fifteen_non_empty_patterns(self):
        patterns = builder.sharing_patterns(4)
        self.assertEqual(len(patterns), 15)
        self.assertNotIn((False, False, False, False), patterns)

    def test_patterns_cover_every_cardinality(self):
        counts = collections.Counter(sum(p) for p in builder.sharing_patterns(4))
        self.assertEqual(counts, {1: 4, 2: 6, 3: 4, 4: 1})

    def test_each_genome_appears_in_eight_patterns(self):
        patterns = builder.sharing_patterns(4)
        for index in range(4):
            self.assertEqual(sum(1 for p in patterns if p[index]), 8)
```

**Step 2: Implement**

Anchor pipeline:
1. Sample candidate positions on Nipponbare Chr1, honouring `min_distance`,
   `flank_acgt`, and each TE group's TSD/context rule (reuse the existing
   `_site_is_usable` logic from `build_multite_panel.py` — do not reimplement).
2. Write candidates as BED in Nipponbare coordinates.
3. For each target genome, run
   `paftools.js liftover -l <min_len> -q <min_mapq> <paf> <candidates.bed>`.
4. Keep only anchors lifting to **exactly one** interval in **every** genome.
5. Reject anchors within `background_te_buffer` of a pre-existing TE of the
   same family in any genome.
6. Record acceptance statistics at each filter step into
   `orthology/anchor_acceptance.tsv`.

**Fallback if fewer than 300 anchors survive:** reduce
`events_per_pattern` rather than relaxing the orthology filter — a weaker
filter produces a silently wrong truth set, which is far worse than a smaller
panel. Report the shortfall explicitly.

---

## Task 7: Sharing matrix and event construction

**Files:**
- Modify: `05_build_pangenome_panel.py`

**Requirements:**

- Enumerate all 15 non-empty subsets of the 4 genomes.
- Assign `events_per_pattern` anchors to each pattern.
- Assign TE groups **round-robin across patterns**, so TE family and sharing
  pattern stay uncorrelated.
- Choose a distinct riceTElib exemplar per event (sample without replacement
  within a group, as `build_multite_panel.py` does).
- Draw TSD length and strand per the group's rule; the TSD is taken from
  **each carrying genome's own local sequence** at its lifted coordinate, not
  copied from Nipponbare — the flanking base composition differs per genome.

**Test:**

```python
def test_te_groups_are_uncorrelated_with_pattern(self):
    events = builder.build_events(settings, anchors, candidates)
    by_pattern = collections.defaultdict(set)
    for event in events:
        by_pattern[event.pattern].add(event.te_group)
    # every pattern should draw from several TE groups
    self.assertTrue(all(len(groups) > 1 for groups in by_pattern.values()))
```

---

## Task 8: Genome augmentation and truth output

**Files:**
- Modify: `05_build_pangenome_panel.py`

**Requirements:**

- For each genome, insert its events **in descending coordinate order** so
  earlier insertions do not shift later coordinates (or accumulate an offset,
  as `_write_component_genomes` does — either way, test it).
- Write `genomes/augmented/<name>.chr1.te.fa`, preserving PanSN naming.
- Write `truth_events.tsv` with the schema in the design doc, and
  `sharing_matrix.tsv` in long form (event_id, genome, present).
- Write `run_metadata.json` with SHA256 of every input and output genome,
  the config, and anchor acceptance statistics.
- Touch `.panel.complete`.

**Critical test — the whole point of the panel:**

```python
def test_shared_insertion_lands_at_orthologous_positions(self):
    """A shared event must sit at each genome's own lifted coordinate."""
    event = next(e for e in events if e.pattern == (True, True, False, False))
    self.assertNotEqual(event.positions["Nipponbare"], event.positions["Azucena"])
    self.assertEqual(
        event.positions["Azucena"], liftover_map["Azucena"][event.anchor]
    )

def test_absent_genomes_have_no_position(self):
    for event in events:
        for genome, present in zip(GENOMES, event.pattern):
            if not present:
                self.assertIsNone(event.positions[genome])

def test_augmented_genome_grew_by_expected_bases(self):
    """Length delta must equal the sum of inserted sequence plus TSDs."""
    expected = sum(e.inserted_length + e.tsd_length for e in nipponbare_events)
    self.assertEqual(augmented_length - original_length, expected)
```

---

## Task 9: Orchestration and documentation

**Files:**
- Create: `$PANEL/pipeline/make_pangenome_panel/submit.sh`
- Create: `$PANEL/pipeline/make_pangenome_panel/README.md`
- Create: `$PANEL/docs/2026-08-08-pangenome-panel-workflow.md`

`submit.sh` chains the stages with `--dependency=afterok`:
download → prepare → orthology → annotation → panel build.

README covers purpose, design, commands, output contract, rerun behavior, and
the MSU_r7-vs-IRGSP-1.0 note from fact 7.

The dated workflow doc follows the house format: date/time, purpose, current
status, commands, failures/issues, decisions/logic, next steps. Record the
anchor acceptance rate — it is the headline diagnostic for this panel.

---

## Task 10: Final verification

```bash
# All tests
cd $PANEL && pixi run --manifest-path $TOOLKIT/pyproject.toml \
  python -m unittest -v pipeline.make_pangenome_panel.test_build_pangenome_panel

# Truth integrity: 300 events, 15 patterns, 20 each
cut -f5 results/pangenome_panel/truth_events.tsv | tail -n +2 | sort | uniq -c

# Each genome carries 160 insertions (8 of 15 patterns x 20)
awk -F'\t' 'NR>1 && $6!="-"' results/pangenome_panel/sharing_matrix.tsv | wc -l

# Augmented genomes exist and are larger than their inputs
ls -la results/pangenome_panel/genomes/augmented/
```

**Historical Phase 1 stop condition:** the original task ended after the four
augmented genomes and truth table. Later tasks built the sample, reads, graph,
long-read index, and Giraffe alignments; RelocaTE3 integration remains future
work.
