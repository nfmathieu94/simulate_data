# Long-read TE Benchmark Panel Implementation Plan

**Current outcome (verified 2026-09-02):** This historical plan was implemented.
The statements below describing broken modules or missing tools record the
pre-implementation state. The operational riceTElib panel is 15 of 18 tasks
complete; all three HiFi 30x tasks reached their 48-hour limit and still need a
safe rerun.

**Goal:** Generate simulated ONT high-accuracy and PacBio HiFi read panels over the existing riceTElib truth set, so long-read support can be added to RelocaTE3 and scored against the same 500 insertion events as the short-read benchmark.

**Architecture:** Two layers. The reusable long-read engine (pbsim3 + ccs wrappers, MAF parsing, gzip pooling, the 8-component mixture loop) goes in the `simulate_data` toolkit, mirroring how `build_multite_panel.py` already imports `te_benchmark_panel._simulate_reads`. The panel-specific driver, TOML config, and SLURM scripts go in `make_simulation_new`. The existing riceTElib catalog is consumed strictly read-only, so the truth set is bit-identical across technologies.

**Tech Stack:** Python 3.12, pixi, PBSIM3 3.0.5, pbccs 6.4.0, pysam, Biopython, SLURM.

**Design doc:** `docs/2026-08-07-longread-te-benchmark-panel.md`

---

## Repository paths

Two repos are involved. Use absolute paths; they are on different filesystems.

- **TOOLKIT** = `/bigdata/stajichlab/nmath020/github/github_tools/data_sim/simulate_data`
  (git repo; implementation originally landed on
  `feat/longread-te-benchmark-panel`, now contained in
  `feat/pangenome-panel-design`)
- **PANEL** = `/bigdata/wesslerlab/shared/Rice/Nathan/rice/make_simulated_genome/make_simulation_new`
  (NOT a git repo — no commits there)

Run all toolkit tests with:

```bash
cd $TOOLKIT && pixi run pytest tests/ -v
```

---

## Critical domain facts

Read these before writing code. Each was verified by execution on 2026-08-07, not assumed.

1. **PBSIM3 model arguments are file paths, not names.** `--qshmm QSHMM-ONT` fails
   with `ERROR: Cannot open file: QSHMM-ONT`. It needs
   `$CONDA_PREFIX/data/QSHMM-ONT.model`. The seven shipped models are
   `QSHMM-{ONT,ONT-HQ,RSII}.model` and `ERRHMM-{ONT,ONT-HQ,RSII,SEQUEL}.model`.
   **This was why `reads_ont.py` and `reads_pacbio.py` did not work before this
   plan was implemented.**

2. **PBSIM3 writes one output set per contig**, numbered from 1:
   `<prefix>_0001.fq.gz`, `<prefix>_0001.maf.gz`, `<prefix>_0001.ref`.
   FASTQ is already gzipped. Reads are single-end. **There is no SAM output.**

3. **The MAF source is literally named `ref`, not the contig name.** A block looks
   like:

   ```
   a
   s ref 100726 58666 + 200000 AGCATTGG-CATCTAGG...
   s S1_1       0 57829 +  57829 AGCATTG-ACA--T--GC...
   ```

   Contig identity must be recovered from the sibling `<prefix>_NNNN.ref` file,
   whose first line is the FASTA header (`>testchr`). Do not assume the index
   order matches anything else; read the `.ref` file.

4. **MAF field layout** is `s <src> <start> <size> <strand> <srcSize> <seq>`, where
   `start` is 0-based and `size` counts non-gap characters. So the reference span
   is `[start, start + size)`, the direct analogue of pysam's
   `reference_start`/`reference_end`. If `strand` is `-`, the true start is
   `srcSize - start - size`.

5. **`--difference-ratio` default `6:55:39` is PacBio RS II and wrong for both our
   platforms.** PBSIM3 documents `39:24:36` for ONT and `22:45:33` for Sequel.

6. **HiFi requires a two-stage chain**: pbsim3 `--pass-num 10` emits a subread
   `.bam`, which `ccs` collapses into HiFi reads. The chain was tested end-to-end
   and works. **ccs yield is ~92.5%**, so requested depth overstates delivered
   depth — record actual yield and fail below a floor.

7. **Read names collide across components** (every component produces `S1_1`,
   `S1_2`, ...). The pooling prefix `<sample>:<component>:h<hap>:s<seed>` is what
   makes them unique, exactly as in the short-read panel.

8. **`src/simulate_data/cli.py` auto-discovers every module in `modules/`** as a
   subcommand. Shared helpers therefore go in `src/simulate_data/longread.py`
   (package root), NOT in `modules/`, or they would become a bogus subcommand.

9. **Compression is a hard requirement.** No uncompressed FASTQ may ever be
   written. Pool gzip-to-gzip. Use `GZIP_LEVEL = 6` to match
   `PANEL/pipeline/fastq_compression.py`.

---

## Task 1: Model path resolution

**Files:**
- Create: `$TOOLKIT/src/simulate_data/longread.py`
- Test: `$TOOLKIT/tests/test_longread.py`

**Step 1: Write the failing test**

```python
import os
import pytest
from pathlib import Path
from simulate_data.longread import resolve_model_path, MODEL_NAMES


def test_resolves_bare_name(tmp_path, monkeypatch):
    data = tmp_path / "data"
    data.mkdir()
    (data / "QSHMM-ONT-HQ.model").write_text("x")
    monkeypatch.setenv("CONDA_PREFIX", str(tmp_path))
    assert resolve_model_path("QSHMM-ONT-HQ") == data / "QSHMM-ONT-HQ.model"


def test_resolves_name_with_suffix(tmp_path, monkeypatch):
    data = tmp_path / "data"
    data.mkdir()
    (data / "ERRHMM-SEQUEL.model").write_text("x")
    monkeypatch.setenv("CONDA_PREFIX", str(tmp_path))
    assert resolve_model_path("ERRHMM-SEQUEL.model") == data / "ERRHMM-SEQUEL.model"


def test_explicit_existing_path_passes_through(tmp_path):
    model = tmp_path / "custom.model"
    model.write_text("x")
    assert resolve_model_path(str(model)) == model


def test_missing_model_lists_available(tmp_path, monkeypatch):
    data = tmp_path / "data"
    data.mkdir()
    (data / "QSHMM-ONT.model").write_text("x")
    monkeypatch.setenv("CONDA_PREFIX", str(tmp_path))
    with pytest.raises(FileNotFoundError, match="QSHMM-ONT"):
        resolve_model_path("NOPE")
```

**Step 2: Run test to verify it fails**

Run: `cd $TOOLKIT && pixi run pytest tests/test_longread.py -v`
Expected: FAIL, `ModuleNotFoundError: No module named 'simulate_data.longread'`

**Step 3: Write minimal implementation**

```python
"""Shared helpers for long-read simulation with PBSIM3 and ccs.

Lives at the package root rather than under ``modules/`` because
``cli.py`` auto-discovers every module in that package as a subcommand.
"""

from __future__ import annotations

import os
from pathlib import Path

MODEL_NAMES = (
    "QSHMM-ONT",
    "QSHMM-ONT-HQ",
    "QSHMM-RSII",
    "ERRHMM-ONT",
    "ERRHMM-ONT-HQ",
    "ERRHMM-RSII",
    "ERRHMM-SEQUEL",
)


def _model_dir() -> Path:
    prefix = os.environ.get("CONDA_PREFIX")
    if not prefix:
        raise RuntimeError(
            "CONDA_PREFIX is unset; run inside the pixi environment so PBSIM3 "
            "model files can be located"
        )
    return Path(prefix) / "data"


def resolve_model_path(model: str) -> Path:
    """Resolve a PBSIM3 model name or path to an existing model file.

    PBSIM3 requires a path to a .model file; a bare name such as
    ``QSHMM-ONT`` fails with 'Cannot open file'.
    """
    candidate = Path(model)
    if candidate.is_file():
        return candidate

    directory = _model_dir()
    name = model if model.endswith(".model") else f"{model}.model"
    resolved = directory / name
    if resolved.is_file():
        return resolved

    available = sorted(p.stem for p in directory.glob("*.model"))
    raise FileNotFoundError(
        f"PBSIM3 model not found: {model}. Looked for {resolved}. "
        f"Available models: {', '.join(available) or '(none)'}"
    )
```

**Step 4: Run test to verify it passes**

Run: `cd $TOOLKIT && pixi run pytest tests/test_longread.py -v`
Expected: 4 passed

**Step 5: Commit**

```bash
cd $TOOLKIT
git add src/simulate_data/longread.py tests/test_longread.py
git commit -m "feat(longread): resolve PBSIM3 model names to file paths"
```

---

## Task 2: Platform presets and pbsim command construction

**Files:**
- Modify: `$TOOLKIT/src/simulate_data/longread.py`
- Test: `$TOOLKIT/tests/test_longread.py`

**Step 1: Write the failing test**

```python
from simulate_data.longread import PLATFORMS, build_pbsim_command


def test_ont_hq_uses_qshmm_and_ont_difference_ratio(tmp_path, monkeypatch):
    data = tmp_path / "data"
    data.mkdir()
    (data / "QSHMM-ONT-HQ.model").write_text("x")
    monkeypatch.setenv("CONDA_PREFIX", str(tmp_path))
    cmd = build_pbsim_command(
        PLATFORMS["ont-hq"], tmp_path / "g.fa", 10.0, tmp_path / "out", seed=7
    )
    assert cmd[0] == "pbsim"
    assert "--method" in cmd and cmd[cmd.index("--method") + 1] == "qshmm"
    assert cmd[cmd.index("--qshmm") + 1] == str(data / "QSHMM-ONT-HQ.model")
    assert cmd[cmd.index("--difference-ratio") + 1] == "39:24:36"
    assert cmd[cmd.index("--seed") + 1] == "7"
    assert "--pass-num" not in cmd


def test_hifi_uses_errhmm_and_multipass(tmp_path, monkeypatch):
    data = tmp_path / "data"
    data.mkdir()
    (data / "ERRHMM-SEQUEL.model").write_text("x")
    monkeypatch.setenv("CONDA_PREFIX", str(tmp_path))
    cmd = build_pbsim_command(
        PLATFORMS["hifi"], tmp_path / "g.fa", 10.0, tmp_path / "out", seed=7
    )
    assert cmd[cmd.index("--method") + 1] == "errhmm"
    assert cmd[cmd.index("--errhmm") + 1] == str(data / "ERRHMM-SEQUEL.model")
    assert cmd[cmd.index("--difference-ratio") + 1] == "22:45:33"
    assert cmd[cmd.index("--pass-num") + 1] == "10"


def test_depth_is_formatted_without_exponent(tmp_path, monkeypatch):
    data = tmp_path / "data"
    data.mkdir()
    (data / "QSHMM-ONT-HQ.model").write_text("x")
    monkeypatch.setenv("CONDA_PREFIX", str(tmp_path))
    cmd = build_pbsim_command(
        PLATFORMS["ont-hq"], tmp_path / "g.fa", 0.0000005, tmp_path / "out", seed=1
    )
    assert "e-" not in cmd[cmd.index("--depth") + 1]
```

**Step 2: Run test to verify it fails**

Run: `cd $TOOLKIT && pixi run pytest tests/test_longread.py -k platform -v`
Expected: FAIL, `ImportError: cannot import name 'PLATFORMS'`

**Step 3: Write minimal implementation**

Append to `longread.py`:

```python
from dataclasses import dataclass


@dataclass(frozen=True)
class Platform:
    """A long-read sequencing platform preset."""

    key: str
    method: str            # "qshmm" or "errhmm"
    model: str
    difference_ratio: str
    length_mean: int
    length_sd: int
    accuracy_mean: float
    pass_num: int          # 1 = single pass; >1 requires ccs
    needs_ccs: bool


# Difference ratios are PBSIM3's documented per-platform values. The tool
# default of 6:55:39 is PacBio RS II and wrong for both platforms here.
PLATFORMS = {
    "ont-hq": Platform(
        key="ont-hq",
        method="qshmm",
        model="QSHMM-ONT-HQ",
        difference_ratio="39:24:36",
        length_mean=12000,
        length_sd=9000,
        accuracy_mean=0.99,
        pass_num=1,
        needs_ccs=False,
    ),
    "hifi": Platform(
        key="hifi",
        method="errhmm",
        model="ERRHMM-SEQUEL",
        difference_ratio="22:45:33",
        length_mean=15000,
        length_sd=2000,
        accuracy_mean=0.999,
        pass_num=10,
        needs_ccs=True,
    ),
}


def build_pbsim_command(
    platform: Platform,
    genome: Path,
    depth: float,
    prefix: Path,
    seed: int,
) -> list[str]:
    """Build the PBSIM3 command for one component genome."""
    model_path = resolve_model_path(platform.model)
    cmd = [
        "pbsim",
        "--strategy", "wgs",
        "--method", platform.method,
        f"--{platform.method}", str(model_path),
        "--genome", str(genome),
        # f-format avoids scientific notation, which PBSIM3 will not parse.
        "--depth", f"{depth:.6f}",
        "--length-mean", str(platform.length_mean),
        "--length-sd", str(platform.length_sd),
        "--accuracy-mean", str(platform.accuracy_mean),
        "--difference-ratio", platform.difference_ratio,
        "--prefix", str(prefix),
        "--seed", str(seed),
    ]
    if platform.pass_num > 1:
        cmd.extend(["--pass-num", str(platform.pass_num)])
    return cmd
```

**Step 4: Run test to verify it passes**

Run: `cd $TOOLKIT && pixi run pytest tests/test_longread.py -v`
Expected: 7 passed

**Step 5: Commit**

```bash
cd $TOOLKIT
git add src/simulate_data/longread.py tests/test_longread.py
git commit -m "feat(longread): add ONT-HQ and HiFi platform presets"
```

---

## Task 3: MAF parsing with contig identity

**Files:**
- Modify: `$TOOLKIT/src/simulate_data/longread.py`
- Test: `$TOOLKIT/tests/test_longread.py`

**Step 1: Write the failing test**

```python
import gzip
from simulate_data.longread import contig_name_for, iter_maf_alignments


def _write_maf(path, blocks):
    with gzip.open(path, "wt") as fh:
        for start, size, strand, src_size in blocks:
            fh.write("a\n")
            fh.write(f"s ref {start} {size} {strand} {src_size} ACGT\n")
            fh.write(f"s S1_1 0 {size} + {size} ACGT\n\n")


def test_contig_name_comes_from_ref_file(tmp_path):
    (tmp_path / "p_0001.ref").write_text(">Chr1 some description\nACGT\n")
    assert contig_name_for(tmp_path / "p_0001.ref") == "Chr1"


def test_iter_maf_yields_forward_spans(tmp_path):
    maf = tmp_path / "p_0001.maf.gz"
    _write_maf(maf, [(100, 50, "+", 1000)])
    assert list(iter_maf_alignments(maf)) == [(100, 150)]


def test_iter_maf_converts_reverse_strand_coordinates(tmp_path):
    maf = tmp_path / "p_0001.maf.gz"
    # start is measured from the reverse strand: true start = 1000 - 100 - 50
    _write_maf(maf, [(100, 50, "-", 1000)])
    assert list(iter_maf_alignments(maf)) == [(850, 900)]


def test_iter_maf_ignores_read_lines(tmp_path):
    maf = tmp_path / "p_0001.maf.gz"
    _write_maf(maf, [(0, 10, "+", 100), (20, 10, "+", 100)])
    assert list(iter_maf_alignments(maf)) == [(0, 10), (20, 30)]
```

**Step 2: Run test to verify it fails**

Run: `cd $TOOLKIT && pixi run pytest tests/test_longread.py -k maf -v`
Expected: FAIL, `ImportError: cannot import name 'iter_maf_alignments'`

**Step 3: Write minimal implementation**

```python
import gzip
from typing import Iterator


def contig_name_for(ref_path: Path) -> str:
    """Return the contig name recorded in a PBSIM3 ``.ref`` file.

    PBSIM3 names the MAF source 'ref' rather than the contig, so contig
    identity has to be recovered from the sibling .ref file's header.
    """
    with ref_path.open() as handle:
        header = handle.readline().strip()
    if not header.startswith(">"):
        raise ValueError(f"Not a FASTA header in {ref_path}: {header!r}")
    return header[1:].split()[0]


def iter_maf_alignments(maf_path: Path) -> Iterator[tuple[int, int]]:
    """Yield ``(start, end)`` reference spans from a PBSIM3 MAF file.

    Fields are ``s <src> <start> <size> <strand> <srcSize> <seq>`` with a
    0-based start and a size that counts only non-gap characters. Only the
    first 's' line of each block (the reference) is used.
    """
    with gzip.open(maf_path, "rt") as handle:
        expect_reference = False
        for line in handle:
            if line.startswith("a"):
                expect_reference = True
                continue
            if not expect_reference or not line.startswith("s "):
                continue
            fields = line.split()
            start, size, strand, src_size = (
                int(fields[2]),
                int(fields[3]),
                fields[4],
                int(fields[5]),
            )
            if strand == "-":
                start = src_size - start - size
            yield start, start + size
            expect_reference = False
```

**Step 4: Run test to verify it passes**

Run: `cd $TOOLKIT && pixi run pytest tests/test_longread.py -v`
Expected: 11 passed

**Step 5: Commit**

```bash
cd $TOOLKIT
git add src/simulate_data/longread.py tests/test_longread.py
git commit -m "feat(longread): parse PBSIM3 MAF spans with contig identity"
```

---

## Task 4: Gzip-to-gzip read pooling

**Files:**
- Modify: `$TOOLKIT/src/simulate_data/longread.py`
- Test: `$TOOLKIT/tests/test_longread.py`

Compression is a hard requirement: this never materializes an uncompressed FASTQ.

**Step 1: Write the failing test**

```python
from simulate_data.longread import GZIP_LEVEL, pool_fastq_gz


def test_pooling_prefixes_names_and_counts_reads(tmp_path):
    source = tmp_path / "in.fq.gz"
    with gzip.open(source, "wt") as fh:
        fh.write("@S1_1\nACGT\n+\nIIII\n@S1_2 note\nTTTT\n+\nIIII\n")
    out = tmp_path / "out.fastq.gz"
    with gzip.open(out, "wt", compresslevel=GZIP_LEVEL) as dest:
        assert pool_fastq_gz(source, dest, "sampleA:baseline:h1:s42") == 2
    text = gzip.open(out, "rt").read()
    assert text.startswith("@sampleA:baseline:h1:s42:S1_1\n")
    assert "@sampleA:baseline:h1:s42:S1_2 note\n" in text
    assert "ACGT\n+\nIIII\n" in text


def test_pooling_rejects_truncated_fastq(tmp_path):
    source = tmp_path / "in.fq.gz"
    with gzip.open(source, "wt") as fh:
        fh.write("@S1_1\nACGT\n+\n")
    out = tmp_path / "out.fastq.gz"
    with gzip.open(out, "wt", compresslevel=GZIP_LEVEL) as dest:
        with pytest.raises(ValueError, match="Truncated FASTQ"):
            pool_fastq_gz(source, dest, "p")
```

**Step 2: Run test to verify it fails**

Run: `cd $TOOLKIT && pixi run pytest tests/test_longread.py -k pool -v`
Expected: FAIL, `ImportError: cannot import name 'pool_fastq_gz'`

**Step 3: Write minimal implementation**

```python
# Matches PANEL/pipeline/fastq_compression.py so both panels compress alike.
GZIP_LEVEL = 6


def pool_fastq_gz(source: Path, destination, prefix: str) -> int:
    """Append a gzipped FASTQ into an open gzip stream, prefixing read names.

    Read names collide across components (every PBSIM3 run emits S1_1,
    S1_2, ...), so the prefix is what makes pooled names unique.
    """
    count = 0
    with gzip.open(source, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq, plus, qual = (
                handle.readline(),
                handle.readline(),
                handle.readline(),
            )
            if not qual:
                raise ValueError(f"Truncated FASTQ: {source}")
            token, *rest = header.rstrip().split(maxsplit=1)
            destination.write(
                f"@{prefix}:{token[1:]}" + (f" {rest[0]}" if rest else "") + "\n"
            )
            destination.write(seq)
            destination.write(plus)
            destination.write(qual)
            count += 1
    return count
```

**Step 4: Run test to verify it passes**

Run: `cd $TOOLKIT && pixi run pytest tests/test_longread.py -v`
Expected: 13 passed

**Step 5: Commit**

```bash
cd $TOOLKIT
git add src/simulate_data/longread.py tests/test_longread.py
git commit -m "feat(longread): pool gzipped FASTQs without decompressing to disk"
```

---

## Task 5: Add pbccs to the environment

**Files:**
- Modify: `$TOOLKIT/pyproject.toml`

**Step 1: Add the dependency**

In `[tool.pixi.dependencies]`, after the `pbsim3` line:

```toml
pbccs = ">=6.4.0,<7"
```

**Step 2: Resolve and verify**

```bash
cd $TOOLKIT && pixi install && pixi run ccs --version
```
Expected: `ccs 6.4.0 (commit v6.4.0)`

**Step 3: Commit**

```bash
cd $TOOLKIT
git add pyproject.toml pixi.lock
git commit -m "build: add pbccs for PacBio HiFi consensus generation"
```

---

## Task 6: ccs wrapper with yield accounting

**Files:**
- Modify: `$TOOLKIT/src/simulate_data/longread.py`
- Test: `$TOOLKIT/tests/test_longread.py`

ccs yield is ~92.5%, so delivered HiFi depth is below requested. Silent
under-coverage would look like a caller recall problem, so it is measured and
gated.

**Step 1: Write the failing test**

```python
from simulate_data.longread import parse_ccs_report

REPORT = """ZMWs input               : 40

ZMWs pass filters        : 37 (92.50%)
ZMWs fail filters        : 3 (7.500%)
"""


def test_parse_ccs_report_extracts_yield(tmp_path):
    path = tmp_path / "r.ccs_report.txt"
    path.write_text(REPORT)
    stats = parse_ccs_report(path)
    assert stats["zmws_input"] == 40
    assert stats["zmws_passed"] == 37
    assert abs(stats["yield_fraction"] - 0.925) < 1e-9


def test_parse_ccs_report_handles_zero_input(tmp_path):
    path = tmp_path / "r.ccs_report.txt"
    path.write_text("ZMWs input               : 0\nZMWs pass filters        : 0 (0%)\n")
    assert parse_ccs_report(path)["yield_fraction"] == 0.0
```

**Step 2: Run test to verify it fails**

Run: `cd $TOOLKIT && pixi run pytest tests/test_longread.py -k ccs -v`
Expected: FAIL, `ImportError: cannot import name 'parse_ccs_report'`

**Step 3: Write minimal implementation**

```python
import re

from simulate_data.utils import check_tool_installed, run_command


def parse_ccs_report(path: Path) -> dict:
    """Extract ZMW yield statistics from a ccs report."""
    text = path.read_text()
    def _find(label: str) -> int:
        match = re.search(rf"{label}\s*:\s*(\d+)", text)
        if match is None:
            raise ValueError(f"Missing '{label}' in ccs report: {path}")
        return int(match.group(1))

    total = _find("ZMWs input")
    passed = _find("ZMWs pass filters")
    return {
        "zmws_input": total,
        "zmws_passed": passed,
        "yield_fraction": passed / total if total else 0.0,
    }


def run_ccs(bam: Path, output_fastq_gz: Path, threads: int, min_yield: float) -> dict:
    """Collapse a PBSIM3 subread BAM into HiFi reads with ccs.

    ccs writes gzipped FASTQ directly, so nothing is decompressed to disk.
    """
    check_tool_installed("ccs")
    run_command(
        [
            "ccs",
            "--num-threads", str(threads),
            str(bam),
            str(output_fastq_gz),
        ]
    )
    report = output_fastq_gz.with_suffix("").with_suffix(".ccs_report.txt")
    if not report.is_file():
        # ccs derives the report name from the output stem.
        candidates = sorted(output_fastq_gz.parent.glob("*.ccs_report.txt"))
        if not candidates:
            raise FileNotFoundError(f"ccs report not found beside {output_fastq_gz}")
        report = candidates[-1]
    stats = parse_ccs_report(report)
    if stats["yield_fraction"] < min_yield:
        raise RuntimeError(
            f"ccs yield {stats['yield_fraction']:.3f} is below the floor "
            f"{min_yield:.3f} for {bam}; delivered HiFi depth would be far "
            "under the requested depth"
        )
    return stats
```

**Step 4: Run test to verify it passes**

Run: `cd $TOOLKIT && pixi run pytest tests/test_longread.py -v`
Expected: 15 passed

**Step 5: Commit**

```bash
cd $TOOLKIT
git add src/simulate_data/longread.py tests/test_longread.py
git commit -m "feat(longread): wrap ccs with yield accounting and a floor"
```

---

## Task 7: Repair the reads-ont module

**Files:**
- Modify: `$TOOLKIT/src/simulate_data/modules/reads_ont.py`
- Test: `$TOOLKIT/tests/test_reads_ont.py`

At the start of this plan, this module was broken because it passed a bare
model name. The implementation and tests were subsequently corrected.

**Step 1: Read the existing tests**

```bash
cd $TOOLKIT && cat tests/test_reads_ont.py
```

**Step 2: Write the failing test**

```python
def test_ont_command_uses_resolved_model_path(tmp_path, monkeypatch):
    data = tmp_path / "data"
    data.mkdir()
    (data / "QSHMM-ONT.model").write_text("x")
    monkeypatch.setenv("CONDA_PREFIX", str(tmp_path))
    cmd = _build_pbsim3_ont_command(
        ref_fasta=tmp_path / "ref.fa",
        coverage=5.0,
        output_prefix=str(tmp_path / "out"),
        error_model="QSHMM-ONT",
    )
    assert cmd[cmd.index("--qshmm") + 1] == str(data / "QSHMM-ONT.model")
    assert cmd[cmd.index("--difference-ratio") + 1] == "39:24:36"
```

**Step 3: Run test to verify it fails**

Run: `cd $TOOLKIT && pixi run pytest tests/test_reads_ont.py -v`
Expected: FAIL — the command still contains the bare string `QSHMM-ONT`

**Step 4: Implement**

In `_build_pbsim3_ont_command`, replace the bare `error_model` with
`str(resolve_model_path(error_model))` and add a `--difference-ratio` argument
defaulting to `39:24:36`. Add a `--difference-ratio` CLI flag in
`register_parser` with the same default.

**Step 5: Run tests**

Run: `cd $TOOLKIT && pixi run pytest tests/test_reads_ont.py -v`
Expected: PASS

**Step 6: Commit**

```bash
cd $TOOLKIT
git add src/simulate_data/modules/reads_ont.py tests/test_reads_ont.py
git commit -m "fix(reads-ont): resolve model path and use ONT difference ratio

The module passed a bare model name, which PBSIM3 rejects with
'Cannot open file'. It had never been run successfully."
```

---

## Task 8: Repair the reads-pacbio module and add a real HiFi path

**Files:**
- Modify: `$TOOLKIT/src/simulate_data/modules/reads_pacbio.py`
- Test: `$TOOLKIT/tests/test_reads_pacbio.py`

Same model-path bug. Additionally, `--read-type HiFi` currently just passes
`--pass-num` and stops, leaving an unusable subread BAM. It must run ccs.

**Step 1: Write the failing tests**

```python
def test_pacbio_command_uses_resolved_model_path(tmp_path, monkeypatch):
    data = tmp_path / "data"
    data.mkdir()
    (data / "QSHMM-RSII.model").write_text("x")
    monkeypatch.setenv("CONDA_PREFIX", str(tmp_path))
    cmd = _build_pbsim3_pacbio_command(
        ref_fasta=tmp_path / "ref.fa",
        coverage=5.0,
        output_prefix=str(tmp_path / "out"),
        error_model="QSHMM-RSII",
    )
    assert cmd[cmd.index("--qshmm") + 1] == str(data / "QSHMM-RSII.model")


def test_hifi_requires_multiple_passes(tmp_path, monkeypatch):
    with pytest.raises(ValueError, match="pass-num"):
        _validate_hifi(pass_num=1)
```

**Step 2: Run to verify failure, then implement**

Resolve the model path; when `read_type == "HiFi"`, require `pass_num > 1`,
select the `errhmm` method with `ERRHMM-SEQUEL`, use difference ratio
`22:45:33`, and after pbsim call `longread.run_ccs` on each per-contig `.bam`.

**Step 3: Run tests and commit**

```bash
cd $TOOLKIT && pixi run pytest tests/test_reads_pacbio.py -v
git add src/simulate_data/modules/reads_pacbio.py tests/test_reads_pacbio.py
git commit -m "fix(reads-pacbio): resolve model path; generate real HiFi via ccs"
```

---

## Task 9: Long-read mixture simulation in the panel module

**Files:**
- Modify: `$TOOLKIT/src/simulate_data/modules/te_benchmark_panel.py`
- Test: `$TOOLKIT/tests/test_te_benchmark_panel.py`

This is the core. It mirrors `_simulate_reads` (line 591) exactly: the same
`COMPONENTS` loop, the same `_seed()` derivation, the same manifest and QC
schema. Differences: single-end, MAF-based QC, gzip output.

**Step 1: Write the failing test**

```python
def test_long_read_support_counting_uses_maf(tmp_path):
    """MAF spans should drive the same present/absent tally as SAM did."""
    counts = defaultdict(lambda: [0, 0])
    rows = [
        {"chrom": "Chr1", "coordinate": 1000, "event_id": "TE000001",
         "present": 1, "haplotype": 1},
    ]
    truth = {"TE000001": {"inserted_length": "500", "tsd_length": "5"}}
    # A read spanning 0-20000 covers both junctions of the insertion.
    _count_origin_support_maf([("Chr1", 0, 20000)], rows, truth, counts)
    assert counts["TE000001"][0] == 1


def test_long_read_seed_matches_short_read_derivation():
    """Seeds must stay comparable across technologies."""
    assert _seed(916, 30.0, 1, 0, 1) == _seed(916, 30.0, 1, 0, 1)
```

**Step 2: Run to verify failure**

Run: `cd $TOOLKIT && pixi run pytest tests/test_te_benchmark_panel.py -k long_read -v`
Expected: FAIL, `NameError: _count_origin_support_maf`

**Step 3: Implement**

Add two functions to `te_benchmark_panel.py`:

`_count_origin_support_maf(spans, coordinate_rows, truth, counts, margin=5)` —
identical bisect logic to `_count_origin_support` (line 483) but taking
`(chrom, start, end)` tuples instead of a SAM path.

`_simulate_long_reads(settings, output, coverage, replicate, platform_key)` —
structure copied from `_simulate_reads`:

1. Assert `.catalog.complete` exists; `check_tool_installed("pbsim")`.
2. `sample = f"cov{coverage:g}x_rep{replicate}"`; sample dir is
   `output / "reads" / platform_key / sample`. Skip on `.complete`; refuse a
   non-empty dir without it.
3. Open `reads.fastq.gz` via `gzip.open(..., "wt", compresslevel=GZIP_LEVEL)`.
4. For each of the 8 `(component, hap)` pairs: depth `coverage * weight / 2`,
   `seed = _seed(...)` (unchanged), run `build_pbsim_command`, then for each
   per-contig output index:
   - `contig = contig_name_for(prefix_NNNN.ref)`
   - HiFi: `run_ccs(prefix_NNNN.bam, hifi.fastq.gz, threads, min_yield)`;
     record yield in the manifest row
   - accumulate MAF spans as `(contig, start, end)` for QC
   - `pool_fastq_gz(fq_gz, dest, f"{sample}:{component}:h{hap}:s{seed}")`
5. Reference-only control into
   `output / "reads" / platform_key / f"{sample}_reference_control"`.
6. Write `component_manifest.tsv` (adding `ccs_yield` and `reads` columns),
   then reuse the existing `_write_observed_support` unchanged.
7. Touch both `.complete` sentinels.

Work in a `tempfile.TemporaryDirectory()` so per-contig intermediates
(`.maf.gz`, `.ref`, `.bam`) never land in the output tree.

**Step 4: Run tests and commit**

```bash
cd $TOOLKIT && pixi run pytest tests/test_te_benchmark_panel.py -v
git add src/simulate_data/modules/te_benchmark_panel.py tests/test_te_benchmark_panel.py
git commit -m "feat(panel): simulate long-read component mixtures with gzip output"
```

---

## Task 10: End-to-end toolkit smoke test

**Files:**
- Create: `$TOOLKIT/tests/test_longread_e2e.py`

**Step 1: Write the test**

Mark it `@pytest.mark.slow`. Build a 100 kb single-contig fixture genome, a
2-event catalog, then run `_simulate_long_reads` for `ont-hq` at 2x and assert:

- `reads.fastq.gz` exists, is non-empty, and gunzips cleanly
- no `*.fastq` or `*.fq` anywhere under the output dir (compression contract)
- `component_manifest.tsv` has 8 rows
- `observed_support.tsv` has one row per truth event
- both `.complete` sentinels exist
- rerunning is a no-op that does not change the file mtime

Then repeat for `hifi` at 2x, asserting a `ccs_yield` column is populated.

**Step 2: Run**

Run: `cd $TOOLKIT && pixi run pytest tests/test_longread_e2e.py -v -m slow`
Expected: PASS (a few minutes)

**Step 3: Commit**

```bash
cd $TOOLKIT
git add tests/test_longread_e2e.py
git commit -m "test(longread): end-to-end panel simulation for both platforms"
```

---

## Task 11: Panel configuration

**Files:**
- Create: `$PANEL/config/riceTElib_lr_benchmark.toml`

```toml
# Long-read companion to riceTElib_benchmark.toml.
# The catalog is REUSED, not rebuilt: the truth set is bit-identical to the
# short-read panel so read technology is the only variable.
[inputs]
catalog = "results/riceTElib_benchmark"

[panel]
coverage = [5, 15, 30]
replicates = 3
seed = 916          # must match the short-read panel

[longread]
platforms = ["ont-hq", "hifi"]
min_ccs_yield = 0.80

[longread.ont-hq]
length_mean = 12000
length_sd = 9000
accuracy_mean = 0.99

[longread.hifi]
length_mean = 15000
length_sd = 2000
pass_num = 10
```

**Step 2: Verify it parses**

```bash
cd $PANEL && python3.12 -c "import tomllib;print(tomllib.load(open('config/riceTElib_lr_benchmark.toml','rb'))['longread']['platforms'])"
```
Expected: `['ont-hq', 'hifi']`

---

## Task 12: Panel driver

**Files:**
- Create: `$PANEL/pipeline/make_riceTElib_lr_benchmark/build_longread_panel.py`
- Test: `$PANEL/pipeline/make_riceTElib_lr_benchmark/test_build_longread_panel.py`

Model it on `$PANEL/pipeline/make_riceTElib_benchmark/build_multite_panel.py`.

**Critical constraint:** this script must NEVER write into the catalog
directory. It validates `<catalog>/.catalog.complete`, reads
`truth_events.tsv`, `component_coordinates.tsv`, and `genomes/`, and writes
only under its own `--output`.

CLI: `--config`, `--output`, `--platform`, `--coverage`, `--replicate`,
`--validate-only`.

`--validate-only` must confirm, without generating data: the catalog is
complete; all 8 component genomes plus `reference.fa` exist; both PBSIM3 models
resolve; and `ccs` is on PATH when `hifi` is among the platforms.

**Test** (run with the toolkit's pixi env):

```bash
cd $PANEL && pixi run --manifest-path $TOOLKIT/pyproject.toml \
  python -m unittest -v pipeline/make_riceTElib_lr_benchmark/test_build_longread_panel.py
```

Cover: refusing an incomplete catalog; refusing to write into the catalog dir;
task-index → (platform, coverage, replicate) mapping.

---

## Task 13: SLURM scripts

**Files:**
- Create: `$PANEL/pipeline/make_riceTElib_lr_benchmark/01_link_catalog.sh`
- Create: `$PANEL/pipeline/make_riceTElib_lr_benchmark/02_simulate_reads.sh`
- Create: `$PANEL/pipeline/make_riceTElib_lr_benchmark/submit.sh`

Follow `$PANEL/pipeline/make_riceTElib_benchmark/02_simulate_reads.sh` for
structure: `#!/usr/bin/bash -l`, `set -euo pipefail`, `THREADS` from
`SLURM_CPUS_PER_TASK`, date-stamped logging, explicit post-run validation.

**Array layout:** 18 tasks = 2 platforms x 3 coverages x 3 replicates, ordered
so `SLURM_ARRAY_TASK_ID` 0-8 is `ont-hq` and 9-17 is `hifi`.

**Resources.** HiFi 30x is the expensive corner — 10 passes over a 380 Mb
genome across 8 components, plus CPU-heavy ccs polishing. Do not use one
allocation for everything:

```bash
#SBATCH -p epyc
#SBATCH --mem=32gb
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH -o logs/riceTElib_lr_benchmark/reads.%A_%a.log
```

**Post-run validation**, mirroring the short-read script:

```bash
for required in \
    "$OUTPUT/reads/$PLATFORM/$SAMPLE/reads.fastq.gz" \
    "$OUTPUT/reads/$PLATFORM/${SAMPLE}_reference_control/reads.fastq.gz"; do
    [[ -s "$required" ]] || { echo "ERROR: missing $required" >&2; exit 1; }
done
if find "$OUTPUT/reads/$PLATFORM" -type f \( -name '*.fastq' -o -name '*.fq' \) \
    -print -quit | grep -q .; then
    echo "ERROR: uncompressed FASTQs remain" >&2
    exit 1
fi
```

`submit.sh` submits `01_link_catalog.sh`, then `02_simulate_reads.sh` as a
dependent array. Throttle with `%6` — 8 concurrent 380 Mb pbsim runs per task
is already heavy I/O.

**Do not submit anything as part of implementation.** Report the commands.

---

## Task 14: Documentation

**Files:**
- Create: `$PANEL/pipeline/make_riceTElib_lr_benchmark/README.md`
- Create: `$PANEL/docs/2026-08-07-longread-benchmark-workflow.md`

README mirrors the short-read one: purpose, design, commands, output contract,
rerun behavior, current status.

The dated doc follows the house format — date/time, purpose, current status,
commands, failures/issues, decisions/logic, next steps. Record: the two
upstream bugs found and fixed, the ccs yield caveat, the storage estimate
(~65 GB per technology), and that the catalog is shared read-only with the
short-read panel.

---

## Task 15: Final verification

**Step 1: Full toolkit test suite**

```bash
cd $TOOLKIT && pixi run pytest tests/ -v
```
Expected: all pass, no regressions in the ART short-read tests.

**Step 2: Lint**

```bash
cd $TOOLKIT && pixi run ruff check src/ tests/ && pixi run black --check src/ tests/
```

**Step 3: Panel dry run**

```bash
cd $PANEL && pixi run --manifest-path $TOOLKIT/pyproject.toml \
  python pipeline/make_riceTElib_lr_benchmark/build_longread_panel.py \
  --config config/riceTElib_lr_benchmark.toml \
  --output results/riceTElib_lr_benchmark \
  --validate-only
```
Expected: catalog complete, 9 genomes found, both models resolve, ccs present.

**Historical Step 4: Report, do not submit**

This was the original stop condition. The panel was later submitted; 15 tasks
completed and three HiFi 30x tasks timed out.
