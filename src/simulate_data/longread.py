"""Shared helpers for long-read simulation with PBSIM3 and ccs.

This module lives at the package root rather than under ``modules/`` because
``cli.py`` auto-discovers every module in that package as a subcommand, and
these are helpers rather than a user-facing simulation module.
"""

from __future__ import annotations

import gzip
import os
from collections.abc import Iterator
from dataclasses import dataclass
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

    PBSIM3 requires a path to a ``.model`` file; passing a bare name such as
    ``QSHMM-ONT`` fails at runtime with ``Cannot open file``.
    """
    candidate = Path(model)
    if candidate.is_file():
        return candidate

    directory = _model_dir()
    name = model if model.endswith(".model") else f"{model}.model"
    resolved = directory / name
    if resolved.is_file():
        return resolved

    available = sorted(path.stem for path in directory.glob("*.model"))
    raise FileNotFoundError(
        f"PBSIM3 model not found: {model}. Looked for {resolved}. "
        f"Available models: {', '.join(available) or '(none)'}"
    )


@dataclass(frozen=True)
class Platform:
    """A long-read sequencing platform preset."""

    key: str
    method: str  # "qshmm" or "errhmm"
    model: str
    difference_ratio: str
    length_mean: int
    length_sd: int
    accuracy_mean: float
    pass_num: int  # 1 = single pass; >1 produces subreads requiring ccs
    needs_ccs: bool


# Difference ratios are PBSIM3's documented per-platform values. The tool's
# own default of 6:55:39 describes PacBio RS II and is wrong for both
# platforms used here.
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
    """Build the PBSIM3 command simulating one component genome."""
    model_path = resolve_model_path(platform.model)
    cmd = [
        "pbsim",
        "--strategy",
        "wgs",
        "--method",
        platform.method,
        f"--{platform.method}",
        str(model_path),
        "--genome",
        str(genome),
        # Fixed-point formatting: PBSIM3 cannot parse scientific notation,
        # which float formatting produces for very small component depths.
        "--depth",
        f"{depth:.6f}",
        "--length-mean",
        str(platform.length_mean),
        "--length-sd",
        str(platform.length_sd),
        "--accuracy-mean",
        str(platform.accuracy_mean),
        "--difference-ratio",
        platform.difference_ratio,
        "--prefix",
        str(prefix),
        "--seed",
        str(seed),
    ]
    if platform.pass_num > 1:
        cmd.extend(["--pass-num", str(platform.pass_num)])
    return cmd


def contig_name_for(ref_path: Path) -> str:
    """Return the contig name recorded in a PBSIM3 ``.ref`` file.

    PBSIM3 names the MAF source ``ref`` rather than the contig, so contig
    identity has to be recovered from the sibling ``.ref`` file's header.
    """
    with ref_path.open() as handle:
        header = handle.readline().strip()
    if not header.startswith(">"):
        raise ValueError(f"Not a FASTA header in {ref_path}: {header!r}")
    return header[1:].split()[0]


def iter_maf_alignments(maf_path: Path) -> Iterator[tuple[int, int]]:
    """Yield ``(start, end)`` reference spans from a PBSIM3 MAF file.

    Fields are ``s <src> <start> <size> <strand> <srcSize> <seq>`` with a
    0-based start and a size counting only non-gap characters, so the span is
    the direct analogue of pysam's ``reference_start``/``reference_end``.
    Only the first ``s`` line of each block -- the reference -- is used.
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


# Matches PANEL/pipeline/fastq_compression.py so both panels compress alike.
GZIP_LEVEL = 6


def pool_fastq_gz(source: Path, destination, prefix: str) -> int:
    """Append a gzipped FASTQ into an open gzip stream, prefixing read names.

    Read names collide across components -- every PBSIM3 run emits ``S1_1``,
    ``S1_2``, ... -- so the prefix is what makes pooled names unique. Nothing
    is decompressed to disk.
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
