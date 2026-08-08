"""Shared helpers for long-read simulation with PBSIM3 and ccs.

This module lives at the package root rather than under ``modules/`` because
``cli.py`` auto-discovers every module in that package as a subcommand, and
these are helpers rather than a user-facing simulation module.
"""

from __future__ import annotations

import os
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
