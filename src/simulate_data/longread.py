"""Shared helpers for long-read simulation with PBSIM3 and ccs.

This module lives at the package root rather than under ``modules/`` because
``cli.py`` auto-discovers every module in that package as a subcommand, and
these are helpers rather than a user-facing simulation module.
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
