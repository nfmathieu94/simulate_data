"""Simulate Oxford Nanopore (ONT) long reads from a reference genome.

Wraps the PBSIM3 (pbsim) tool to generate realistic ONT reads
using quality score hidden Markov models (QSHMM).
"""

import logging
from pathlib import Path

from simulate_data.utils import (
    check_tool_installed,
    ensure_output_dir,
    run_command,
    validate_fasta,
)

logger = logging.getLogger(__name__)


def register_parser(parser):
    """Define the reads-ont subcommand's CLI arguments."""
    parser.add_argument(
        "--ref",
        required=True,
        help="Reference genome FASTA",
    )
    parser.add_argument(
        "--coverage",
        type=float,
        default=20.0,
        help="Coverage depth (default: 20)",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output directory for FASTQ files",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducibility",
    )
    parser.add_argument(
        "--read-length",
        type=int,
        default=9000,
        help="Mean read length in bp (default: 9000)",
    )
    parser.add_argument(
        "--read-std",
        type=int,
        default=7000,
        help="Read length standard deviation (default: 7000)",
    )
    parser.add_argument(
        "--error-model",
        default="QSHMM-ONT",
        help="PBSIM3 error model (default: QSHMM-ONT)",
    )
    parser.add_argument(
        "--qscore-model",
        default=None,
        help="PBSIM3 quality score model",
    )


def _build_pbsim3_ont_command(
    ref_fasta: Path,
    coverage: float,
    output_prefix: str,
    read_length: int = 9000,
    read_std: int = 7000,
    seed: int | None = None,
    error_model: str = "QSHMM-ONT",
    qscore_model: str | None = None,
) -> list[str]:
    """Build the PBSIM3 (pbsim) command for ONT read simulation."""
    cmd = [
        "pbsim",
        "--strategy",
        "wgs",
        "--method",
        "qshmm",
        "--qshmm",
        error_model,
        "--depth",
        str(coverage),
        "--length-mean",
        str(read_length),
        "--length-sd",
        str(read_std),
        "--genome",
        str(ref_fasta),
        "--prefix",
        output_prefix,
    ]

    if qscore_model is not None:
        cmd.extend(["--qscore-model", qscore_model])

    if seed is not None:
        cmd.extend(["--seed", str(seed)])

    return cmd


def main(args):
    """Simulate ONT reads from a reference genome using PBSIM3."""
    # Validate inputs
    ref_path = Path(args.ref)
    validate_fasta(ref_path)

    if args.coverage <= 0:
        raise ValueError(f"--coverage must be positive, got {args.coverage}")

    if args.read_length <= 0:
        raise ValueError(f"--read-length must be positive, got {args.read_length}")

    if args.read_std < 0:
        raise ValueError(f"--read-std must be non-negative, got {args.read_std}")

    # Check tool installation
    check_tool_installed("pbsim")

    # Create output directory
    output_dir = ensure_output_dir(Path(args.output))

    # Build and run PBSIM3 command
    output_prefix = str(output_dir / "ont_reads")
    cmd = _build_pbsim3_ont_command(
        ref_fasta=ref_path,
        coverage=args.coverage,
        output_prefix=output_prefix,
        read_length=args.read_length,
        read_std=args.read_std,
        seed=args.seed,
        error_model=args.error_model,
        qscore_model=args.qscore_model,
    )

    logger.info("Simulating ONT reads with PBSIM3")
    run_command(cmd)

    logger.info("ONT read simulation complete. Output: %s", output_dir)
