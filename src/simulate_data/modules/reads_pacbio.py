"""Simulate PacBio long reads from a reference genome.

Wraps the PBSIM3 (pbsim) tool to generate realistic PacBio
continuous long reads (CLR) or HiFi circular consensus reads.
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
    """Define the reads-pacbio subcommand's CLI arguments."""
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
        help="Output directory for FASTQ/BAM files",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducibility",
    )
    parser.add_argument(
        "--read-type",
        choices=["CLR", "HiFi"],
        default="CLR",
        help="Read type: CLR or HiFi (default: CLR)",
    )
    parser.add_argument(
        "--read-length",
        type=int,
        default=15000,
        help="Mean read length in bp (default: 15000)",
    )
    parser.add_argument(
        "--read-std",
        type=int,
        default=13000,
        help="Read length standard deviation (default: 13000)",
    )
    parser.add_argument(
        "--pass-num",
        type=int,
        default=10,
        help="Number of passes for HiFi (default: 10)",
    )
    parser.add_argument(
        "--error-model",
        default="QSHMM-RSII",
        help="PBSIM3 error model (default: QSHMM-RSII)",
    )
    parser.add_argument(
        "--qscore-model",
        default=None,
        help="PBSIM3 quality score model",
    )


def _build_pbsim3_pacbio_command(
    ref_fasta: Path,
    coverage: float,
    output_prefix: str,
    read_type: str = "CLR",
    read_length: int = 15000,
    read_std: int = 13000,
    pass_num: int = 10,
    seed: int | None = None,
    error_model: str = "QSHMM-RSII",
    qscore_model: str | None = None,
) -> list[str]:
    """Build the PBSIM3 (pbsim) command for PacBio read simulation."""
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

    if read_type == "HiFi":
        cmd.extend(["--pass-num", str(pass_num)])

    if qscore_model is not None:
        cmd.extend(["--qscore-model", qscore_model])

    if seed is not None:
        cmd.extend(["--seed", str(seed)])

    return cmd


def main(args):
    """Simulate PacBio reads from a reference genome using PBSIM3."""
    # Validate inputs
    ref_path = Path(args.ref)
    validate_fasta(ref_path)

    if args.coverage <= 0:
        raise ValueError(f"--coverage must be positive, got {args.coverage}")

    if args.read_length <= 0:
        raise ValueError(f"--read-length must be positive, got {args.read_length}")

    if args.read_std < 0:
        raise ValueError(f"--read-std must be non-negative, got {args.read_std}")

    if args.read_type == "HiFi" and args.pass_num <= 0:
        raise ValueError(f"--pass-num must be positive for HiFi, got {args.pass_num}")

    # Check tool installation
    check_tool_installed("pbsim")

    # Create output directory
    output_dir = ensure_output_dir(Path(args.output))

    # Build and run PBSIM3 command
    output_prefix = str(output_dir / "pacbio_reads")
    cmd = _build_pbsim3_pacbio_command(
        ref_fasta=ref_path,
        coverage=args.coverage,
        output_prefix=output_prefix,
        read_type=args.read_type,
        read_length=args.read_length,
        read_std=args.read_std,
        pass_num=args.pass_num,
        seed=args.seed,
        error_model=args.error_model,
        qscore_model=args.qscore_model,
    )

    logger.info("Simulating PacBio %s reads with PBSIM3", args.read_type)
    run_command(cmd)

    logger.info("PacBio read simulation complete. Output: %s", output_dir)
