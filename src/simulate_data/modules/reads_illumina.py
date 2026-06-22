"""Simulate Illumina short reads from a reference genome.

Wraps the ART (art_illumina) tool to generate realistic Illumina
paired-end or single-end reads with empirical error profiles.
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
    """Define the reads-illumina subcommand's CLI arguments."""
    parser.add_argument(
        "--ref",
        required=True,
        help="Reference genome FASTA",
    )
    parser.add_argument(
        "--read-length",
        type=int,
        default=150,
        help="Read length in bp (default: 150)",
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
        "--paired",
        action="store_true",
        default=True,
        help="Generate paired-end reads (default: True)",
    )
    parser.add_argument(
        "--fragment-size",
        type=int,
        default=300,
        help="Mean fragment size for paired-end (default: 300)",
    )
    parser.add_argument(
        "--fragment-std",
        type=int,
        default=30,
        help="Fragment size standard deviation (default: 30)",
    )
    parser.add_argument(
        "--profile",
        default=None,
        help="ART error profile path (e.g., HiSeq, NovaSeq)",
    )


def _build_art_command(
    ref_fasta: Path,
    read_length: int,
    coverage: float,
    output_prefix: str,
    paired: bool = True,
    fragment_size: int = 300,
    fragment_std: int = 30,
    seed: int | None = None,
    profile: str | None = None,
) -> list[str]:
    """Build the art_illumina command."""
    cmd = [
        "art_illumina",
        "-i",
        str(ref_fasta),
        "-l",
        str(read_length),
        "-f",
        str(coverage),
        "-o",
        output_prefix,
    ]

    if paired:
        cmd.extend(["--paired", "-m", str(fragment_size), "-s", str(fragment_std)])

    if seed is not None:
        cmd.extend(["-rs", str(seed)])

    if profile is not None:
        cmd.extend(["-sp", profile])

    return cmd


def main(args):
    """Simulate Illumina reads from a reference genome using ART."""
    # Validate inputs
    ref_path = Path(args.ref)
    validate_fasta(ref_path)

    if args.read_length <= 0:
        raise ValueError(f"--read-length must be positive, got {args.read_length}")

    if args.coverage <= 0:
        raise ValueError(f"--coverage must be positive, got {args.coverage}")

    if args.fragment_size <= 0:
        raise ValueError(f"--fragment-size must be positive, got {args.fragment_size}")

    if args.fragment_std < 0:
        raise ValueError(
            f"--fragment-std must be non-negative, got {args.fragment_std}"
        )

    # Check tool installation
    check_tool_installed("art_illumina")

    # Create output directory
    output_dir = ensure_output_dir(Path(args.output))

    # Build and run ART command
    output_prefix = str(output_dir / "illumina_reads")
    cmd = _build_art_command(
        ref_fasta=ref_path,
        read_length=args.read_length,
        coverage=args.coverage,
        output_prefix=output_prefix,
        paired=args.paired,
        fragment_size=args.fragment_size,
        fragment_std=args.fragment_std,
        seed=args.seed,
        profile=args.profile,
    )

    logger.info("Simulating Illumina reads with ART")
    run_command(cmd)

    logger.info("Illumina read simulation complete. Output: %s", output_dir)
