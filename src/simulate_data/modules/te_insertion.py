"""Insert transposable elements (TEs) into a reference genome.

Wraps the TEvarSim tool to insert TEs from a TE consensus FASTA
into specified chromosomes of a reference genome.
"""

import logging
import tempfile
from pathlib import Path

from simulate_data.utils import (
    check_tool_installed,
    ensure_output_dir,
    extract_chromosomes,
    get_chromosomes_from_fasta,
    merge_fasta,
    parse_chromosome_spec,
    run_command,
    validate_fasta,
    validate_file_exists,
)

logger = logging.getLogger(__name__)


def register_parser(parser):
    """Define the te-insertion subcommand's CLI arguments."""
    parser.add_argument(
        "--ref",
        required=True,
        help="Reference genome FASTA",
    )
    parser.add_argument(
        "--te",
        required=True,
        help="TE consensus FASTA (TE pool)",
    )
    parser.add_argument(
        "--num",
        type=int,
        required=True,
        help="Number of TE insertions to simulate",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output directory for modified genome and VCF",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducibility",
    )
    parser.add_argument(
        "--bed",
        default=None,
        help="BED file of TE positions (if pre-generated)",
    )
    parser.add_argument(
        "--chroms",
        default="all",
        help=(
            "Chromosomes to insert TEs into. "
            "Examples: 'Chr1', 'Chr2-5', 'Chr1,Chr3', 'all'"
        ),
    )


def _build_tevarsim_command(
    ref_fasta: Path,
    te_fasta: Path,
    num_insertions: int,
    output_dir: Path,
    seed: int | None = None,
    bed: Path | None = None,
) -> list[str]:
    """Build the TEvarSim Simulate command."""
    cmd = [
        "tevarsim",
        "Simulate",
        "--ref",
        str(ref_fasta),
        "--te-pool",
        str(te_fasta),
        "--num",
        str(num_insertions),
        "--output",
        str(output_dir),
    ]
    if seed is not None:
        cmd.extend(["--seed", str(seed)])
    if bed is not None:
        cmd.extend(["--bed", str(bed)])
    return cmd


def main(args):
    """Insert TEs into a reference genome using TEvarSim."""
    # Validate inputs
    ref_path = Path(args.ref)
    te_path = Path(args.te)
    validate_fasta(ref_path)
    validate_fasta(te_path)

    if args.num <= 0:
        raise ValueError(f"--num must be positive, got {args.num}")

    if args.bed:
        validate_file_exists(Path(args.bed))

    # Check tool installation
    check_tool_installed("tevarsim")

    # Parse chromosome specification
    available_chroms = get_chromosomes_from_fasta(ref_path)
    target_chroms = parse_chromosome_spec(args.chroms, available_chroms)

    logger.info(
        "Inserting %d TEs into chromosomes: %s",
        args.num,
        ", ".join(target_chroms),
    )

    # Create output directory
    output_dir = ensure_output_dir(Path(args.output))

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)

        if args.chroms.strip().lower() == "all":
            # Use full reference genome
            ref_for_tool = ref_path
        else:
            # Extract target chromosomes into a temporary FASTA
            ref_for_tool = tmpdir_path / "ref_subset.fa"
            extract_chromosomes(ref_path, target_chroms, ref_for_tool)
            logger.info(
                "Extracted %d chromosomes to %s",
                len(target_chroms),
                ref_for_tool,
            )

        # Build and run TEvarSim command
        cmd = _build_tevarsim_command(
            ref_fasta=ref_for_tool,
            te_fasta=te_path,
            num_insertions=args.num,
            output_dir=output_dir,
            seed=args.seed,
            bed=Path(args.bed) if args.bed else None,
        )

        run_command(cmd)

        # If we used a subset, merge the modified chromosomes back
        # with the unmodified chromosomes from the original reference
        if args.chroms.strip().lower() != "all":
            modified_genome = output_dir / "modified_genome.fa"
            merged_output = output_dir / "final_genome.fa"

            if modified_genome.exists():
                merge_fasta(
                    base_fasta=ref_path,
                    modified_fasta=modified_genome,
                    output_path=merged_output,
                    modified_chroms=set(target_chroms),
                )
                logger.info(
                    "Merged modified chromosomes with reference. Final genome: %s",
                    merged_output,
                )

    logger.info("TE insertion complete. Output: %s", output_dir)
