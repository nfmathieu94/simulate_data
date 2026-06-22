"""Place structural variants (SVs) in a reference genome.

Wraps the SURVIVOR tool to introduce deletions, duplications,
inversions, and translocations into specified chromosomes.
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

VALID_SV_TYPES = {"DEL", "DUP", "INV", "TRA"}


def register_parser(parser):
    """Define the sv-placement subcommand's CLI arguments."""
    parser.add_argument(
        "--ref",
        required=True,
        help="Reference genome FASTA",
    )
    parser.add_argument(
        "--num-sv",
        type=int,
        required=True,
        help="Number of SVs to simulate",
    )
    parser.add_argument(
        "--sv-types",
        default="DEL,DUP,INV,TRA",
        help="Comma-separated SV types (DEL, DUP, INV, TRA). Default: all",
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
        "--config",
        default=None,
        help="SURVIVOR config file (if pre-generated)",
    )
    parser.add_argument(
        "--chroms",
        default="all",
        help=(
            "Chromosomes to place SVs in. "
            "Examples: 'Chr1', 'Chr2-5', 'Chr1,Chr3', 'all'"
        ),
    )


def _generate_survivor_config(
    config_path: Path,
    num_sv: int,
    sv_types: list[str],
    seed: int | None = None,
) -> Path:
    """Generate a SURVIVOR config file for simSV.

    The config file format for SURVIVOR simSV is a text file with
    one parameter per line.
    """
    lines = [
        "# SURVIVOR simSV configuration file",
        f"# SV types: {','.join(sv_types)}",
        f"# Number of SVs: {num_sv}",
        f"SV_TYPES={','.join(sv_types)}",
        f"NUM_SV={num_sv}",
        "SV_LENGTH_MIN=100",
        "SV_LENGTH_MAX=100000",
        "BALANCED=1",
    ]
    if seed is not None:
        lines.append(f"SEED={seed}")

    with open(config_path, "w") as f:
        f.write("\n".join(lines) + "\n")

    return config_path


def _build_survivor_command(
    ref_fasta: Path,
    config_file: Path,
    output_prefix: str,
) -> list[str]:
    """Build the SURVIVOR simSV command."""
    return [
        "SURVIVOR",
        "simSV",
        str(ref_fasta),
        str(config_file),
        output_prefix,
    ]


def main(args):
    """Place SVs in a reference genome using SURVIVOR."""
    # Validate inputs
    ref_path = Path(args.ref)
    validate_fasta(ref_path)

    if args.num_sv <= 0:
        raise ValueError(f"--num-sv must be positive, got {args.num_sv}")

    # Validate SV types
    sv_types = [t.strip().upper() for t in args.sv_types.split(",") if t.strip()]
    invalid_types = [t for t in sv_types if t not in VALID_SV_TYPES]
    if invalid_types:
        raise ValueError(
            f"Invalid SV types: {invalid_types}. Valid types: {sorted(VALID_SV_TYPES)}"
        )

    if args.config:
        validate_file_exists(Path(args.config))

    # Parse chromosome specification before checking tools
    available_chroms = get_chromosomes_from_fasta(ref_path)
    target_chroms = parse_chromosome_spec(args.chroms, available_chroms)

    logger.info(
        "Placing %d SVs (types: %s) in chromosomes: %s",
        args.num_sv,
        ",".join(sv_types),
        ", ".join(target_chroms),
    )

    # Check tool installation after all validation passes
    check_tool_installed("SURVIVOR")

    # Create output directory
    output_dir = ensure_output_dir(Path(args.output))

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)

        if args.chroms.strip().lower() == "all":
            ref_for_tool = ref_path
        else:
            ref_for_tool = tmpdir_path / "ref_subset.fa"
            extract_chromosomes(ref_path, target_chroms, ref_for_tool)
            logger.info(
                "Extracted %d chromosomes to %s",
                len(target_chroms),
                ref_for_tool,
            )

        # Generate or use config file
        if args.config:
            config_file = Path(args.config)
        else:
            config_file = tmpdir_path / "survivor_config.txt"
            _generate_survivor_config(
                config_path=config_file,
                num_sv=args.num_sv,
                sv_types=sv_types,
                seed=args.seed,
            )

        # Build and run SURVIVOR command
        output_prefix = str(output_dir / "sv_output")
        cmd = _build_survivor_command(
            ref_fasta=ref_for_tool,
            config_file=config_file,
            output_prefix=output_prefix,
        )

        run_command(cmd)

        # If we used a subset, merge the modified chromosomes back
        if args.chroms.strip().lower() != "all":
            modified_genome = output_dir / "sv_output_modified_genome.fa"
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

    logger.info("SV placement complete. Output: %s", output_dir)
