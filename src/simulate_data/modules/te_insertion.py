"""Insert transposable elements (TEs) into a reference genome.

Wraps the TEvarSim tool to insert TEs from a TE consensus FASTA
into specified chromosomes of a reference genome.

The TEvarSim workflow has two steps:
1. TErandom: Generates a TE pool FASTA and BED file with random TE
   positions from a TE consensus FASTA and a known TE deletion file.
2. Simulate: Uses the TE pool and BED file to create a modified
   genome FASTA and VCF file.
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
        "--known-del",
        required=True,
        help="Known TE deletion file (RepeatMasker .out or UCSC .txt)",
    )
    parser.add_argument(
        "--num",
        type=int,
        required=True,
        help="Number of TE events to simulate (insertions + deletions)",
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
        help="Pre-generated BED file of TE positions (skips TErandom step)",
    )
    parser.add_argument(
        "--chroms",
        default="all",
        help=(
            "Chromosomes to insert TEs into. "
            "Examples: 'Chr1', 'Chr2-5', 'Chr1,Chr3', 'all'"
        ),
    )
    parser.add_argument(
        "--num-genomes",
        type=int,
        default=1,
        help="Number of genomes to simulate (default: 1)",
    )
    parser.add_argument(
        "--ins-ratio",
        type=float,
        default=0.6,
        help=(
            "Proportion of insertion events among all simulated pTE (0-1, default: 0.6)"
        ),
    )


def _build_terandom_command(
    te_fasta: Path,
    known_del: Path,
    chrom: str,
    num_te: int,
    outprefix: str,
    seed: int | None = None,
    ins_ratio: float = 0.6,
) -> list[str]:
    """Build the TEvarSim TErandom command.

    TErandom generates a TE pool FASTA (``{outprefix}.fa``) and a BED
    file with random TE insertion/deletion positions
    (``{outprefix}.bed``).
    """
    cmd = [
        "tevarsim",
        "TErandom",
        "--consensus",
        str(te_fasta),
        "--knownDEL",
        str(known_del),
        "--CHR",
        chrom,
        "--nTE",
        str(num_te),
        "--outprefix",
        outprefix,
        "--ins-ratio",
        str(ins_ratio),
    ]
    if seed is not None:
        cmd.extend(["--seed", str(seed)])
    return cmd


def _build_simulate_command(
    ref_fasta: Path,
    te_pool: Path,
    bed_file: Path,
    num_genomes: int,
    outprefix: str,
    seed: int | None = None,
) -> list[str]:
    """Build the TEvarSim Simulate command.

    Simulate reads the TE pool and BED file produced by TErandom and
    writes a modified genome FASTA (``{outprefix}.fa``) and a VCF file
    (``{outprefix}.vcf``).
    """
    cmd = [
        "tevarsim",
        "Simulate",
        "--ref",
        str(ref_fasta),
        "--pool",
        str(te_pool),
        "--bed",
        str(bed_file),
        "--num",
        str(num_genomes),
        "--outprefix",
        outprefix,
    ]
    if seed is not None:
        cmd.extend(["--seed", str(seed)])
    return cmd


def main(args):
    """Insert TEs into a reference genome using TEvarSim.

    Workflow per target chromosome:
      1. Extract the chromosome into a single-sequence FASTA.
      2. Run ``tevarsim TErandom`` to build a TE pool and BED file.
      3. Run ``tevarsim Simulate`` to produce the modified chromosome
         FASTA and VCF.

    After processing all chromosomes, modified chromosomes are merged
    back with the unmodified chromosomes from the original reference.
    """
    ref_path = Path(args.ref)
    te_path = Path(args.te)
    known_del_path = Path(args.known_del)
    validate_fasta(ref_path)
    validate_fasta(te_path)
    validate_file_exists(known_del_path)

    if args.num <= 0:
        raise ValueError(f"--num must be positive, got {args.num}")

    if args.bed:
        validate_file_exists(Path(args.bed))

    check_tool_installed("tevarsim")

    available_chroms = get_chromosomes_from_fasta(ref_path)
    target_chroms = parse_chromosome_spec(args.chroms, available_chroms)

    logger.info(
        "Inserting %d TEs into chromosomes: %s",
        args.num,
        ", ".join(target_chroms),
    )

    output_dir = ensure_output_dir(Path(args.output))

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)

        if args.bed:
            cmd = _build_simulate_command(
                ref_fasta=ref_path,
                te_pool=te_path,
                bed_file=Path(args.bed),
                num_genomes=args.num_genomes,
                outprefix=str(output_dir / "Sim"),
                seed=args.seed,
            )
            run_command(cmd)
        else:
            for chrom in target_chroms:
                logger.info("Processing chromosome %s", chrom)

                chrom_fasta = tmpdir_path / f"{chrom}.fa"
                extract_chromosomes(ref_path, [chrom], chrom_fasta)

                terandom_prefix = str(tmpdir_path / f"terandom_{chrom}")
                terandom_cmd = _build_terandom_command(
                    te_fasta=te_path,
                    known_del=known_del_path,
                    chrom=chrom,
                    num_te=args.num,
                    outprefix=terandom_prefix,
                    seed=args.seed,
                    ins_ratio=args.ins_ratio,
                )
                run_command(terandom_cmd)

                te_pool = Path(f"{terandom_prefix}.fa")
                bed_file = Path(f"{terandom_prefix}.bed")
                sim_prefix = str(output_dir / f"Sim_{chrom}")
                simulate_cmd = _build_simulate_command(
                    ref_fasta=chrom_fasta,
                    te_pool=te_pool,
                    bed_file=bed_file,
                    num_genomes=args.num_genomes,
                    outprefix=sim_prefix,
                    seed=args.seed,
                )
                run_command(simulate_cmd)

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
