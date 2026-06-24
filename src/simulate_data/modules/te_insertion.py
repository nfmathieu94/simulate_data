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
        help=(
            "Known TE deletion file. Supports TEvarSim-compatible "
            "RepeatMasker .out/UCSC .txt, or RepeatMasker GFF3."
        ),
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
    parser.add_argument(
        "--te-type",
        action="append",
        default=None,
        help=(
            "TE family/superfamily to extract from the known deletion file. "
            "Can be repeated or comma-separated. Example: --te-type Harbinger"
        ),
    )
    parser.add_argument(
        "--snp-rate",
        type=float,
        default=0.02,
        help="SNP mutation rate per base for inserted TE copies (default: 0.02)",
    )
    parser.add_argument(
        "--indel-rate",
        type=float,
        default=0.005,
        help="Indel mutation rate per base for inserted TE copies (default: 0.005)",
    )
    parser.add_argument(
        "--indel-ins",
        type=float,
        default=0.4,
        help="Proportion of TE-copy indels that are insertions (0-1, default: 0.4)",
    )
    parser.add_argument(
        "--indel-geom-p",
        type=float,
        default=0.7,
        help="Geometric distribution p for TE-copy indel lengths (default: 0.7)",
    )
    parser.add_argument(
        "--truncated-ratio",
        type=float,
        default=0.3,
        help="Proportion of inserted TE copies to truncate (0-1, default: 0.3)",
    )
    parser.add_argument(
        "--truncated-max-length",
        type=float,
        default=0.5,
        help=(
            "Maximum proportion of each TE copy that can be truncated "
            "(0-1, default: 0.5)"
        ),
    )
    parser.add_argument(
        "--polyA-ratio",
        type=float,
        default=0.8,
        help="Proportion of inserted TE copies to add a polyA tail (0-1, default: 0.8)",
    )
    parser.add_argument(
        "--polyA-min",
        type=int,
        default=5,
        help="Minimum polyA tail length for inserted TE copies (default: 5)",
    )
    parser.add_argument(
        "--polyA-max",
        type=int,
        default=20,
        help="Maximum polyA tail length for inserted TE copies (default: 20)",
    )
    parser.add_argument(
        "--sense-strand-ratio",
        type=float,
        default=0.5,
        help=(
            "Proportion of TE insertions simulated on the sense strand "
            "(0-1, default: 0.5)"
        ),
    )


def _split_te_types(te_types: list[str] | None) -> list[str]:
    """Normalize repeated or comma-separated --te-type values."""
    if not te_types:
        return []

    result: list[str] = []
    for value in te_types:
        for item in value.split(","):
            item = item.strip()
            if item:
                result.append(item)
    return result


def _parse_gff_attributes(attributes: str) -> dict[str, str]:
    """Parse a simple GFF3 attributes field."""
    parsed: dict[str, str] = {}
    for item in attributes.rstrip(";").split(";"):
        if not item or "=" not in item:
            continue
        key, value = item.split("=", 1)
        parsed[key] = value
    return parsed


def _repeatmasker_gff_to_out(gff_path: Path, output_path: Path) -> Path:
    """Convert RepeatMasker GFF3 to the .out-like columns TEvarSim reads."""
    converted = 0

    with open(gff_path, "r") as fin, open(output_path, "w") as fout:
        for line in fin:
            if line.startswith("#") or not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue

            chrom, _source, _feature, start, end, score, strand, _phase, attrs = fields
            attr = _parse_gff_attributes(attrs)
            repeat_class = attr.get("Class")
            target = attr.get("Target")

            if not repeat_class or "/" not in repeat_class or not target:
                continue

            repeat_name = target.split()[0]
            rm_strand = "C" if strand == "-" else "+"
            repeat_end = max(int(end) - int(start) + 1, 1)

            fout.write(
                " ".join(
                    [
                        score if score != "." else "0",
                        attr.get("PercDiv", "0"),
                        attr.get("PercDel", "0"),
                        attr.get("PercIns", "0"),
                        chrom,
                        start,
                        end,
                        "(0)",
                        rm_strand,
                        repeat_name,
                        repeat_class,
                        "1",
                        str(repeat_end),
                        "(0)",
                        attr.get("ID", str(converted + 1)),
                    ]
                )
                + "\n"
            )
            converted += 1

    if converted == 0:
        raise ValueError(f"No RepeatMasker records could be converted from {gff_path}")

    logger.info(
        "Converted %d RepeatMasker GFF records to TEvarSim .out format: %s",
        converted,
        output_path,
    )
    return output_path


def _prepare_known_del_file(known_del: Path, tmpdir: Path) -> Path:
    """Return a TEvarSim-compatible known deletion file."""
    if known_del.suffix.lower() in {".out", ".txt"}:
        return known_del
    if known_del.suffix.lower() in {".gff", ".gff3"}:
        return _repeatmasker_gff_to_out(known_del, tmpdir / f"{known_del.stem}.out")

    raise ValueError(
        "Known deletion file must end with .out, .txt, .gff, or .gff3; "
        f"got {known_del}"
    )


def _build_terandom_command(
    te_fasta: Path,
    known_del: Path,
    chrom: str,
    num_te: int,
    outprefix: str,
    seed: int | None = None,
    ins_ratio: float = 0.6,
    te_types: list[str] | None = None,
    snp_rate: float = 0.02,
    indel_rate: float = 0.005,
    indel_ins: float = 0.4,
    indel_geom_p: float = 0.7,
    truncated_ratio: float = 0.3,
    truncated_max_length: float = 0.5,
    polya_ratio: float = 0.8,
    polya_min: int = 5,
    polya_max: int = 20,
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
        "--snp-rate",
        str(snp_rate),
        "--indel-rate",
        str(indel_rate),
        "--indel-ins",
        str(indel_ins),
        "--indel-geom-p",
        str(indel_geom_p),
        "--truncated-ratio",
        str(truncated_ratio),
        "--truncated-max-length",
        str(truncated_max_length),
        "--polyA-ratio",
        str(polya_ratio),
        "--polyA-min",
        str(polya_min),
        "--polyA-max",
        str(polya_max),
    ]
    if seed is not None:
        cmd.extend(["--seed", str(seed)])
    for te_type in te_types or []:
        cmd.extend(["--TEtype", te_type])
    return cmd


def _build_simulate_command(
    ref_fasta: Path,
    te_pool: Path,
    bed_file: Path,
    num_genomes: int,
    outprefix: str,
    seed: int | None = None,
    sense_strand_ratio: float = 0.5,
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
        "--sense-strand-ratio",
        str(sense_strand_ratio),
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

    _validate_tevarsim_rates(args)

    if args.bed:
        validate_file_exists(Path(args.bed))

    available_chroms = get_chromosomes_from_fasta(ref_path)
    target_chroms = parse_chromosome_spec(args.chroms, available_chroms)

    check_tool_installed("tevarsim")

    logger.info(
        "Inserting %d TEs into chromosomes: %s",
        args.num,
        ", ".join(target_chroms),
    )

    output_dir = ensure_output_dir(Path(args.output))

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        known_del_for_tool = _prepare_known_del_file(known_del_path, tmpdir_path)
        te_types = _split_te_types(getattr(args, "te_type", None))

        if args.bed:
            cmd = _build_simulate_command(
                ref_fasta=ref_path,
                te_pool=te_path,
                bed_file=Path(args.bed),
                num_genomes=args.num_genomes,
                outprefix=str(output_dir / "Sim"),
                seed=args.seed,
                sense_strand_ratio=getattr(args, "sense_strand_ratio", 0.5),
            )
            run_command(cmd)
        else:
            if args.num_genomes != 1:
                raise ValueError(
                    "--num-genomes > 1 is only supported without chromosome "
                    "subsetting in the current TE wrapper"
                )

            modified_subset_fasta = tmpdir_path / "modified_chromosomes.fa"
            for chrom in target_chroms:
                logger.info("Processing chromosome %s", chrom)

                chrom_fasta = tmpdir_path / f"{chrom}.fa"
                extract_chromosomes(ref_path, [chrom], chrom_fasta)

                terandom_prefix = str(tmpdir_path / f"terandom_{chrom}")
                terandom_cmd = _build_terandom_command(
                    te_fasta=te_path,
                    known_del=known_del_for_tool,
                    chrom=chrom,
                    num_te=args.num,
                    outprefix=terandom_prefix,
                    seed=args.seed,
                    ins_ratio=args.ins_ratio,
                    te_types=te_types,
                    snp_rate=getattr(args, "snp_rate", 0.02),
                    indel_rate=getattr(args, "indel_rate", 0.005),
                    indel_ins=getattr(args, "indel_ins", 0.4),
                    indel_geom_p=getattr(args, "indel_geom_p", 0.7),
                    truncated_ratio=getattr(args, "truncated_ratio", 0.3),
                    truncated_max_length=getattr(args, "truncated_max_length", 0.5),
                    polya_ratio=getattr(args, "polyA_ratio", 0.8),
                    polya_min=getattr(args, "polyA_min", 5),
                    polya_max=getattr(args, "polyA_max", 20),
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
                    sense_strand_ratio=getattr(args, "sense_strand_ratio", 0.5),
                )
                run_command(simulate_cmd)

                sim_fasta = Path(f"{sim_prefix}.fa")
                if not sim_fasta.exists():
                    raise FileNotFoundError(
                        f"Expected TEvarSim output FASTA was not created: {sim_fasta}"
                    )
                _append_renamed_first_fasta_record(
                    input_fasta=sim_fasta,
                    output_fasta=modified_subset_fasta,
                    record_name=chrom,
                )

            merged_output = output_dir / "final_genome.fa"

            merge_fasta(
                base_fasta=ref_path,
                modified_fasta=modified_subset_fasta,
                output_path=merged_output,
                modified_chroms=set(target_chroms),
            )
            logger.info(
                "Merged modified chromosomes with reference. Final genome: %s",
                merged_output,
            )

    logger.info("TE insertion complete. Output: %s", output_dir)


def _validate_tevarsim_rates(args) -> None:
    """Validate TEvarSim rate-style arguments before launching jobs."""
    rate_fields = [
        "ins_ratio",
        "snp_rate",
        "indel_rate",
        "indel_ins",
        "indel_geom_p",
        "truncated_ratio",
        "truncated_max_length",
        "polyA_ratio",
        "sense_strand_ratio",
    ]
    defaults = {
        "snp_rate": 0.02,
        "indel_rate": 0.005,
        "indel_ins": 0.4,
        "indel_geom_p": 0.7,
        "truncated_ratio": 0.3,
        "truncated_max_length": 0.5,
        "polyA_ratio": 0.8,
        "sense_strand_ratio": 0.5,
    }
    for field in rate_fields:
        value = getattr(args, field, defaults.get(field))
        if value is None:
            continue
        if not 0 <= value <= 1:
            raise ValueError(f"--{field.replace('_', '-')} must be between 0 and 1")

    polya_min = getattr(args, "polyA_min", 5)
    polya_max = getattr(args, "polyA_max", 20)
    if polya_min < 0:
        raise ValueError("--polyA-min must be non-negative")
    if polya_max < polya_min:
        raise ValueError("--polyA-max must be greater than or equal to --polyA-min")


def _append_renamed_first_fasta_record(
    input_fasta: Path,
    output_fasta: Path,
    record_name: str,
) -> None:
    """Append the first FASTA record from input_fasta using record_name."""
    sequence_lines: list[str] = []
    seen_header = False

    with open(input_fasta, "r") as fin:
        for line in fin:
            if line.startswith(">"):
                if seen_header:
                    break
                seen_header = True
                continue
            if seen_header:
                sequence_lines.append(line)

    if not seen_header or not sequence_lines:
        raise ValueError(f"No FASTA record found in TEvarSim output: {input_fasta}")

    with open(output_fasta, "a") as fout:
        fout.write(f">{record_name}\n")
        fout.writelines(sequence_lines)
