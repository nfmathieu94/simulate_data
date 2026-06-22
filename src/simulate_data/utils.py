"""Shared utilities for simulation modules."""

import logging
import re
import shutil
import subprocess
import sys
from pathlib import Path

logger = logging.getLogger(__name__)


def setup_logging(verbose: bool = False) -> None:
    """Configure logging for the application."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def run_command(cmd: list[str], check: bool = True) -> subprocess.CompletedProcess:
    """Run a command and return the result.

    Args:
        cmd: Command as a list of strings.
        check: If True, raise RuntimeError on non-zero exit.

    Returns:
        CompletedProcess instance.
    """
    logger.info("Running command: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if check and result.returncode != 0:
        logger.error("Command failed (exit %d): %s", result.returncode, " ".join(cmd))
        if result.stderr:
            logger.error("stderr: %s", result.stderr)
        raise RuntimeError(
            f"Command failed with exit code {result.returncode}:\n"
            f"  Command: {' '.join(cmd)}\n"
            f"  stderr: {result.stderr}"
        )
    return result


def validate_file_exists(path: Path) -> None:
    """Check that a file exists and is not a directory.

    Args:
        path: Path to the file.

    Raises:
        FileNotFoundError: If the file does not exist.
    """
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")
    if not path.is_file():
        raise FileNotFoundError(f"Path is not a file: {path}")


def validate_fasta(path: Path) -> None:
    """Validate that a file is a proper FASTA file.

    Checks that the file exists, starts with a header line (">"),
    and contains at least one sequence.

    Args:
        path: Path to the FASTA file.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the file is not a valid FASTA.
    """
    validate_file_exists(path)

    with open(path, "r") as f:
        first_line = f.readline()
        if not first_line.startswith(">"):
            raise ValueError(f"FASTA file must start with '>' header: {path}")

        has_sequence = False
        for line in f:
            stripped = line.strip()
            if stripped and not stripped.startswith(">"):
                has_sequence = True
                break

        if not has_sequence:
            raise ValueError(f"FASTA file has no sequence data: {path}")


def check_tool_installed(tool_name: str) -> str:
    """Check that an external tool is installed and on PATH.

    Args:
        tool_name: Name of the executable to find.

    Returns:
        Path to the tool executable.

    Raises:
        FileNotFoundError: If the tool is not found on PATH.
    """
    tool_path = shutil.which(tool_name)
    if tool_path is None:
        raise FileNotFoundError(
            f"Tool '{tool_name}' not found on PATH. "
            f"Please install it via pixi/conda or add it to your PATH."
        )
    return tool_path


def ensure_output_dir(path: Path) -> Path:
    """Create output directory if it does not exist.

    Args:
        path: Path to the output directory.

    Returns:
        The resolved Path object.
    """
    path.mkdir(parents=True, exist_ok=True)
    return path


def get_chromosomes_from_fasta(fasta_path: Path) -> list[str]:
    """Extract chromosome names from a FASTA file.

    Reads only the header lines (lines starting with ">") to get
    chromosome names efficiently without loading sequences.

    Args:
        fasta_path: Path to the FASTA file.

    Returns:
        List of chromosome name strings in file order.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the file is empty or has no headers.
    """
    validate_file_exists(fasta_path)

    chrom_names: list[str] = []
    with open(fasta_path, "r") as f:
        for line in f:
            if line.startswith(">"):
                # Extract the chromosome name (first token after ">")
                name = line[1:].strip().split()[0]
                chrom_names.append(name)

    if not chrom_names:
        raise ValueError(f"No chromosome headers found in FASTA: {fasta_path}")

    return chrom_names


def parse_chromosome_spec(
    spec: str,
    available_chroms: list[str],
) -> list[str]:
    """Parse a chromosome specification string into a list of chromosome names.

    Supports the following formats:
        "all"             -> all available chromosomes
        "Chr1"            -> ["Chr1"]
        "Chr1,Chr3,Chr5"  -> ["Chr1", "Chr3", "Chr5"]
        "Chr2-5"          -> ["Chr2", "Chr3", "Chr4", "Chr5"]
        "Chr2-Chr5"       -> ["Chr2", "Chr3", "Chr4", "Chr5"]
        "Chr1-3,Chr8"     -> ["Chr1", "Chr2", "Chr3", "Chr8"]

    Args:
        spec: Chromosome specification string.
        available_chroms: List of available chromosome names from the reference.

    Returns:
        List of chromosome name strings.

    Raises:
        ValueError: If the spec is invalid or references unknown chromosomes.
    """
    spec = spec.strip()

    if spec.lower() == "all":
        return list(available_chroms)

    available_set = set(available_chroms)
    result: list[str] = []

    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue

        if "-" in part:
            # Parse range: "Chr2-5" or "Chr2-Chr5"
            tokens = part.split("-")
            if len(tokens) != 2:
                raise ValueError(
                    f"Invalid chromosome range '{part}'. "
                    f"Expected format like 'Chr2-5' or 'Chr2-Chr5'."
                )
            left = tokens[0].strip()
            right = tokens[1].strip()

            # Extract the non-numeric prefix from the left side
            # e.g., "Chr2" -> prefix="Chr", num="2"
            prefix_match = re.match(r"^([^\d]+)", left)
            prefix = prefix_match.group(1) if prefix_match else ""

            # If right side doesn't have the prefix, add it
            if not right.startswith(prefix):
                right = prefix + right

            # Extract numeric parts
            left_num_match = re.search(r"(\d+)$", left)
            right_num_match = re.search(r"(\d+)$", right)

            if not left_num_match or not right_num_match:
                raise ValueError(
                    f"Cannot parse numeric range from '{part}'. "
                    f"Expected format like 'Chr2-5' or 'Chr2-Chr5'."
                )

            left_num = int(left_num_match.group(1))
            right_num = int(right_num_match.group(1))

            if left_num > right_num:
                raise ValueError(
                    f"Invalid chromosome range '{part}': "
                    f"start ({left_num}) > end ({right_num})."
                )

            for i in range(left_num, right_num + 1):
                chrom = f"{prefix}{i}"
                result.append(chrom)
        else:
            # Single chromosome
            result.append(part)

    # Validate all chromosomes exist
    invalid = [c for c in result if c not in available_set]
    if invalid:
        raise ValueError(
            f"Chromosomes not found in reference genome: {invalid}. "
            f"Available chromosomes: {available_chroms}"
        )

    # Remove duplicates while preserving order
    seen: set[str] = set()
    unique_result: list[str] = []
    for chrom in result:
        if chrom not in seen:
            seen.add(chrom)
            unique_result.append(chrom)

    if not unique_result:
        raise ValueError(
            "No chromosomes specified. Use 'all' or specify chromosome names."
        )

    return unique_result


def extract_chromosomes(
    fasta_path: Path,
    chroms: list[str],
    output_path: Path,
) -> Path:
    """Extract specific chromosomes from a FASTA file into a new file.

    Args:
        fasta_path: Path to the input FASTA file.
        chroms: List of chromosome names to extract.
        output_path: Path to write the extracted FASTA.

    Returns:
        Path to the output FASTA file.

    Raises:
        FileNotFoundError: If input file does not exist.
        ValueError: If requested chromosomes are not found in the FASTA.
    """
    validate_file_exists(fasta_path)

    chrom_set = set(chroms)
    found: set[str] = set()
    records: list[str] = []

    current_name: str | None = None
    current_lines: list[str] = []

    def _flush() -> None:
        nonlocal current_name, current_lines
        if current_name is not None and current_name in chrom_set:
            found.add(current_name)
            records.extend(current_lines)

    with open(fasta_path, "r") as fin:
        for line in fin:
            if line.startswith(">"):
                _flush()
                current_name = line[1:].strip().split()[0]
                current_lines = [line]
            else:
                if current_name is not None:
                    current_lines.append(line)

        _flush()

    missing = chrom_set - found
    if missing:
        raise ValueError(
            f"Chromosomes not found in FASTA {fasta_path}: {sorted(missing)}"
        )

    with open(output_path, "w") as fout:
        fout.writelines(records)

    return output_path


def merge_fasta(
    base_fasta: Path,
    modified_fasta: Path,
    output_path: Path,
    modified_chroms: set[str],
) -> Path:
    """Merge a modified FASTA with the original, replacing specified chromosomes.

    Takes chromosomes from `modified_fasta` if they are in `modified_chroms`,
    otherwise takes chromosomes from `base_fasta`.

    Args:
        base_fasta: Path to the original (unmodified) FASTA.
        modified_fasta: Path to the FASTA with modified chromosomes.
        output_path: Path to write the merged FASTA.
        modified_chroms: Set of chromosome names to take from modified_fasta.

    Returns:
        Path to the merged output FASTA.
    """
    validate_file_exists(base_fasta)
    validate_file_exists(modified_fasta)

    # Read modified chromosomes
    modified_records: dict[str, list[str]] = {}
    current_name: str | None = None
    current_lines: list[str] = []

    with open(modified_fasta, "r") as fin:
        for line in fin:
            if line.startswith(">"):
                if current_name is not None:
                    modified_records[current_name] = current_lines
                current_name = line[1:].strip().split()[0]
                current_lines = [line]
            else:
                if current_name is not None:
                    current_lines.append(line)

        if current_name is not None:
            modified_records[current_name] = current_lines

    # Write merged FASTA
    with open(base_fasta, "r") as fin, open(output_path, "w") as fout:
        current_name = None
        current_lines = []

        def _flush_to(fout) -> None:
            if current_name is not None:
                if current_name in modified_chroms and current_name in modified_records:
                    fout.writelines(modified_records[current_name])
                else:
                    fout.writelines(current_lines)

        for line in fin:
            if line.startswith(">"):
                _flush_to(fout)
                current_name = line[1:].strip().split()[0]
                current_lines = [line]
            else:
                if current_name is not None:
                    current_lines.append(line)

        _flush_to(fout)

    return output_path


def write_fasta_record(
    file_handle,
    name: str,
    sequence: str,
    line_width: int = 60,
) -> None:
    """Write a single FASTA record to a file handle.

    Args:
        file_handle: Open writable file handle.
        name: Sequence name (goes after ">").
        sequence: Nucleotide sequence string.
        line_width: Number of characters per sequence line.
    """
    file_handle.write(f">{name}\n")
    for i in range(0, len(sequence), line_width):
        file_handle.write(sequence[i : i + line_width] + "\n")


def eprint(*args, **kwargs) -> None:
    """Print to stderr."""
    print(*args, file=sys.stderr, **kwargs)
