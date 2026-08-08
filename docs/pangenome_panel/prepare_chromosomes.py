#!/usr/bin/env python3
"""Extract one chromosome per assembly and rename it to PanSN form.

MAGIC16 assemblies name sequences by GenBank accession (``CM020633.1``), so
chromosome identity is resolved from each assembly's ``*_assembly_report.txt``
rather than by parsing FASTA description text, which is not a stable contract.

Output records are named ``<sample>#<haplotype>#<contig>`` (PanSN), the naming
pangenome graph tools expect. Applying it here rather than before graph
construction avoids invalidating every coordinate later.
"""

from __future__ import annotations

import argparse
import gzip
import tomllib
from pathlib import Path

PANSN_SEPARATOR = "#"
LINE_WIDTH = 60


def pansn_name(sample: str, contig: str, haplotype: int = 1) -> str:
    """Build a PanSN sequence name: ``sample#haplotype#contig``."""
    for part, label in ((sample, "sample"), (contig, "contig")):
        if PANSN_SEPARATOR in part:
            raise ValueError(
                f"PanSN {label} name may not contain {PANSN_SEPARATOR!r}: {part!r}"
            )
    return f"{sample}{PANSN_SEPARATOR}{haplotype}{PANSN_SEPARATOR}{contig}"


def accession_for_chromosome(report_path: Path, chromosome: str) -> str:
    """Return the GenBank accession for an assembled chromosome.

    Reads the NCBI assembly report, matching on Assigned-Molecule and
    requiring the sequence role to be ``assembled-molecule`` so unplaced
    scaffolds can never satisfy the lookup.
    """
    with report_path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            sequence_role, assigned_molecule, genbank = (
                fields[1],
                fields[2],
                fields[4],
            )
            if sequence_role != "assembled-molecule":
                continue
            if assigned_molecule == chromosome:
                return genbank
    raise ValueError(
        f"No assembled molecule named {chromosome!r} in {report_path}"
    )


def extract_record(fasta_gz: Path, accession: str) -> str:
    """Return the sequence of one record from a gzipped FASTA."""
    chunks: list[str] = []
    capturing = False
    with gzip.open(fasta_gz, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if capturing:
                    break
                capturing = line[1:].split()[0] == accession
                continue
            if capturing:
                chunks.append(line.strip())
    if not chunks:
        raise ValueError(f"Record {accession} not found in {fasta_gz}")
    return "".join(chunks).upper()


def write_fasta(path: Path, name: str, sequence: str) -> None:
    with path.open("w") as handle:
        handle.write(f">{name}\n")
        for start in range(0, len(sequence), LINE_WIDTH):
            handle.write(sequence[start : start + LINE_WIDTH] + "\n")


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    with args.config.open("rb") as handle:
        config = tomllib.load(handle)

    chrom = config["panel"]["chrom"]
    raw_dir = args.output / "genomes" / "raw"
    chr_dir = args.output / "genomes" / chrom.lower()
    chr_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for genome in config["genomes"]:
        name = genome["name"]
        directory = genome["directory"]
        destination = chr_dir / f"{name}.{chrom.lower()}.fa"
        sentinel = chr_dir / f".{name}.complete"

        if sentinel.exists() and destination.is_file():
            sequence_length = sum(
                len(line.strip())
                for line in destination.open()
                if not line.startswith(">")
            )
            print(f"{name}\tcached\t{sequence_length}")
            rows.append((name, destination, sequence_length))
            continue

        report = raw_dir / name / f"{directory}_assembly_report.txt"
        fasta_gz = raw_dir / name / f"{directory}_genomic.fna.gz"
        for required in (report, fasta_gz):
            if not required.is_file():
                raise FileNotFoundError(f"Missing download for {name}: {required}")

        accession = accession_for_chromosome(report, genome["chromosome"])
        sequence = extract_record(fasta_gz, accession)
        record_name = pansn_name(name, chrom)
        write_fasta(destination, record_name, sequence)
        sentinel.touch()

        print(f"{name}\t{accession}\t{record_name}\t{len(sequence)}")
        rows.append((name, destination, len(sequence)))

    with (chr_dir / "chromosome_manifest.tsv").open("w") as handle:
        handle.write("genome\tpansn_name\tpath\tlength\n")
        for name, path, length in rows:
            handle.write(f"{name}\t{pansn_name(name, chrom)}\t{path}\t{length}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
