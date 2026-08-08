"""Build mixed germline and somatic TE benchmark read panels."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import logging
import os
import random
import shutil
import tempfile
from bisect import bisect_left, bisect_right
from collections import defaultdict
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path

import pysam
import tomllib
from Bio import SeqIO
from Bio.Seq import Seq

from simulate_data.longread import (
    GZIP_LEVEL,
    build_pbsim_command,
    contig_name_for,
    iter_maf_alignments,
    pool_fastq_gz,
    run_ccs,
)
from simulate_data.longread import (
    PLATFORMS as LONG_READ_PLATFORMS,
)
from simulate_data.modules.reads_illumina import _build_art_command
from simulate_data.modules.te_insertion import (
    _build_terandom_command,
    _prepare_known_del_file,
    _split_te_types,
)
from simulate_data.utils import (
    check_tool_installed,
    parse_chromosome_spec,
    run_command,
    validate_fasta,
    validate_file_exists,
)

logger = logging.getLogger(__name__)

CLASSES = ("homozygous", "heterozygous", "somatic_10", "somatic_20", "somatic_40")
COMPONENTS = (
    ("baseline", 0.60, frozenset()),
    ("clone40", 0.20, frozenset({"somatic_40"})),
    ("clone20", 0.10, frozenset({"somatic_40", "somatic_20"})),
    ("clone10", 0.10, frozenset({"somatic_40", "somatic_20", "somatic_10"})),
)
SOMATIC_CCF = {"somatic_10": 0.10, "somatic_20": 0.20, "somatic_40": 0.40}
TRUTH_FIELDS = (
    "event_id",
    "chrom",
    "position",
    "te_id",
    "te_family",
    "biological_class",
    "cellular_fraction",
    "expected_vaf",
    "hap1",
    "hap2",
    "strand",
    "tsd",
    "tsd_length",
    "consensus_length",
    "inserted_length",
    "sequence_identity",
    "snp_count",
    "indel_count",
    "length_change",
)


@dataclass
class Event:
    event_id: str
    chrom: str
    position: int
    te_id: str
    te_family: str
    biological_class: str
    cellular_fraction: float
    expected_vaf: float
    hap1: int
    hap2: int
    strand: str
    tsd: str
    tsd_length: int
    consensus_length: int
    inserted_length: int
    sequence_identity: float
    snp_count: int
    indel_count: int
    length_change: int
    sequence: str


def register_parser(parser):
    parser.add_argument("--config", help="Optional TOML panel configuration")
    parser.add_argument("--ref", help="Reference genome FASTA")
    parser.add_argument("--te", help="TE consensus FASTA")
    parser.add_argument(
        "--known-del", help="TEvarSim-compatible RepeatMasker annotation"
    )
    parser.add_argument("--output", required=True, help="Panel output directory")
    parser.add_argument("--stage", choices=("all", "catalog", "reads"), default="all")
    parser.add_argument(
        "--chroms", default=None, help="Chromosomes used for insertion placement"
    )
    parser.add_argument("--events-per-class", type=int, default=None)
    parser.add_argument("--coverage", type=float, action="append", default=None)
    parser.add_argument("--replicates", type=int, default=None)
    parser.add_argument(
        "--replicate", type=int, default=None, help="One replicate for an array task"
    )
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--read-length", type=int, default=None)
    parser.add_argument("--fragment-size", type=int, default=None)
    parser.add_argument("--fragment-std", type=int, default=None)
    parser.add_argument("--sequencing-system", default=None)
    parser.add_argument("--te-type", action="append", default=None)
    parser.add_argument("--snp-rate", type=float, default=None)
    parser.add_argument("--indel-rate", type=float, default=None)
    parser.add_argument("--tsd-min", type=int, default=None)
    parser.add_argument("--tsd-max", type=int, default=None)
    parser.add_argument("--sense-strand-ratio", type=float, default=None)


def _load_settings(args) -> dict:
    settings = {
        "chroms": "all",
        "events_per_class": 100,
        "coverage": [5.0, 15.0, 30.0],
        "replicates": 3,
        "seed": 916,
        "read_length": 150,
        "fragment_size": 300,
        "fragment_std": 30,
        "sequencing_system": "HSXn",
        "snp_rate": 0.0,
        "indel_rate": 0.0,
        "tsd_min": 3,
        "tsd_max": 5,
        "sense_strand_ratio": 0.5,
        "te_type": [],
    }
    if args.config:
        validate_file_exists(Path(args.config))
        with open(args.config, "rb") as handle:
            config = tomllib.load(handle)
        settings.update(config.get("panel", {}))
        settings.update(config.get("illumina", {}))
        settings.update(config.get("te", {}))
        inputs = config.get("inputs", {})
    else:
        inputs = {}
    for key in settings:
        value = getattr(args, key, None)
        if value is not None:
            settings[key] = value
    for key in ("ref", "te", "known_del"):
        settings[key] = getattr(args, key, None) or inputs.get(key)
    settings["te_type"] = _split_te_types(settings.get("te_type"))
    _validate_settings(settings, args.stage, args.replicate)
    return settings


def _validate_settings(s: dict, stage: str, replicate: int | None) -> None:
    if stage in {"all", "catalog"}:
        for key in ("ref", "te", "known_del"):
            if not s.get(key):
                raise ValueError(
                    f"--{key.replace('_', '-')} is required for the catalog stage"
                )
    if s["events_per_class"] <= 0 or s["replicates"] <= 0:
        raise ValueError("event and replicate counts must be positive")
    if not s["coverage"] or any(float(x) <= 0 for x in s["coverage"]):
        raise ValueError("coverage values must be positive")
    if not 0 <= s["sense_strand_ratio"] <= 1:
        raise ValueError("sense_strand_ratio must be between 0 and 1")
    for rate in ("snp_rate", "indel_rate"):
        if not 0 <= s[rate] <= 1:
            raise ValueError(f"{rate} must be between 0 and 1")
    if s["tsd_min"] <= 0 or s["tsd_max"] < s["tsd_min"]:
        raise ValueError("TSD bounds must be positive and ordered")
    if replicate is not None and not 1 <= replicate <= s["replicates"]:
        raise ValueError("--replicate is outside the configured replicate range")


def _read_fasta(path: Path) -> dict[str, str]:
    return {record.id: str(record.seq).upper() for record in SeqIO.parse(path, "fasta")}


def _te_family(te_id: str) -> str:
    return te_id.split("#", 1)[0].split("_", 1)[0]


def _sequence_stats(inserted: str, consensus: str) -> tuple[int, int, float]:
    length = min(len(inserted), len(consensus))
    snps = sum(a != b for a, b in zip(inserted[:length], consensus[:length]))
    indels = abs(len(inserted) - len(consensus))
    identity = 100.0 * max(length - snps, 0) / max(len(inserted), len(consensus), 1)
    return snps, indels, identity


def _allocate_counts(total: int, names: list[str]) -> dict[str, int]:
    base, remainder = divmod(total, len(names))
    return {name: base + (index < remainder) for index, name in enumerate(names)}


def _generate_catalog(settings: dict, output: Path) -> None:
    if (output / ".catalog.complete").exists():
        logger.info("Catalog already complete: %s", output)
        return
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(
            f"Refusing non-empty incomplete output directory: {output}"
        )
    output.mkdir(parents=True, exist_ok=True)
    ref_path, te_path, known_path = map(
        Path, (settings["ref"], settings["te"], settings["known_del"])
    )
    validate_fasta(ref_path)
    validate_fasta(te_path)
    validate_file_exists(known_path)
    check_tool_installed("tevarsim")
    reference = _read_fasta(ref_path)
    consensus = _read_fasta(te_path)
    chroms = parse_chromosome_spec(settings["chroms"], list(reference))
    total = settings["events_per_class"] * len(CLASSES)
    counts = _allocate_counts(total, chroms)
    raw_events: list[tuple[str, int, str, str]] = []
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        known = _prepare_known_del_file(known_path, tmpdir)
        for index, chrom in enumerate(chroms):
            prefix = str(tmpdir / f"events_{chrom}")
            cmd = _build_terandom_command(
                te_fasta=te_path,
                known_del=known,
                chrom=chrom,
                num_te=counts[chrom],
                outprefix=prefix,
                seed=settings["seed"] + index,
                ins_ratio=1.0,
                te_types=settings["te_type"],
                snp_rate=settings["snp_rate"],
                indel_rate=settings["indel_rate"],
                truncated_ratio=0.0,
                polya_ratio=0.0,
            )
            run_command(cmd)
            pool = _read_fasta(Path(prefix + ".fa"))
            with open(prefix + ".bed") as handle:
                for line in handle:
                    if not line.strip() or line.startswith("#"):
                        continue
                    fields = line.rstrip().split("\t")
                    pos, end, te_id = int(fields[1]), int(fields[2]), fields[3]
                    if end != pos + 1 or te_id not in pool:
                        continue
                    raw_events.append((chrom, pos, te_id, pool[te_id]))
    if len(raw_events) != total:
        raise RuntimeError(
            f"TEvarSim produced {len(raw_events)} insertions; expected exactly {total}"
        )
    rng = random.Random(settings["seed"])
    rng.shuffle(raw_events)
    class_labels = [
        label for label in CLASSES for _ in range(settings["events_per_class"])
    ]
    events: list[Event] = []
    for number, ((chrom, pos, te_id, sequence), label) in enumerate(
        zip(raw_events, class_labels), 1
    ):
        tsd_len = rng.randint(settings["tsd_min"], settings["tsd_max"])
        if pos < tsd_len or pos >= len(reference[chrom]):
            raise ValueError(
                f"Insertion position cannot support its TSD: {chrom}:{pos}"
            )
        strand = "+" if rng.random() < settings["sense_strand_ratio"] else "-"
        inserted = (
            sequence if strand == "+" else str(Seq(sequence).reverse_complement())
        )
        family = _te_family(te_id)
        source = next(
            (seq for name, seq in consensus.items() if _te_family(name) == family),
            sequence,
        )
        snps, indels, identity = _sequence_stats(sequence, source)
        if label == "homozygous":
            ccf, vaf, h1, h2 = 1.0, 1.0, 1, 1
        elif label == "heterozygous":
            ccf, vaf, h1, h2 = 1.0, 0.5, 1, 0
        else:
            ccf, vaf, h1, h2 = SOMATIC_CCF[label], SOMATIC_CCF[label] / 2, 1, 0
        events.append(
            Event(
                f"TE{number:06d}",
                chrom,
                pos,
                te_id,
                family,
                "somatic_insertion" if label.startswith("somatic") else label,
                ccf,
                vaf,
                h1,
                h2,
                strand,
                reference[chrom][pos - tsd_len : pos],
                tsd_len,
                len(source),
                len(inserted),
                identity,
                snps,
                indels,
                len(inserted) - len(source),
                inserted,
            )
        )
    _write_truth(events, output)
    _write_component_genomes(events, reference, output)
    shutil.copyfile(ref_path, output / "genomes" / "reference.fa")
    _write_metadata(settings, output, "catalog")
    (output / ".catalog.complete").touch()


def _event_group(event: Event) -> str:
    if event.biological_class != "somatic_insertion":
        return event.biological_class
    return f"somatic_{round(event.cellular_fraction * 100)}"


def _present(event: Event, component_classes: frozenset[str], hap: int) -> bool:
    group = _event_group(event)
    if group == "homozygous":
        return True
    if group == "heterozygous":
        return hap == 1
    return hap == 1 and group in component_classes


def _write_component_genomes(
    events: list[Event], reference: dict[str, str], output: Path
) -> None:
    genome_dir = output / "genomes"
    genome_dir.mkdir()
    coord_rows = []
    for component, weight, classes in COMPONENTS:
        for hap in (1, 2):
            path = genome_dir / f"{component}.hap{hap}.fa"
            with open(path, "w") as handle:
                for chrom, ref_seq in reference.items():
                    chrom_events = sorted(
                        (e for e in events if e.chrom == chrom),
                        key=lambda e: e.position,
                    )
                    chunks: list[str] = []
                    cursor = 0
                    offset = 0
                    for event in chrom_events:
                        coordinate = event.position + offset
                        is_present = _present(event, classes, hap)
                        coord_rows.append(
                            (
                                component,
                                weight,
                                hap,
                                event.event_id,
                                chrom,
                                coordinate,
                                int(is_present),
                            )
                        )
                        if is_present:
                            addition = event.sequence + event.tsd
                            chunks.append(ref_seq[cursor : event.position])
                            chunks.append(addition)
                            cursor = event.position
                            offset += len(addition)
                    chunks.append(ref_seq[cursor:])
                    sequence = "".join(chunks)
                    handle.write(f">{chrom}\n")
                    for start in range(0, len(sequence), 80):
                        handle.write(sequence[start : start + 80] + "\n")
    with open(output / "component_coordinates.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            (
                "component",
                "weight",
                "haplotype",
                "event_id",
                "chrom",
                "coordinate",
                "present",
            )
        )
        writer.writerows(coord_rows)


def _write_truth(events: list[Event], output: Path) -> None:
    with open(output / "truth_events.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TRUTH_FIELDS, delimiter="\t")
        writer.writeheader()
        for event in sorted(events, key=lambda e: (e.chrom, e.position)):
            row = asdict(event)
            row.pop("sequence")
            writer.writerow(row)
    with open(output / "truth_events.vcf", "w") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write('##INFO=<ID=TE,Number=1,Type=String,Description="TE family">\n')
        handle.write(
            "##INFO=<ID=CLASS,Number=1,Type=String,"
            'Description="Biological truth class">\n'
        )
        handle.write(
            '##INFO=<ID=CCF,Number=1,Type=Float,Description="Cellular fraction">\n'
        )
        handle.write(
            "##INFO=<ID=VAF,Number=1,Type=Float,"
            'Description="Expected insertion allele fraction">\n'
        )
        handle.write(
            "##INFO=<ID=TSD,Number=1,Type=String,"
            'Description="Target-site duplication">\n'
        )
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for e in sorted(events, key=lambda item: (item.chrom, item.position)):
            info = (
                f"TE={e.te_family};CLASS={e.biological_class};"
                f"CCF={e.cellular_fraction};VAF={e.expected_vaf};TSD={e.tsd}"
            )
            handle.write(
                f"{e.chrom}\t{e.position}\t{e.event_id}\tN\t<INS:{e.te_family}>\t.\tPASS\t{info}\n"
            )
    vcf = output / "truth_events.vcf"
    compressed = output / "truth_events.vcf.gz"
    pysam.tabix_compress(str(vcf), str(compressed), force=True)
    pysam.tabix_index(str(compressed), preset="vcf", force=True)


def _seed(
    base: int, coverage: float, replicate: int, component_index: int, hap: int
) -> int:
    value = f"{base}|{coverage:g}|{replicate}|{component_index}|{hap}".encode()
    return int(hashlib.sha256(value).hexdigest()[:8], 16) % 2_147_483_647 or 1


def _pool_fastq(source: Path, destination, prefix: str) -> int:
    count = 0
    with open(source) as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq, plus, qual = handle.readline(), handle.readline(), handle.readline()
            if not qual:
                raise ValueError(f"Truncated FASTQ: {source}")
            token, *rest = header.rstrip().split(maxsplit=1)
            destination.write(
                f"@{prefix}:{token[1:]}" + (f" {rest[0]}" if rest else "") + "\n"
            )
            destination.write(seq)
            destination.write(plus)
            destination.write(qual)
            count += 1
    return count


def _load_truth_events(path: Path) -> dict[str, dict]:
    with open(path) as handle:
        return {row["event_id"]: row for row in csv.DictReader(handle, delimiter="\t")}


def _load_component_coordinates(path: Path) -> dict[tuple[str, int], list[dict]]:
    result: dict[tuple[str, int], list[dict]] = defaultdict(list)
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            row["coordinate"] = int(row["coordinate"])
            row["present"] = int(row["present"])
            row["haplotype"] = int(row["haplotype"])
            result[(row["component"], row["haplotype"])].append(row)
    return result


def _build_coordinate_index(
    coordinate_rows: list[dict], truth: dict[str, dict]
) -> dict[str, tuple[list[int], list[tuple[int, str, bool]]]]:
    """Index insertion junctions per contig for span queries."""
    entries: dict[str, list[tuple[int, str, bool]]] = defaultdict(list)
    for row in coordinate_rows:
        coordinate = row["coordinate"]
        event_id = row["event_id"]
        is_present = bool(row["present"])
        entries[row["chrom"]].append((coordinate, event_id, is_present))
        if is_present:
            inserted = int(truth[event_id]["inserted_length"])
            tsd = int(truth[event_id]["tsd_length"])
            entries[row["chrom"]].append((coordinate + inserted + tsd, event_id, True))
    by_chrom = {}
    for chrom, rows in entries.items():
        rows.sort(key=lambda item: item[0])
        by_chrom[chrom] = ([item[0] for item in rows], rows)
    return by_chrom


def _tally_spans(
    spans,
    by_chrom: dict[str, tuple[list[int], list[tuple[int, str, bool]]]],
    counts: dict[str, list[int]],
    margin: int = 5,
) -> None:
    """Tally insertion-supporting and reference-spanning reads over spans.

    Each event is counted at most once per read, so a long read crossing
    both junctions of the same insertion still contributes a single count.
    """
    for chrom, start, end in spans:
        chrom_index = by_chrom.get(chrom)
        if chrom_index is None:
            continue
        positions, indexed = chrom_index
        left = bisect_left(positions, start + margin)
        right = bisect_right(positions, end - margin)
        seen: set[str] = set()
        for _coordinate, event_id, is_present in indexed[left:right]:
            if event_id in seen:
                continue
            counts[event_id][0 if is_present else 1] += 1
            seen.add(event_id)


def _count_origin_support(
    sam_path: Path,
    coordinate_rows: list[dict],
    truth: dict[str, dict],
    counts: dict[str, list[int]],
    margin: int = 5,
) -> None:
    """Count exact ART-origin insertion junctions and reference spanners."""
    by_chrom = _build_coordinate_index(coordinate_rows, truth)

    def _spans():
        with pysam.AlignmentFile(sam_path, "r") as sam:
            for record in sam.fetch(until_eof=True):
                if record.is_unmapped or record.reference_end is None:
                    continue
                yield (
                    sam.get_reference_name(record.reference_id),
                    record.reference_start,
                    record.reference_end,
                )

    _tally_spans(_spans(), by_chrom, counts, margin)


def _count_origin_support_maf(
    spans,
    coordinate_rows: list[dict],
    truth: dict[str, dict],
    counts: dict[str, list[int]],
    margin: int = 5,
) -> None:
    """Count PBSIM3-origin junctions and reference spanners from MAF spans.

    PBSIM3 emits no SAM, so long-read origin support comes from MAF
    ``(contig, start, end)`` tuples rather than alignment records. The
    tallying logic is otherwise identical to the short-read path.
    """
    by_chrom = _build_coordinate_index(coordinate_rows, truth)
    _tally_spans(spans, by_chrom, counts, margin)


def _write_observed_support(
    output_path: Path,
    truth: dict[str, dict],
    counts: dict[str, list[int]],
) -> None:
    fields = (
        "event_id",
        "biological_class",
        "cellular_fraction",
        "expected_vaf",
        "insertion_supporting_reads",
        "reference_spanning_reads",
        "observed_vaf",
    )
    summary: dict[tuple[str, str], list[tuple[float, float | None]]] = defaultdict(list)
    with open(output_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for event_id in sorted(truth):
            insertion, reference = counts[event_id]
            total = insertion + reference
            observed = insertion / total if total else None
            key = (
                truth[event_id]["biological_class"],
                truth[event_id]["cellular_fraction"],
            )
            summary[key].append((float(truth[event_id]["expected_vaf"]), observed))
            writer.writerow(
                {
                    "event_id": event_id,
                    "biological_class": truth[event_id]["biological_class"],
                    "cellular_fraction": truth[event_id]["cellular_fraction"],
                    "expected_vaf": truth[event_id]["expected_vaf"],
                    "insertion_supporting_reads": insertion,
                    "reference_spanning_reads": reference,
                    "observed_vaf": observed if observed is not None else "NA",
                }
            )
    with open(output_path.with_name("panel_summary.tsv"), "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            (
                "biological_class",
                "cellular_fraction",
                "event_count",
                "expected_vaf",
                "events_with_support",
                "mean_observed_vaf",
            )
        )
        for (label, ccf), values in sorted(summary.items()):
            observed = [value for _expected, value in values if value is not None]
            writer.writerow(
                (
                    label,
                    ccf,
                    len(values),
                    values[0][0],
                    len(observed),
                    sum(observed) / len(observed) if observed else "NA",
                )
            )


def _simulate_reads(
    settings: dict, output: Path, coverage: float, replicate: int
) -> None:
    if not (output / ".catalog.complete").exists():
        raise FileNotFoundError("Catalog is incomplete; run --stage catalog first")
    check_tool_installed("art_illumina")
    sample = f"cov{coverage:g}x_rep{replicate}"
    sample_dir = output / "reads" / sample
    if (sample_dir / ".complete").exists():
        logger.info("Read set already complete: %s", sample_dir)
        return
    if sample_dir.exists() and any(sample_dir.iterdir()):
        raise FileExistsError(f"Refusing incomplete read-set directory: {sample_dir}")
    sample_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    truth = _load_truth_events(output / "truth_events.tsv")
    coordinates = _load_component_coordinates(output / "component_coordinates.tsv")
    support_counts = defaultdict(lambda: [0, 0])
    with (
        tempfile.TemporaryDirectory() as tmp,
        open(sample_dir / "R1.fastq", "w") as r1_out,
        open(sample_dir / "R2.fastq", "w") as r2_out,
    ):
        tmpdir = Path(tmp)
        for component_index, (component, weight, _classes) in enumerate(COMPONENTS):
            for hap in (1, 2):
                component_cov = coverage * weight / 2
                seed = _seed(
                    settings["seed"], coverage, replicate, component_index, hap
                )
                prefix = tmpdir / f"{component}.hap{hap}."
                cmd = _build_art_command(
                    output / "genomes" / f"{component}.hap{hap}.fa",
                    settings["read_length"],
                    component_cov,
                    str(prefix),
                    True,
                    settings["fragment_size"],
                    settings["fragment_std"],
                    seed,
                    sequencing_system=settings["sequencing_system"],
                    sam=True,
                )
                run_command(cmd)
                sam_path = Path(str(prefix) + ".sam")
                _count_origin_support(
                    sam_path,
                    coordinates[(component, hap)],
                    truth,
                    support_counts,
                )
                name_prefix = f"{sample}:{component}:h{hap}:s{seed}"
                n1 = _pool_fastq(Path(str(prefix) + "1.fq"), r1_out, name_prefix)
                n2 = _pool_fastq(Path(str(prefix) + "2.fq"), r2_out, name_prefix)
                if n1 != n2:
                    raise RuntimeError(
                        f"ART produced unequal mate counts for {component} hap{hap}"
                    )
                rows.append((sample, component, weight, hap, component_cov, seed, n1))
                sam_path.unlink()
    control_dir = output / "reads" / f"{sample}_reference_control"
    control_dir.mkdir(exist_ok=False)
    control_seed = _seed(settings["seed"], coverage, replicate, 99, 1)
    prefix = control_dir / "control."
    run_command(
        _build_art_command(
            output / "genomes" / "reference.fa",
            settings["read_length"],
            coverage,
            str(prefix),
            True,
            settings["fragment_size"],
            settings["fragment_std"],
            control_seed,
            sequencing_system=settings["sequencing_system"],
            sam=True,
        )
    )
    Path(str(prefix) + "1.fq").rename(control_dir / "R1.fastq")
    Path(str(prefix) + "2.fq").rename(control_dir / "R2.fastq")
    Path(str(prefix) + ".sam").unlink(missing_ok=True)
    Path(str(prefix) + "1.aln").unlink(missing_ok=True)
    Path(str(prefix) + "2.aln").unlink(missing_ok=True)
    with open(sample_dir / "component_manifest.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            (
                "sample",
                "component",
                "cell_fraction",
                "haplotype",
                "art_coverage",
                "seed",
                "read_pairs",
            )
        )
        writer.writerows(rows)
    _write_observed_support(sample_dir / "observed_support.tsv", truth, support_counts)
    (sample_dir / ".complete").touch()
    (control_dir / ".complete").touch()


def _pbsim_contig_outputs(prefix: Path) -> list[tuple[Path, Path, Path | None]]:
    """Collect PBSIM3's per-contig outputs in contig order.

    PBSIM3 writes one output set per contig as ``<prefix>_0001.*``. Returns
    ``(ref, maf, reads_or_bam)`` triples; the third element is the .fq.gz for
    single-pass runs and the subread .bam for multi-pass runs.
    """
    outputs = []
    for ref in sorted(prefix.parent.glob(f"{prefix.name}_*.ref")):
        # Build sibling names by string, not Path.with_suffix: component
        # prefixes contain dots (e.g. "baseline.hap1"), so with_suffix would
        # strip part of the prefix rather than the ".ref" extension.
        stem = str(ref)[: -len(".ref")]
        maf = Path(stem + ".maf.gz")
        if not maf.is_file():
            raise FileNotFoundError(f"PBSIM3 MAF missing for {ref}")
        fastq = Path(stem + ".fq.gz")
        bam = Path(stem + ".bam")
        payload = fastq if fastq.is_file() else (bam if bam.is_file() else None)
        if payload is None:
            raise FileNotFoundError(f"PBSIM3 produced no reads for {ref}")
        outputs.append((ref, maf, payload))
    if not outputs:
        raise FileNotFoundError(f"PBSIM3 produced no output for prefix {prefix}")
    return outputs


def _simulate_long_reads(
    settings: dict,
    output: Path,
    coverage: float,
    replicate: int,
    platform_key: str,
) -> None:
    """Simulate one long-read sample as a mixture of the component genomes.

    Mirrors :func:`_simulate_reads`: same components, same seed derivation,
    same manifest and QC schema. Differs in being single-end, deriving origin
    support from MAF rather than SAM, and writing gzip throughout.
    """
    if not (output / ".catalog.complete").exists():
        raise FileNotFoundError("Catalog is incomplete; run --stage catalog first")
    if platform_key not in LONG_READ_PLATFORMS:
        raise ValueError(
            f"Unknown long-read platform: {platform_key}; "
            f"choose from {', '.join(sorted(LONG_READ_PLATFORMS))}"
        )
    platform = LONG_READ_PLATFORMS[platform_key]
    check_tool_installed("pbsim")
    if platform.needs_ccs:
        check_tool_installed("ccs")

    threads = int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))
    min_yield = float(settings.get("min_ccs_yield", 0.5))

    sample = f"cov{coverage:g}x_rep{replicate}"
    sample_dir = output / "reads" / platform_key / sample
    if (sample_dir / ".complete").exists():
        logger.info("Read set already complete: %s", sample_dir)
        return
    if sample_dir.exists() and any(sample_dir.iterdir()):
        raise FileExistsError(f"Refusing incomplete read-set directory: {sample_dir}")
    sample_dir.mkdir(parents=True, exist_ok=True)

    truth = _load_truth_events(output / "truth_events.tsv")
    coordinates = _load_component_coordinates(output / "component_coordinates.tsv")
    support_counts = defaultdict(lambda: [0, 0])
    rows = []

    with (
        tempfile.TemporaryDirectory() as tmp,
        gzip.open(
            sample_dir / "reads.fastq.gz", "wt", compresslevel=GZIP_LEVEL
        ) as pooled,
    ):
        tmpdir = Path(tmp)
        for component_index, (component, weight, _classes) in enumerate(COMPONENTS):
            for hap in (1, 2):
                component_cov = coverage * weight / 2
                seed = _seed(
                    settings["seed"], coverage, replicate, component_index, hap
                )
                prefix = tmpdir / f"{component}.hap{hap}"
                run_command(
                    build_pbsim_command(
                        platform,
                        output / "genomes" / f"{component}.hap{hap}.fa",
                        component_cov,
                        prefix,
                        seed,
                    )
                )

                name_prefix = f"{sample}:{component}:h{hap}:s{seed}"
                reads_written = 0
                yields = []
                for ref, maf, payload in _pbsim_contig_outputs(prefix):
                    contig = contig_name_for(ref)
                    _count_origin_support_maf(
                        (
                            (contig, start, end)
                            for start, end in iter_maf_alignments(maf)
                        ),
                        coordinates[(component, hap)],
                        truth,
                        support_counts,
                    )
                    if payload.suffix == ".bam":
                        fastq = Path(str(payload)[: -len(".bam")] + ".fastq.gz")
                        stats = run_ccs(payload, fastq, threads, min_yield)
                        yields.append(stats["yield_fraction"])
                        payload.unlink()
                    else:
                        fastq = payload
                    reads_written += pool_fastq_gz(fastq, pooled, name_prefix)

                rows.append(
                    (
                        sample,
                        platform_key,
                        component,
                        weight,
                        hap,
                        f"{component_cov:g}",
                        seed,
                        reads_written,
                        f"{sum(yields) / len(yields):.4f}" if yields else "NA",
                    )
                )

    _simulate_long_read_control(
        settings, output, coverage, replicate, platform, threads, min_yield
    )

    with open(sample_dir / "component_manifest.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            (
                "sample",
                "platform",
                "component",
                "cell_fraction",
                "haplotype",
                "pbsim_coverage",
                "seed",
                "reads",
                "ccs_yield",
            )
        )
        writer.writerows(rows)
    _write_observed_support(sample_dir / "observed_support.tsv", truth, support_counts)
    (sample_dir / ".complete").touch()


def _simulate_long_read_control(
    settings: dict,
    output: Path,
    coverage: float,
    replicate: int,
    platform,
    threads: int,
    min_yield: float,
) -> None:
    """Simulate the matched reference-only control for a long-read sample."""
    sample = f"cov{coverage:g}x_rep{replicate}"
    control_dir = output / "reads" / platform.key / f"{sample}_reference_control"
    if (control_dir / ".complete").exists():
        return
    if control_dir.exists() and any(control_dir.iterdir()):
        raise FileExistsError(f"Refusing incomplete control directory: {control_dir}")
    control_dir.mkdir(parents=True, exist_ok=True)

    control_seed = _seed(settings["seed"], coverage, replicate, 99, 1)
    with (
        tempfile.TemporaryDirectory() as tmp,
        gzip.open(
            control_dir / "reads.fastq.gz", "wt", compresslevel=GZIP_LEVEL
        ) as pooled,
    ):
        prefix = Path(tmp) / "control"
        run_command(
            build_pbsim_command(
                platform,
                output / "genomes" / "reference.fa",
                coverage,
                prefix,
                control_seed,
            )
        )
        name_prefix = f"{sample}:control:s{control_seed}"
        for _ref, _maf, payload in _pbsim_contig_outputs(prefix):
            if payload.suffix == ".bam":
                fastq = Path(str(payload)[: -len(".bam")] + ".fastq.gz")
                run_ccs(payload, fastq, threads, min_yield)
                payload.unlink()
            else:
                fastq = payload
            pool_fastq_gz(fastq, pooled, name_prefix)
    (control_dir / ".complete").touch()


def _write_panel_manifest(settings: dict, output: Path) -> None:
    path = output / "panel_manifest.tsv"
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            ("sample", "coverage", "replicate", "r1", "r2", "control_r1", "control_r2")
        )
        for coverage in settings["coverage"]:
            for rep in range(1, settings["replicates"] + 1):
                sample = f"cov{float(coverage):g}x_rep{rep}"
                writer.writerow(
                    (
                        sample,
                        coverage,
                        rep,
                        f"reads/{sample}/R1.fastq",
                        f"reads/{sample}/R2.fastq",
                        f"reads/{sample}_reference_control/R1.fastq",
                        f"reads/{sample}_reference_control/R2.fastq",
                    )
                )


def _write_metadata(settings: dict, output: Path, stage: str) -> None:
    checksums = {}
    for key in ("ref", "te", "known_del"):
        path = Path(settings[key]) if settings.get(key) else None
        if path is not None and path.is_file():
            digest = hashlib.sha256()
            with open(path, "rb") as handle:
                for block in iter(lambda: handle.read(1024 * 1024), b""):
                    digest.update(block)
            checksums[key] = {"path": str(path), "sha256": digest.hexdigest()}
    metadata = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "stage": stage,
        "settings": settings,
        "input_checksums": checksums,
        "tools": {
            name: shutil.which(name)
            for name in ("tevarsim", "art_illumina")
            if shutil.which(name)
        },
    }
    with open(output / "run_metadata.json", "w") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
        handle.write("\n")


def main(args):
    settings = _load_settings(args)
    output = Path(args.output)
    if args.stage in {"all", "catalog"}:
        _generate_catalog(settings, output)
        _write_panel_manifest(settings, output)
    if args.stage in {"all", "reads"}:
        coverages = settings["coverage"]
        replicates = (
            [args.replicate] if args.replicate else range(1, settings["replicates"] + 1)
        )
        for coverage in coverages:
            for replicate in replicates:
                _simulate_reads(settings, output, float(coverage), replicate)
