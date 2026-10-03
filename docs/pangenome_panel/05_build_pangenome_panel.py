#!/usr/bin/env python3
"""Build the multi-genome TE panel: anchors, sharing matrix, augmented genomes.

Ground truth for RelocaTE3 multi-FASTA support. Synthetic riceTElib insertions
are distributed across all 15 non-empty sharing patterns of four rice genomes,
so every "is this TE in the references, and which ones?" answer is known.

The central constraint is orthology. The genomes differ by megabases of
indels, so a shared insertion is only meaningful at orthologous coordinates.
Anchors are chosen in the reference genome and lifted into every other genome
with ``paftools.js liftover``; an anchor is rejected unless it lifts to
exactly one position in every genome. A weaker filter yields a silently wrong
truth set, which is far worse than a smaller panel.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
import re
import subprocess
import tomllib
from bisect import bisect_left
from dataclasses import dataclass, field
from datetime import datetime, timezone
from itertools import product
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

LINE_WIDTH = 60

IUPAC = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "S": frozenset("GC"),
    "W": frozenset("AT"),
    "K": frozenset("GT"),
    "M": frozenset("AC"),
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}

TRUTH_FIELDS = (
    "event_id",
    "te_id",
    "te_family",
    "te_group",
    "pattern",
    "n_genomes",
    "genomes",
    "anchor_id",
    "strand",
    "tsd_length",
    "inserted_length",
)


# --------------------------------------------------------------------------
# Sharing patterns
# --------------------------------------------------------------------------


def sharing_patterns(n_genomes: int) -> list[tuple[bool, ...]]:
    """All non-empty presence/absence patterns over ``n_genomes`` genomes.

    For four genomes this is 15 patterns: 4 private, 6 pairs, 4 triples, and
    1 core. Ordered by cardinality then by genome index for readability.
    """
    patterns = [
        pattern
        for pattern in product((False, True), repeat=n_genomes)
        if any(pattern)
    ]
    return sorted(patterns, key=lambda p: (sum(p), tuple(not x for x in p)))


def pattern_label(pattern: tuple[bool, ...]) -> str:
    """Render a pattern as a bitmask string, e.g. ``1010``."""
    return "".join("1" if present else "0" for present in pattern)


# --------------------------------------------------------------------------
# Liftover
# --------------------------------------------------------------------------


def parse_liftover_bed(path: Path) -> dict[str, list[int]]:
    """Parse a ``paftools.js liftover`` BED into {source_key: [start, ...]}.

    paftools.js **discards the input BED name** and rewrites column 4 as
    ``<query_chrom>_<start>_<end>``, so anchors are keyed by their original
    reference coordinate rather than by the name we supplied. See
    ``source_key`` for the matching key builder.

    More than one interval for a key means the site is ambiguous in that
    genome; the caller rejects those.
    """
    lifted: dict[str, list[int]] = {}
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            start, name = int(fields[1]), fields[3]
            lifted.setdefault(name, []).append(start)
    return lifted


def source_key(record: str, start: int) -> str:
    """Rebuild the key paftools.js writes into column 4 of its output."""
    return f"{record}_{start}_{start + 1}"


def anchor_is_usable(
    lifted: dict[str, list[int]], targets: list[str]
) -> bool:
    """True when an anchor lifts to exactly one position in every genome."""
    for target in targets:
        positions = lifted.get(target)
        if not positions or len(positions) != 1:
            return False
    return True


def run_liftover(
    paf: Path,
    bed: Path,
    output: Path,
    min_length: int,
    min_quality: int,
) -> Path:
    """Lift a BED of reference anchors into one target genome."""
    with output.open("w") as handle:
        subprocess.run(
            [
                "paftools.js",
                "liftover",
                "-l",
                str(min_length),
                "-q",
                str(min_quality),
                str(paf),
                str(bed),
            ],
            stdout=handle,
            check=True,
        )
    return output


# --------------------------------------------------------------------------
# Site selection
# --------------------------------------------------------------------------


def iupac_matches(sequence: str, pattern: str) -> bool:
    sequence, pattern = sequence.upper(), pattern.upper()
    return len(sequence) == len(pattern) and all(
        base in IUPAC.get(code, frozenset())
        for base, code in zip(sequence, pattern)
    )


def site_is_usable(
    reference: str,
    position: int,
    tsd_length: int,
    group: dict,
    flank_acgt: int,
) -> bool:
    """Reference-side site check: unambiguous flanks, TSD motif, context.

    Mirrors the equivalent check in build_multite_panel.py so the two panels
    place insertions under the same rules.
    """
    context = group.get("insertion_context")
    left_context, right_context = context.split("|", 1) if context else ("", "")
    if position < max(tsd_length, flank_acgt, len(left_context)):
        return False
    if position + max(flank_acgt, len(right_context)) >= len(reference):
        return False
    if flank_acgt and set(
        reference[position - flank_acgt : position + flank_acgt]
    ) - set("ACGT"):
        return False

    tsd = reference[position - tsd_length : position] if tsd_length else ""
    motif = group.get("tsd_motif")
    if motif and not iupac_matches(tsd, motif):
        return False

    if context:
        if not iupac_matches(
            reference[position - len(left_context) : position], left_context
        ):
            return False
        if not iupac_matches(
            reference[position : position + len(right_context)], right_context
        ):
            return False
    return True


def far_enough(positions: list[int], position: int, distance: int) -> bool:
    index = bisect_left(positions, position)
    if index and position - positions[index - 1] < distance:
        return False
    if index < len(positions) and positions[index] - position < distance:
        return False
    return True


# --------------------------------------------------------------------------
# Background TE interference
# --------------------------------------------------------------------------


def load_background_bed(path: Path) -> list[tuple[str, int, int]]:
    """Load a background TE BED as (family, start, end) tuples."""
    intervals = []
    if not path.is_file():
        return intervals
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            intervals.append((fields[3], int(fields[1]), int(fields[2])))
    return intervals


def site_clear_of_background(
    positions: dict[str, int],
    te_family: str,
    background: dict[str, list[tuple[str, int, int]]],
    buffer_bp: int,
) -> bool:
    """True when no genome already carries this TE family near the site.

    A pre-existing element of the same family close to a synthetic insertion
    would make a correct RelocaTE3 call look like a false positive, so a
    collision in ANY genome disqualifies the site everywhere.
    """
    for genome, position in positions.items():
        if position is None:
            continue
        for family, start, end in background.get(genome, ()):
            if family != te_family:
                continue
            if start - buffer_bp <= position <= end + buffer_bp:
                return False
    return True


# --------------------------------------------------------------------------
# Insertion
# --------------------------------------------------------------------------


def insert_all(sequence: str, events: list[tuple[int, str, str]]) -> str:
    """Insert (position, te_sequence, tsd) events into a sequence.

    Applied in descending coordinate order so an insertion never shifts the
    coordinates of one not yet applied. Each insertion contributes the TE
    followed by a duplicate of the target site.
    """
    result = sequence
    for position, te_sequence, tsd in sorted(events, reverse=True):
        result = result[:position] + te_sequence + tsd + result[position:]
    return result


# --------------------------------------------------------------------------
# Panel construction
# --------------------------------------------------------------------------


@dataclass
class Candidate:
    te_id: str
    te_family: str
    classification: str
    sequence: str


@dataclass
class Event:
    event_id: str
    te_id: str
    te_family: str
    te_group: str
    pattern: tuple[bool, ...]
    anchor_id: str
    strand: str
    tsd_length: int
    inserted_length: int
    sequence: str
    positions: dict[str, int | None] = field(default_factory=dict)
    tsds: dict[str, str] = field(default_factory=dict)


def read_fasta_single(path: Path) -> tuple[str, str]:
    record = next(SeqIO.parse(path, "fasta"))
    return record.id, str(record.seq).upper()


def read_te_library(path: Path) -> list[Candidate]:
    records = []
    for record in SeqIO.parse(path, "fasta"):
        if "#" not in record.id:
            continue
        family, classification = record.id.split("#", 1)
        sequence = str(record.seq).upper()
        records.append(Candidate(record.id, family, classification, sequence))
    return records


def candidates_for_group(
    records: list[Candidate], group: dict, selection: dict
) -> list[Candidate]:
    wanted = set(group["classifications"])
    usable = []
    for record in records:
        if record.classification not in wanted:
            continue
        if not (
            selection["min_te_length"]
            <= len(record.sequence)
            <= selection["max_te_length"]
        ):
            continue
        if set(record.sequence) - set("ACGT"):
            continue
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", record.te_family):
            continue
        usable.append(record)
    return usable


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_fasta(path: Path, name: str, sequence: str) -> None:
    with path.open("w") as handle:
        handle.write(f">{name}\n")
        for start in range(0, len(sequence), LINE_WIDTH):
            handle.write(sequence[start : start + LINE_WIDTH] + "\n")


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Report anchor yield without writing augmented genomes.",
    )
    return parser


def load_settings(config_path: Path) -> dict:
    with config_path.open("rb") as handle:
        config = tomllib.load(handle)

    group_config_path = Path(config["inputs"]["te_group_config"])
    with group_config_path.open("rb") as handle:
        group_config = tomllib.load(handle)

    genomes = config["genomes"]
    reference = next(g["name"] for g in genomes if g.get("reference"))
    return {
        "te_library": Path(config["inputs"]["te_library"]),
        "te_groups": group_config["te_groups"],
        "chrom": config["panel"]["chrom"],
        "events_per_pattern": config["panel"]["events_per_pattern"],
        "seed": config["panel"]["seed"],
        "genomes": [g["name"] for g in genomes],
        "reference": reference,
        "targets": [g["name"] for g in genomes if not g.get("reference")],
        "orthology": config["orthology"],
        "selection": config["selection"],
        "config_path": str(config_path),
    }


def build_anchors(settings: dict, output: Path) -> tuple[dict, dict]:
    """Sample reference anchors and lift them into every other genome.

    Returns (accepted anchors, acceptance statistics).
    """
    chrom = settings["chrom"].lower()
    chr_dir = output / "genomes" / chrom
    ortho_dir = output / "orthology"
    ortho_dir.mkdir(parents=True, exist_ok=True)

    reference = settings["reference"]
    _name, ref_sequence = read_fasta_single(chr_dir / f"{reference}.{chrom}.fa")

    selection = settings["selection"]
    rng = random.Random(settings["seed"])
    groups = settings["te_groups"]
    patterns = sharing_patterns(len(settings["genomes"]))
    needed = settings["events_per_pattern"] * len(patterns)

    # Oversample: the orthology filter rejects a large fraction, and the
    # rejection rate is not known until the liftover has run.
    oversample = int(needed * 20)

    candidates: list[tuple[str, int, int, dict]] = []
    taken: list[int] = []
    attempts = 0
    max_attempts = oversample * 50
    while len(candidates) < oversample and attempts < max_attempts:
        attempts += 1
        group = groups[len(candidates) % len(groups)]
        tsd_length = rng.randint(group["tsd_min"], group["tsd_max"])
        position = rng.randrange(
            max(tsd_length, selection["flank_acgt"], 1),
            len(ref_sequence) - selection["flank_acgt"] - 1,
        )
        if not far_enough(taken, position, selection["min_distance"]):
            continue
        if not site_is_usable(
            ref_sequence, position, tsd_length, group, selection["flank_acgt"]
        ):
            continue
        index = bisect_left(taken, position)
        taken.insert(index, position)
        anchor_id = f"anchor_{len(candidates):06d}"
        candidates.append((anchor_id, position, tsd_length, group))

    bed_path = ortho_dir / "candidate_anchors.bed"
    ref_record = f"{reference}#1#{settings['chrom']}"
    with bed_path.open("w") as handle:
        for anchor_id, position, _tsd, _group in candidates:
            handle.write(
                f"{ref_record}\t{position}\t{position + 1}\t{anchor_id}\n"
            )

    # paftools.js rewrites the BED name column, so anchors are recovered by
    # their original reference coordinate.
    key_to_anchor = {
        source_key(ref_record, position): anchor_id
        for anchor_id, position, _tsd, _group in candidates
    }

    lifted_by_target: dict[str, dict[str, list[int]]] = {}
    for target in settings["targets"]:
        paf = ortho_dir / f"{reference}_to_{target}.paf"
        if not paf.is_file():
            raise FileNotFoundError(f"Orthology PAF is missing: {paf}")
        out_bed = ortho_dir / f"lifted_{target}.bed"
        run_liftover(
            paf,
            bed_path,
            out_bed,
            settings["orthology"]["min_alignment_length"],
            settings["orthology"]["min_mapping_quality"],
        )
        by_source = parse_liftover_bed(out_bed)
        lifted_by_target[target] = {
            key_to_anchor[key]: starts
            for key, starts in by_source.items()
            if key in key_to_anchor
        }

    background = {}
    annotation_dir = output / "annotation"
    for genome in settings["genomes"]:
        background[genome] = load_background_bed(
            annotation_dir / f"{genome}.background_te.bed"
        )

    stats = {
        "candidates_sampled": len(candidates),
        "rejected_unlifted": 0,
        "rejected_multimapping": 0,
        "rejected_background_te": 0,
        "accepted": 0,
        "needed": needed,
    }

    accepted = {}
    for anchor_id, position, tsd_length, group in candidates:
        per_target = {
            target: lifted_by_target[target].get(anchor_id, [])
            for target in settings["targets"]
        }
        if any(len(v) == 0 for v in per_target.values()):
            stats["rejected_unlifted"] += 1
            continue
        if any(len(v) > 1 for v in per_target.values()):
            stats["rejected_multimapping"] += 1
            continue
        if not anchor_is_usable(per_target, settings["targets"]):
            stats["rejected_unlifted"] += 1
            continue

        positions = {reference: position}
        positions.update({t: per_target[t][0] for t in settings["targets"]})
        accepted[anchor_id] = {
            "positions": positions,
            "tsd_length": tsd_length,
            "group": group,
        }
        stats["accepted"] += 1

    stats["background"] = background
    return accepted, stats


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    settings = load_settings(args.config)
    output = args.output

    anchors, stats = build_anchors(settings, output)

    rate = stats["accepted"] / max(stats["candidates_sampled"], 1)
    print(f"candidates_sampled\t{stats['candidates_sampled']}")
    print(f"rejected_unlifted\t{stats['rejected_unlifted']}")
    print(f"rejected_multimapping\t{stats['rejected_multimapping']}")
    print(f"accepted\t{stats['accepted']}")
    print(f"acceptance_rate\t{rate:.4f}")
    print(f"needed\t{stats['needed']}")

    if args.validate_only:
        if stats["accepted"] < stats["needed"]:
            print(
                f"WARNING: only {stats['accepted']} anchors survived; "
                f"{stats['needed']} needed. Reduce events_per_pattern rather "
                "than relaxing the orthology filter."
            )
        return 0

    build_panel(settings, output, anchors, stats)

    # Written after the panel build so the background-TE rejections, which are
    # only counted during TE assignment, are included.
    with (output / "orthology" / "anchor_acceptance.tsv").open(
        "w", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("metric", "count"))
        for key in (
            "candidates_sampled",
            "rejected_unlifted",
            "rejected_multimapping",
            "rejected_background_te",
            "accepted",
            "needed",
        ):
            writer.writerow((key, stats[key]))
    print(f"rejected_background_te\t{stats['rejected_background_te']}")
    return 0


def build_panel(settings: dict, output: Path, anchors: dict, stats: dict) -> None:
    """Assign anchors to sharing patterns and write augmented genomes."""
    complete = output / ".panel.complete"
    if complete.exists():
        print(f"Panel already complete: {output}")
        return

    rng = random.Random(settings["seed"] + 1)
    chrom = settings["chrom"]
    chrom_key = chrom.lower()
    chr_dir = output / "genomes" / chrom_key
    augmented_dir = output / "genomes" / "augmented"
    augmented_dir.mkdir(parents=True, exist_ok=True)

    sequences = {
        genome: read_fasta_single(chr_dir / f"{genome}.{chrom_key}.fa")
        for genome in settings["genomes"]
    }

    library = read_te_library(settings["te_library"])
    by_group = {
        group["name"]: candidates_for_group(library, group, settings["selection"])
        for group in settings["te_groups"]
    }

    patterns = sharing_patterns(len(settings["genomes"]))
    per_pattern = settings["events_per_pattern"]

    anchor_ids = sorted(anchors)
    rng.shuffle(anchor_ids)
    available = min(len(anchor_ids), per_pattern * len(patterns))
    per_pattern = min(per_pattern, available // len(patterns))
    if per_pattern == 0:
        raise RuntimeError(
            f"Only {len(anchor_ids)} anchors survived the orthology filter; "
            f"not enough for {len(patterns)} sharing patterns"
        )

    background = stats["background"]
    buffer_bp = settings["selection"]["background_te_buffer"]

    used_exemplars: dict[str, set[str]] = {name: set() for name in by_group}
    events: list[Event] = []
    cursor = 0
    number = 0
    for pattern_index, pattern in enumerate(patterns):
        for slot in range(per_pattern):
            # Round-robin TE groups across patterns so TE family and sharing
            # pattern stay uncorrelated.
            group = settings["te_groups"][
                (pattern_index * per_pattern + slot) % len(settings["te_groups"])
            ]
            pool = [
                candidate
                for candidate in by_group[group["name"]]
                if candidate.te_id not in used_exemplars[group["name"]]
            ]
            if not pool:
                pool = by_group[group["name"]]
                used_exemplars[group["name"]].clear()
            candidate = rng.choice(pool)
            used_exemplars[group["name"]].add(candidate.te_id)

            # Advance past anchors where this TE family already sits nearby in
            # some genome: a real element of the same family next to a
            # synthetic one would make a correct call look like a false
            # positive. The family is only known now, which is why this filter
            # lives here rather than in anchor selection.
            anchor_id = None
            while cursor < len(anchor_ids):
                trial_id = anchor_ids[cursor]
                cursor += 1
                carrying = {
                    genome: anchors[trial_id]["positions"][genome]
                    for index, genome in enumerate(settings["genomes"])
                    if pattern[index]
                }
                if site_clear_of_background(
                    carrying, candidate.te_family, background, buffer_bp
                ):
                    anchor_id = trial_id
                    break
                stats["rejected_background_te"] += 1
            if anchor_id is None:
                raise RuntimeError(
                    "Ran out of anchors while avoiding background TEs; "
                    "reduce events_per_pattern or increase oversampling"
                )
            anchor = anchors[anchor_id]

            tsd_length = rng.randint(group["tsd_min"], group["tsd_max"])
            strand = (
                "+"
                if rng.random() < settings["selection"]["sense_strand_ratio"]
                else "-"
            )
            inserted = (
                candidate.sequence
                if strand == "+"
                else str(Seq(candidate.sequence).reverse_complement())
            )

            number += 1
            event = Event(
                event_id=f"PGTE{number:06d}",
                te_id=candidate.te_id,
                te_family=candidate.te_family,
                te_group=group["name"],
                pattern=pattern,
                anchor_id=anchor_id,
                strand=strand,
                tsd_length=tsd_length,
                inserted_length=len(inserted),
                sequence=inserted,
            )
            for index, genome in enumerate(settings["genomes"]):
                if not pattern[index]:
                    event.positions[genome] = None
                    continue
                position = anchor["positions"][genome]
                event.positions[genome] = position
                # The TSD is taken from each genome's own local sequence: the
                # flanking bases differ between genomes.
                _name, sequence = sequences[genome]
                event.tsds[genome] = (
                    sequence[position - tsd_length : position]
                    if tsd_length
                    else ""
                )
            events.append(event)

    for index, genome in enumerate(settings["genomes"]):
        name, sequence = sequences[genome]
        insertions = [
            (
                event.positions[genome],
                event.sequence,
                event.tsds.get(genome, ""),
            )
            for event in events
            if event.pattern[index]
        ]
        augmented = insert_all(sequence, insertions)
        write_fasta(
            augmented_dir / f"{genome}.{chrom_key}.te.fa", name, augmented
        )
        print(
            f"{genome}\tinsertions={len(insertions)}\t"
            f"length {len(sequence)} -> {len(augmented)}"
        )

    _write_truth(settings, output, events)
    _write_metadata(settings, output, events, stats, per_pattern)
    complete.touch()


def _write_truth(settings: dict, output: Path, events: list[Event]) -> None:
    genomes = settings["genomes"]
    with (output / "truth_events.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            list(TRUTH_FIELDS) + [f"{genome}_pos" for genome in genomes]
        )
        for event in events:
            carriers = [
                genome
                for genome, present in zip(genomes, event.pattern)
                if present
            ]
            writer.writerow(
                [
                    event.event_id,
                    event.te_id,
                    event.te_family,
                    event.te_group,
                    pattern_label(event.pattern),
                    sum(event.pattern),
                    ",".join(carriers),
                    event.anchor_id,
                    event.strand,
                    event.tsd_length,
                    event.inserted_length,
                ]
                + [
                    event.positions[genome]
                    if event.positions[genome] is not None
                    else "-"
                    for genome in genomes
                ]
            )

    with (output / "sharing_matrix.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("event_id", "genome", "present", "position"))
        for event in events:
            for index, genome in enumerate(genomes):
                present = event.pattern[index]
                writer.writerow(
                    (
                        event.event_id,
                        genome,
                        int(present),
                        event.positions[genome] if present else "-",
                    )
                )


def _write_metadata(
    settings: dict,
    output: Path,
    events: list[Event],
    stats: dict,
    per_pattern: int,
) -> None:
    chrom_key = settings["chrom"].lower()
    checksums = {}
    for genome in settings["genomes"]:
        for label, path in (
            ("input", output / "genomes" / chrom_key / f"{genome}.{chrom_key}.fa"),
            (
                "augmented",
                output / "genomes" / "augmented" / f"{genome}.{chrom_key}.te.fa",
            ),
        ):
            if path.is_file():
                checksums[f"{genome}.{label}"] = {
                    "path": str(path),
                    "sha256": sha256(path),
                }

    metadata = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "generator": "pipeline/make_pangenome_panel/05_build_pangenome_panel.py",
        "phase": 1,
        "settings": {
            key: value
            for key, value in settings.items()
            if key not in {"te_groups"}
        },
        "te_groups": [group["name"] for group in settings["te_groups"]],
        "anchor_acceptance": {
            key: stats[key]
            for key in (
                "candidates_sampled",
                "rejected_unlifted",
                "rejected_multimapping",
                "accepted",
                "needed",
            )
        },
        "events": {
            "total": len(events),
            "patterns": len(sharing_patterns(len(settings["genomes"]))),
            "events_per_pattern": per_pattern,
        },
        "checksums": checksums,
    }
    with (output / "run_metadata.json").open("w") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True, default=str)
        handle.write("\n")


if __name__ == "__main__":
    raise SystemExit(main())
