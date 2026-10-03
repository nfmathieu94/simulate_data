#!/usr/bin/env python3
"""Build the Phase 2 sample genome whose reads are what a caller actually sees.

The multi-genome panel supplies the *references*. This builds the *sample*:
a genome derived from the reference background carrying three deliberately
different classes of insertion, so both halves of the RelocaTE3 report have
ground truth.

    reference     an insertion the sample shares with >=1 panel genome.
                  Expected call: "reference", attributed to exactly the
                  panel genomes carrying that event.
    non_reference an insertion at a locus empty in every panel genome.
                  Expected call: "non-reference".
    absent        a panel event the sample does NOT carry. Expected: no call
                  at all. This is the over-calling control, and without it a
                  caller that reports every panel event scores perfectly.

Nothing here is rice-specific: genome names, insertion counts and the TE
library all come from config, so the logic can migrate into simulate_data as
a general pangenome-simulation module.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import random
import sys
import tomllib
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

from Bio.Seq import Seq

HERE = Path(__file__).resolve().parent


def _load_builder():
    spec = importlib.util.spec_from_file_location(
        "build_pangenome_panel", HERE / "05_build_pangenome_panel.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


builder = _load_builder()

SAMPLE_TRUTH_FIELDS = (
    "event_id",
    "sample_class",
    "expected_call",
    "expected_genomes",
    "te_id",
    "te_family",
    "te_group",
    "panel_pattern",
    "sample_position",
    "reference_position",
    "strand",
    "tsd",
    "tsd_length",
    "inserted_length",
)


@dataclass
class SampleEvent:
    event_id: str
    sample_class: str
    expected_call: str
    expected_genomes: str
    te_id: str
    te_family: str
    te_group: str
    panel_pattern: str
    reference_position: int
    strand: str
    tsd: str
    tsd_length: int
    inserted_length: int
    sequence: str
    sample_position: int | None = None


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser


def load_sample_settings(config_path: Path) -> dict:
    settings = builder.load_settings(config_path)
    with config_path.open("rb") as handle:
        config = tomllib.load(handle)
    sample = config.get("sample", {})
    settings.update(
        {
            "sample_name": sample.get("name", "SampleA"),
            "n_reference": sample.get("n_reference", 150),
            "n_non_reference": sample.get("n_non_reference", 100),
            "sample_seed": sample.get("seed", settings["seed"] + 100),
        }
    )
    return settings


def load_panel_events(path: Path) -> list[dict]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def choose_reference_events(
    panel_events: list[dict], n_reference: int, rng: random.Random
) -> tuple[list[dict], list[dict]]:
    """Split panel events into those the sample carries and those it does not.

    Sampling is stratified by sharing pattern so every pattern contributes
    both carried and absent events; otherwise whole patterns could vanish and
    attribution would go untested for them.
    """
    by_pattern: dict[str, list[dict]] = {}
    for event in panel_events:
        by_pattern.setdefault(event["pattern"], []).append(event)

    patterns = sorted(by_pattern)
    per_pattern = n_reference // len(patterns)
    remainder = n_reference - per_pattern * len(patterns)

    carried, absent = [], []
    for index, pattern in enumerate(patterns):
        events = by_pattern[pattern][:]
        rng.shuffle(events)
        take = per_pattern + (1 if index < remainder else 0)
        carried.extend(events[:take])
        absent.extend(events[take:])
    return carried, absent


def used_anchor_ids(panel_events: list[dict]) -> set[str]:
    return {event["anchor_id"] for event in panel_events}


def reference_position_for(event: dict, anchors: dict, reference: str) -> int:
    """Reference-background coordinate for a panel event.

    truth_events.tsv records "-" for a genome that does not carry the event,
    so an Azucena-only event has no Nipponbare position there. The anchor it
    was placed on still has a coordinate in every genome -- that is the whole
    point of the orthology filter -- so the position comes from the anchor
    map. This is what lets the sample carry a TE that is "reference" by virtue
    of another genome while being absent from its own background.
    """
    recorded = event.get(f"{reference}_pos", "-")
    if recorded not in ("-", ""):
        return int(recorded)
    anchor = anchors.get(event["anchor_id"])
    if anchor is None:
        raise KeyError(
            f"Anchor {event['anchor_id']} for {event['event_id']} is not in "
            "the accepted anchor set"
        )
    return anchor["positions"][reference]


def build_sample_events(
    settings: dict,
    panel_events: list[dict],
    anchors: dict,
    reference_sequence: str,
) -> list[SampleEvent]:
    rng = random.Random(settings["sample_seed"])
    reference = settings["reference"]
    genomes = settings["genomes"]

    carried, absent = choose_reference_events(
        panel_events, settings["n_reference"], rng
    )

    events: list[SampleEvent] = []

    # 1. Reference-class: shared with >=1 panel genome, at the same locus.
    for event in carried:
        carriers = event["genomes"]
        events.append(
            SampleEvent(
                event_id=event["event_id"],
                sample_class="reference",
                expected_call="reference",
                expected_genomes=carriers,
                te_id=event["te_id"],
                te_family=event["te_family"],
                te_group=event["te_group"],
                panel_pattern=event["pattern"],
                reference_position=reference_position_for(
                    event, anchors, reference
                ),
                strand=event["strand"],
                tsd="",
                tsd_length=int(event["tsd_length"]),
                inserted_length=int(event["inserted_length"]),
                sequence="",
            )
        )

    # 2. Absent: in the panel but not the sample. No sequence is inserted;
    #    these exist purely so over-calling is measurable.
    for event in absent:
        events.append(
            SampleEvent(
                event_id=event["event_id"],
                sample_class="absent",
                expected_call="none",
                expected_genomes=event["genomes"],
                te_id=event["te_id"],
                te_family=event["te_family"],
                te_group=event["te_group"],
                panel_pattern=event["pattern"],
                reference_position=reference_position_for(
                    event, anchors, reference
                ),
                strand=event["strand"],
                tsd="",
                tsd_length=int(event["tsd_length"]),
                inserted_length=int(event["inserted_length"]),
                sequence="",
            )
        )

    # 3. Non-reference: novel insertions at accepted anchors the panel never
    #    used, so the locus is empty in all four panel genomes by construction.
    library = builder.read_te_library(settings["te_library"])
    by_group = {
        group["name"]: builder.candidates_for_group(
            library, group, settings["selection"]
        )
        for group in settings["te_groups"]
    }
    spare = sorted(set(anchors) - used_anchor_ids(panel_events))
    rng.shuffle(spare)
    if len(spare) < settings["n_non_reference"]:
        raise RuntimeError(
            f"Only {len(spare)} unused anchors available; "
            f"{settings['n_non_reference']} non-reference events requested"
        )

    for index in range(settings["n_non_reference"]):
        anchor_id = spare[index]
        anchor = anchors[anchor_id]
        group = settings["te_groups"][index % len(settings["te_groups"])]
        candidate = rng.choice(by_group[group["name"]])
        tsd_length = rng.randint(group["tsd_min"], group["tsd_max"])
        strand = (
            "+"
            if rng.random() < settings["selection"]["sense_strand_ratio"]
            else "-"
        )
        sequence = (
            candidate.sequence
            if strand == "+"
            else str(Seq(candidate.sequence).reverse_complement())
        )
        events.append(
            SampleEvent(
                event_id=f"SMPL{index + 1:06d}",
                sample_class="non_reference",
                expected_call="non_reference",
                expected_genomes="",
                te_id=candidate.te_id,
                te_family=candidate.te_family,
                te_group=group["name"],
                panel_pattern="0" * len(genomes),
                reference_position=anchor["positions"][reference],
                strand=strand,
                tsd="",
                tsd_length=tsd_length,
                inserted_length=len(sequence),
                sequence=sequence,
            )
        )

    # Fill in TE sequence and TSD for everything the sample actually carries.
    panel_sequences = {
        event["event_id"]: event for event in panel_events
    }
    for event in events:
        if event.sample_class == "absent":
            continue
        position = event.reference_position
        event.tsd = (
            reference_sequence[position - event.tsd_length : position]
            if event.tsd_length
            else ""
        )
        if event.sample_class == "reference":
            # Rebuild the exact sequence used in the panel so the sample and
            # the panel genomes carry an identical element at this locus.
            record = panel_sequences[event.event_id]
            library_hit = next(
                (c for c in library if c.te_id == record["te_id"]), None
            )
            if library_hit is None:
                raise RuntimeError(f"TE not found in library: {record['te_id']}")
            event.sequence = (
                library_hit.sequence
                if event.strand == "+"
                else str(Seq(library_hit.sequence).reverse_complement())
            )
            event.inserted_length = len(event.sequence)
    return events


def write_sample_genome(
    settings: dict,
    output: Path,
    events: list[SampleEvent],
    record_name: str,
    reference_sequence: str,
) -> tuple[str, int]:
    """Insert carried events and record each one's sample coordinate."""
    carried = [e for e in events if e.sample_class != "absent"]
    insertions = [
        (e.reference_position, e.sequence, e.tsd) for e in carried
    ]
    augmented = builder.insert_all(reference_sequence, insertions)

    # Coordinates shift by everything inserted at or before a site.
    ordered = sorted(carried, key=lambda e: e.reference_position)
    offset = 0
    for event in ordered:
        event.sample_position = event.reference_position + offset
        offset += len(event.sequence) + len(event.tsd)

    sample_dir = output / "sample" / settings["sample_name"]
    sample_dir.mkdir(parents=True, exist_ok=True)
    path = sample_dir / f"{settings['sample_name']}.chr1.fa"
    builder.write_fasta(path, record_name, augmented)
    return str(path), len(augmented)


def write_sample_truth(
    settings: dict, output: Path, events: list[SampleEvent]
) -> None:
    sample_dir = output / "sample" / settings["sample_name"]
    with (sample_dir / "sample_truth.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=SAMPLE_TRUTH_FIELDS, delimiter="\t"
        )
        writer.writeheader()
        for event in sorted(events, key=lambda e: e.reference_position):
            writer.writerow(
                {
                    "event_id": event.event_id,
                    "sample_class": event.sample_class,
                    "expected_call": event.expected_call,
                    "expected_genomes": event.expected_genomes or "-",
                    "te_id": event.te_id,
                    "te_family": event.te_family,
                    "te_group": event.te_group,
                    "panel_pattern": event.panel_pattern,
                    "sample_position": (
                        event.sample_position
                        if event.sample_position is not None
                        else "-"
                    ),
                    "reference_position": event.reference_position,
                    "strand": event.strand,
                    "tsd": event.tsd or "NONE",
                    "tsd_length": event.tsd_length,
                    "inserted_length": event.inserted_length,
                }
            )


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    settings = load_sample_settings(args.config)
    output = args.output

    sample_dir = output / "sample" / settings["sample_name"]
    if (sample_dir / ".complete").exists():
        print(f"Sample already complete: {sample_dir}")
        return 0

    if not (output / ".panel.complete").exists():
        raise FileNotFoundError(
            "Panel is incomplete; run the Phase 1 panel build first"
        )

    chrom_key = settings["chrom"].lower()
    reference = settings["reference"]
    record_name, reference_sequence = builder.read_fasta_single(
        output / "genomes" / chrom_key / f"{reference}.{chrom_key}.fa"
    )
    sample_record = record_name.replace(reference, settings["sample_name"], 1)

    panel_events = load_panel_events(output / "truth_events.tsv")
    anchors, _stats = builder.build_anchors(settings, output)

    events = build_sample_events(
        settings, panel_events, anchors, reference_sequence
    )
    path, length = write_sample_genome(
        settings, output, events, sample_record, reference_sequence
    )
    write_sample_truth(settings, output, events)

    counts = {
        label: sum(1 for e in events if e.sample_class == label)
        for label in ("reference", "non_reference", "absent")
    }
    metadata = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "generator": "pipeline/make_pangenome_panel/06_build_sample_genome.py",
        "phase": 2,
        "sample_name": settings["sample_name"],
        "record": sample_record,
        "derived_from": reference,
        "genome_path": path,
        "genome_length": length,
        "reference_length": len(reference_sequence),
        "event_counts": counts,
        "expected_calls": {
            "reference": counts["reference"],
            "non_reference": counts["non_reference"],
            "none": counts["absent"],
        },
    }
    with (sample_dir / "run_metadata.json").open("w") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
        handle.write("\n")

    (sample_dir / ".complete").touch()

    for label, count in counts.items():
        print(f"{label}\t{count}")
    print(f"genome\t{len(reference_sequence)} -> {length}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
