#!/usr/bin/env python3
"""Derive each event's coordinate in the TE-augmented genomes.

`truth_events.tsv` stores positions in the ORIGINAL genomes, because that is
the coordinate system anchors were chosen in. But the pangenome graph, and any
caller a user runs against the augmented FASTAs, works in AUGMENTED
coordinates -- where every insertion a genome carries shifts everything
downstream of it.

This emits `truth_events_augmented.tsv` with both systems side by side, so
graph output can be compared to truth without recomputing offsets by hand.

Derived rather than baked into the panel build: the genomes themselves are
correct, only the convenience column was missing.
"""

from __future__ import annotations

import argparse
import csv
import tomllib
from pathlib import Path


def augmented_positions(
    events: list[dict], genome: str, index: int
) -> dict[str, int]:
    """Map event_id -> augmented coordinate for one genome.

    Insertions are applied in ascending coordinate order, each contributing
    its TE plus its TSD, so an event's shift is the total inserted length of
    every event placed strictly before it.
    """
    carried = [
        event
        for event in events
        if event["pattern"][index] == "1" and event[f"{genome}_pos"] != "-"
    ]
    carried.sort(key=lambda event: int(event[f"{genome}_pos"]))

    positions, offset = {}, 0
    for event in carried:
        original = int(event[f"{genome}_pos"])
        positions[event["event_id"]] = original + offset
        offset += int(event["inserted_length"]) + int(event["tsd_length"])
    return positions


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    with args.config.open("rb") as handle:
        config = tomllib.load(handle)
    genomes = [genome["name"] for genome in config["genomes"]]

    with (args.output / "truth_events.tsv").open() as handle:
        events = list(csv.DictReader(handle, delimiter="\t"))

    augmented = {
        genome: augmented_positions(events, genome, index)
        for index, genome in enumerate(genomes)
    }

    fields = (
        ["event_id", "te_id", "te_family", "te_group", "pattern", "genomes"]
        + [f"{genome}_pos" for genome in genomes]
        + [f"{genome}_aug_pos" for genome in genomes]
        + ["tsd_length", "inserted_length"]
    )
    out_path = args.output / "truth_events_augmented.tsv"
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for event in events:
            row = {
                "event_id": event["event_id"],
                "te_id": event["te_id"],
                "te_family": event["te_family"],
                "te_group": event["te_group"],
                "pattern": event["pattern"],
                "genomes": event["genomes"],
                "tsd_length": event["tsd_length"],
                "inserted_length": event["inserted_length"],
            }
            for genome in genomes:
                row[f"{genome}_pos"] = event[f"{genome}_pos"]
                row[f"{genome}_aug_pos"] = augmented[genome].get(
                    event["event_id"], "-"
                )
            writer.writerow(row)

    for genome in genomes:
        shifts = [
            augmented[genome][event["event_id"]] - int(event[f"{genome}_pos"])
            for event in events
            if event["event_id"] in augmented[genome]
        ]
        print(
            f"{genome}\tevents={len(shifts)}\t"
            f"max_shift={max(shifts) if shifts else 0}"
        )
    print(f"Written: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
