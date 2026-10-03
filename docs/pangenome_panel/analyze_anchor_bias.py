#!/usr/bin/env python3
"""Quantify how unrepresentative the accepted anchor set is.

The orthology filter keeps only sites that lift 1:1 into every genome, which
systematically discards non-syntenic and repeat-dense regions -- exactly where
real TE insertions concentrate. This script measures that bias so the panel's
optimism can be reported as a number rather than asserted as a caveat.

Compares accepted vs rejected anchors on:
  - distance to the nearest background TE
  - background TE density in a window
  - fraction of the window masked by TEs
  - position along the chromosome (centromere proximity)
"""

from __future__ import annotations

import argparse
import bisect
import csv
import statistics
import tomllib
from pathlib import Path

WINDOW = 10_000


def load_background(path: Path) -> tuple[list[int], list[int]]:
    """Return sorted (starts, ends) of background TE intervals."""
    intervals = []
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            intervals.append((int(fields[1]), int(fields[2])))
    intervals.sort()
    return [s for s, _ in intervals], [e for _, e in intervals]


def nearest_distance(starts: list[int], ends: list[int], position: int) -> int:
    """Distance to the nearest background TE; 0 when inside one."""
    index = bisect.bisect_right(starts, position)
    best = None
    if index:
        # Scan back a few: intervals can nest and overlap.
        for probe in range(max(0, index - 8), index):
            if starts[probe] <= position <= ends[probe]:
                return 0
            best = min(best or 10**9, position - ends[probe])
    if index < len(starts):
        best = min(best or 10**9, starts[index] - position)
    return max(best or 10**9, 0)


def window_stats(
    starts: list[int], ends: list[int], position: int, window: int
) -> tuple[int, float]:
    """Count of TEs overlapping a window, and the masked fraction of it."""
    low, high = position - window, position + window
    index = bisect.bisect_left(starts, low)
    # Back up: an interval starting before `low` may still overlap it.
    index = max(0, index - 50)
    count, masked = 0, 0
    for probe in range(index, len(starts)):
        if starts[probe] > high:
            break
        overlap = min(ends[probe], high) - max(starts[probe], low)
        if overlap > 0:
            count += 1
            masked += overlap
    return count, masked / (2 * window)


def summarize(label: str, values: list[float]) -> dict:
    if not values:
        return {"set": label, "n": 0}
    ordered = sorted(values)
    return {
        "set": label,
        "n": len(values),
        "mean": statistics.mean(values),
        "median": statistics.median(ordered),
        "p25": ordered[len(ordered) // 4],
        "p75": ordered[3 * len(ordered) // 4],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    with args.config.open("rb") as handle:
        config = tomllib.load(handle)
    genomes = [g["name"] for g in config["genomes"]]
    reference = next(g["name"] for g in config["genomes"] if g.get("reference"))
    targets = [g for g in genomes if g != reference]

    ortho = args.output / "orthology"
    chrom = config["panel"]["chrom"]
    ref_record = f"{reference}#1#{chrom}"

    candidates = {}
    with (ortho / "candidate_anchors.bed").open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            candidates[f"{ref_record}_{fields[1]}_{fields[2]}"] = int(fields[1])

    lifted_counts: dict[str, int] = dict.fromkeys(candidates, 0)
    for target in targets:
        seen: dict[str, int] = {}
        with (ortho / f"lifted_{target}.bed").open() as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 4:
                    continue
                seen[fields[3]] = seen.get(fields[3], 0) + 1
        for key in candidates:
            if seen.get(key) == 1:
                lifted_counts[key] += 1

    accepted = [
        position
        for key, position in candidates.items()
        if lifted_counts[key] == len(targets)
    ]
    rejected = [
        position
        for key, position in candidates.items()
        if lifted_counts[key] != len(targets)
    ]

    starts, ends = load_background(
        args.output / "annotation" / f"{reference}.background_te.bed"
    )

    rows = []
    metrics: dict[str, dict[str, list[float]]] = {
        "distance_to_nearest_te": {"accepted": [], "rejected": []},
        "te_count_20kb_window": {"accepted": [], "rejected": []},
        "te_masked_fraction": {"accepted": [], "rejected": []},
        "position_mb": {"accepted": [], "rejected": []},
    }

    for label, positions in (("accepted", accepted), ("rejected", rejected)):
        for position in positions:
            distance = nearest_distance(starts, ends, position)
            count, masked = window_stats(starts, ends, position, WINDOW)
            metrics["distance_to_nearest_te"][label].append(distance)
            metrics["te_count_20kb_window"][label].append(count)
            metrics["te_masked_fraction"][label].append(masked)
            metrics["position_mb"][label].append(position / 1e6)

    for metric, sets in metrics.items():
        for label in ("accepted", "rejected"):
            summary = summarize(label, sets[label])
            summary["metric"] = metric
            rows.append(summary)

    out_path = ortho / "anchor_bias.tsv"
    fields = ("metric", "set", "n", "mean", "median", "p25", "p75")
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fields})

    print(f"accepted\t{len(accepted)}")
    print(f"rejected\t{len(rejected)}")
    print()
    print(f"{'metric':<28}{'accepted':>14}{'rejected':>14}{'ratio':>10}")
    for metric in metrics:
        acc = statistics.median(metrics[metric]["accepted"])
        rej = statistics.median(metrics[metric]["rejected"])
        ratio = (acc / rej) if rej else float("nan")
        print(f"{metric:<28}{acc:>14.3f}{rej:>14.3f}{ratio:>10.2f}")
    print()
    print(f"Written: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
