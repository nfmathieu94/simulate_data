#!/usr/bin/env python3
"""Check that the pangenome graph recovers the known TE sharing patterns.

This closes the loop on the whole panel. We know, by construction, which of
the four genomes carries each of the 300 TE insertions. Minigraph-Cactus
deconstructs the graph into a VCF with one column per genome, so the graph's
own genotypes can be compared against that truth.

Two things are being tested at once:

  detection    does a variant of roughly the right size appear at the
               expected locus at all?
  attribution  do the per-genome genotypes match the sharing pattern?

Attribution is the one that matters for RelocaTE3's "which genome(s)" report:
a graph that finds every insertion but assigns it to the wrong genomes would
be useless for the feature this panel exists to test.

Caveats, stated rather than hidden: the Cactus VCF is normalized and
`--vcfbub` flattens nested sites, so a TE inside another variant may be
reported differently than inserted. Matching is therefore positional with a
tolerance rather than exact, and unmatched events are reported separately
from mis-attributed ones.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import importlib.util
import sys
import tomllib
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent

MATCH_WINDOW = 200      # bp around the expected locus
LENGTH_TOLERANCE = 0.2  # fractional length agreement


def _load_builder():
    spec = importlib.util.spec_from_file_location(
        "build_pangenome_panel", HERE / "05_build_pangenome_panel.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


builder = _load_builder()


def open_maybe_gzip(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def parse_vcf(path: Path) -> tuple[list[str], list[dict]]:
    """Return (sample names, variant records) from a Cactus-deconstructed VCF."""
    samples: list[str] = []
    records: list[dict] = []
    with open_maybe_gzip(path) as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = fields[9:]
                continue
            if len(fields) < 10:
                continue
            ref, alts = fields[3], fields[4].split(",")
            # Length change of the largest alt relative to ref: positive for
            # insertion-like, negative for deletion-like.
            deltas = [len(alt) - len(ref) for alt in alts]
            records.append(
                {
                    "chrom": fields[0],
                    "pos": int(fields[1]),
                    "ref_len": len(ref),
                    "deltas": deltas,
                    "max_abs_delta": max(abs(d) for d in deltas),
                    "genotypes": {
                        sample: fields[9 + index].split(":")[0]
                        for index, sample in enumerate(samples)
                    },
                }
            )
    return samples, records


def genotype_carries_alt(genotype: str) -> bool:
    """True when a genotype contains any non-reference allele."""
    if genotype in (".", "./.", ".|.", ""):
        return False
    alleles = genotype.replace("|", "/").split("/")
    return any(allele not in ("0", ".") for allele in alleles)


def index_by_position(records: list[dict]) -> tuple[list[int], list[dict]]:
    ordered = sorted(records, key=lambda record: record["pos"])
    return [record["pos"] for record in ordered], ordered


def find_match(
    positions: list[int],
    ordered: list[dict],
    target: int,
    expected_length: int,
) -> dict | None:
    """Nearest variant within the window whose size roughly agrees."""
    import bisect

    low = bisect.bisect_left(positions, target - MATCH_WINDOW)
    high = bisect.bisect_right(positions, target + MATCH_WINDOW)
    best, best_distance = None, None
    for record in ordered[low:high]:
        if expected_length:
            ratio = record["max_abs_delta"] / expected_length
            if not (1 - LENGTH_TOLERANCE) <= ratio <= (1 + LENGTH_TOLERANCE):
                continue
        distance = abs(record["pos"] - target)
        if best_distance is None or distance < best_distance:
            best, best_distance = record, distance
    return best


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--vcf", type=Path, help="Override the graph VCF path")
    args = parser.parse_args()

    settings = builder.load_settings(args.config)
    with args.config.open("rb") as handle:
        config = tomllib.load(handle)
    genomes = [genome["name"] for genome in config["genomes"]]
    reference = settings["reference"]

    vcf_path = args.vcf
    if vcf_path is None:
        candidates = sorted((args.output / "graph").glob("*.vcf.gz"))
        if not candidates:
            raise SystemExit(
                f"No VCF found in {args.output / 'graph'}; pass --vcf"
            )
        vcf_path = candidates[0]
    print(f"VCF\t{vcf_path}")

    samples, records = parse_vcf(vcf_path)
    print(f"vcf_samples\t{','.join(samples)}")
    print(f"vcf_records\t{len(records)}")

    positions, ordered = index_by_position(records)

    # Augmented reference coordinates for every event, including those the
    # reference itself does not carry (their anchor still exists in it).
    anchors, _stats = builder.build_anchors(settings, args.output)
    with (args.output / "truth_events.tsv").open() as handle:
        events = list(csv.DictReader(handle, delimiter="\t"))

    ref_index = genomes.index(reference)
    carried = sorted(
        (event for event in events if event["pattern"][ref_index] == "1"),
        key=lambda event: int(event[f"{reference}_pos"]),
    )
    # shift(p) = total sequence inserted at positions strictly before p.
    shift_points = [int(event[f"{reference}_pos"]) for event in carried]
    cumulative = [0]
    for event in carried:
        cumulative.append(
            cumulative[-1]
            + int(event["inserted_length"])
            + int(event["tsd_length"])
        )

    def augmented_reference_position(original: int) -> int:
        import bisect

        return original + cumulative[bisect.bisect_left(shift_points, original)]

    rows = []
    stats = Counter()
    for event in events:
        original = event.get(f"{reference}_pos", "-")
        if original != "-":
            reference_original = int(original)
        else:
            anchor = anchors.get(event["anchor_id"])
            if anchor is None:
                stats["no_anchor"] += 1
                continue
            reference_original = anchor["positions"][reference]

        target = augmented_reference_position(reference_original)
        expected_length = int(event["inserted_length"])
        match = find_match(positions, ordered, target, expected_length)

        expected = {
            genome: event["pattern"][index] == "1"
            for index, genome in enumerate(genomes)
        }
        if match is None:
            stats["undetected"] += 1
            rows.append(
                {
                    "event_id": event["event_id"],
                    "pattern": event["pattern"],
                    "te_group": event["te_group"],
                    "detected": 0,
                    "attribution": "NA",
                    "vcf_pos": "-",
                }
            )
            continue

        stats["detected"] += 1
        observed = {
            genome: genotype_carries_alt(match["genotypes"].get(genome, "."))
            for genome in genomes
            if genome in match["genotypes"]
        }
        # The reference genome has no VCF column; its state is implied by
        # whether the TE sits on the reference path.
        comparable = {g: v for g, v in observed.items() if g != reference}
        agree = all(
            comparable[genome] == (not expected[genome])
            if expected[reference]
            else comparable[genome] == expected[genome]
            for genome in comparable
        )
        stats["attribution_ok" if agree else "attribution_mismatch"] += 1
        rows.append(
            {
                "event_id": event["event_id"],
                "pattern": event["pattern"],
                "te_group": event["te_group"],
                "detected": 1,
                "attribution": "match" if agree else "mismatch",
                "vcf_pos": match["pos"],
            }
        )

    out_path = args.output / "graph_validation.tsv"
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "event_id",
                "pattern",
                "te_group",
                "detected",
                "attribution",
                "vcf_pos",
            ),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    total = len(rows)
    print()
    print(f"events\t{total}")
    print(f"detected\t{stats['detected']}\t{stats['detected'] / max(total,1):.3f}")
    print(f"undetected\t{stats['undetected']}")
    print(f"attribution_ok\t{stats['attribution_ok']}")
    print(f"attribution_mismatch\t{stats['attribution_mismatch']}")
    print()
    print("detection by pattern cardinality:")
    by_card = Counter()
    seen = Counter()
    for row in rows:
        card = sum(1 for c in row["pattern"] if c == "1")
        seen[card] += 1
        by_card[card] += row["detected"]
    for card in sorted(seen):
        print(f"  {card} genome(s)\t{by_card[card]}/{seen[card]}")
    print()
    print(f"Written: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
