"""End-to-end long-read panel simulation against real PBSIM3 and ccs.

These tests actually run the external tools, so they are marked slow.
"""

import csv
import gzip
import random
from pathlib import Path

import pytest

from simulate_data.modules import te_benchmark_panel as panel

pytestmark = pytest.mark.slow

COMPONENT_NAMES = [component for component, _weight, _classes in panel.COMPONENTS]


def _write_fasta(path: Path, name: str, sequence: str) -> None:
    with path.open("w") as handle:
        handle.write(f">{name}\n")
        for start in range(0, len(sequence), 80):
            handle.write(sequence[start : start + 80] + "\n")


@pytest.fixture
def catalog(tmp_path):
    """A minimal but structurally complete catalog: 1 contig, 2 events."""
    rng = random.Random(1)
    reference = "".join(rng.choice("ACGT") for _ in range(120_000))
    insertion = "".join(rng.choice("ACGT") for _ in range(600))

    output = tmp_path / "panel"
    genomes = output / "genomes"
    genomes.mkdir(parents=True)

    events = [
        {"event_id": "TE000001", "position": 30_000, "class": "homozygous"},
        {"event_id": "TE000002", "position": 80_000, "class": "somatic_insertion"},
    ]

    with (output / "truth_events.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=panel.TRUTH_FIELDS, delimiter="\t")
        writer.writeheader()
        for event in events:
            writer.writerow(
                {
                    "event_id": event["event_id"],
                    "chrom": "Chr1",
                    "position": event["position"],
                    "te_id": "mPing",
                    "te_family": "mPing",
                    "biological_class": event["class"],
                    "cellular_fraction": 1.0,
                    "expected_vaf": 1.0,
                    "hap1": 1,
                    "hap2": 1,
                    "strand": "+",
                    "tsd": "TTA",
                    "tsd_length": 3,
                    "consensus_length": len(insertion),
                    "inserted_length": len(insertion),
                    "sequence_identity": 100.0,
                    "snp_count": 0,
                    "indel_count": 0,
                    "length_change": 0,
                }
            )

    rows = []
    for component, weight, _classes in panel.COMPONENTS:
        for hap in (1, 2):
            # Insertions are present in every component so the tiny fixture
            # produces support at low depth; mixture logic is unit-tested.
            chunks, cursor, offset = [], 0, 0
            for event in events:
                position = event["position"]
                chunks.append(reference[cursor:position])
                chunks.append(insertion + "TTA")
                cursor = position
                rows.append(
                    (
                        component,
                        weight,
                        hap,
                        event["event_id"],
                        "Chr1",
                        position + offset,
                        1,
                    )
                )
                offset += len(insertion) + 3
            chunks.append(reference[cursor:])
            _write_fasta(genomes / f"{component}.hap{hap}.fa", "Chr1", "".join(chunks))

    _write_fasta(genomes / "reference.fa", "Chr1", reference)

    with (output / "component_coordinates.tsv").open("w", newline="") as handle:
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
        writer.writerows(rows)

    (output / ".catalog.complete").touch()
    return output


def _assert_no_uncompressed_fastq(root: Path) -> None:
    leftovers = [
        str(path)
        for path in root.rglob("*")
        if path.suffix in {".fastq", ".fq"} and path.is_file()
    ]
    assert not leftovers, f"uncompressed FASTQ must never be written: {leftovers}"


def _read_count(path: Path) -> int:
    with gzip.open(path, "rt") as handle:
        return sum(1 for index, _line in enumerate(handle) if index % 4 == 0)


@pytest.mark.parametrize("platform", ["ont-hq", "hifi"])
def test_long_read_panel_end_to_end(catalog, platform):
    settings = {"seed": 916, "min_ccs_yield": 0.3}

    panel._simulate_long_reads(settings, catalog, 2.0, 1, platform)

    sample_dir = catalog / "reads" / platform / "cov2x_rep1"
    control_dir = catalog / "reads" / platform / "cov2x_rep1_reference_control"

    reads = sample_dir / "reads.fastq.gz"
    assert reads.is_file() and reads.stat().st_size > 0
    assert _read_count(reads) > 0
    assert _read_count(control_dir / "reads.fastq.gz") > 0

    # Compression contract.
    _assert_no_uncompressed_fastq(catalog)

    # One manifest row per component/haplotype pair.
    with (sample_dir / "component_manifest.tsv").open() as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    assert len(manifest) == len(panel.COMPONENTS) * 2
    assert {row["component"] for row in manifest} == set(COMPONENT_NAMES)
    assert sum(int(row["reads"]) for row in manifest) == _read_count(reads)

    # HiFi records its ccs yield; ONT has none to record.
    if platform == "hifi":
        assert all(row["ccs_yield"] != "NA" for row in manifest)
    else:
        assert all(row["ccs_yield"] == "NA" for row in manifest)

    # QC table covers every truth event.
    with (sample_dir / "observed_support.tsv").open() as handle:
        support = list(csv.DictReader(handle, delimiter="\t"))
    assert {row["event_id"] for row in support} == {"TE000001", "TE000002"}

    assert (sample_dir / ".complete").exists()
    assert (control_dir / ".complete").exists()


def test_rerun_is_a_noop(catalog):
    settings = {"seed": 916, "min_ccs_yield": 0.3}
    panel._simulate_long_reads(settings, catalog, 2.0, 1, "ont-hq")

    reads = catalog / "reads" / "ont-hq" / "cov2x_rep1" / "reads.fastq.gz"
    before = reads.stat().st_mtime_ns

    panel._simulate_long_reads(settings, catalog, 2.0, 1, "ont-hq")

    assert reads.stat().st_mtime_ns == before


def test_unknown_platform_is_rejected(catalog):
    with pytest.raises(ValueError, match="Unknown long-read platform"):
        panel._simulate_long_reads({"seed": 916}, catalog, 2.0, 1, "nanopore")


def test_incomplete_catalog_is_refused(tmp_path):
    with pytest.raises(FileNotFoundError, match="Catalog is incomplete"):
        panel._simulate_long_reads({"seed": 916}, tmp_path, 2.0, 1, "ont-hq")
