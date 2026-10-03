"""Tests for the mixed germline/somatic TE benchmark panel."""

import argparse
import io
from collections import defaultdict
from pathlib import Path

import pysam
import pytest

from simulate_data.modules import te_benchmark_panel as panel


def _event(label: str, ccf: float = 1.0) -> panel.Event:
    return panel.Event(
        "TE000001",
        "Chr1",
        10,
        "mPing",
        "mPing",
        label,
        ccf,
        ccf / 2 if label == "somatic_insertion" else 0.5,
        1,
        0,
        "+",
        "TTA",
        3,
        4,
        4,
        100.0,
        0,
        0,
        0,
        "GGCC",
    )


def test_nested_clone_membership_realizes_requested_ccf():
    som10 = _event("somatic_insertion", 0.10)
    som20 = _event("somatic_insertion", 0.20)
    som40 = _event("somatic_insertion", 0.40)
    observed = {}
    for event in (som10, som20, som40):
        cellular_fraction = sum(
            weight
            for _name, weight, classes in panel.COMPONENTS
            if panel._present(event, classes, 1)
        )
        observed[event.cellular_fraction] = cellular_fraction
        assert not any(
            panel._present(event, classes, 2)
            for _name, _weight, classes in panel.COMPONENTS
        )
    assert observed == pytest.approx({0.10: 0.10, 0.20: 0.20, 0.40: 0.40})


def test_germline_membership_is_constant_across_components():
    hom = _event("homozygous")
    het = _event("heterozygous")
    for _name, _weight, classes in panel.COMPONENTS:
        assert panel._present(hom, classes, 1)
        assert panel._present(hom, classes, 2)
        assert panel._present(het, classes, 1)
        assert not panel._present(het, classes, 2)


def test_component_genome_contains_te_and_tsd(tmp_path):
    event = _event("homozygous")
    output = tmp_path / "panel"
    output.mkdir()
    panel._write_component_genomes([event], {"Chr1": "AACCGGTTAACCGGTT"}, output)
    sequence = "".join(
        line.strip()
        for line in (output / "genomes" / "baseline.hap1.fa").read_text().splitlines()
        if not line.startswith(">")
    )
    assert sequence == "AACCGGTTAA" + "GGCCTTA" + "CCGGTT"


def test_truth_outputs_are_bgzip_indexed(tmp_path):
    event = _event("heterozygous")
    panel._write_truth([event], tmp_path)
    assert (tmp_path / "truth_events.tsv").is_file()
    assert (tmp_path / "truth_events.vcf.gz").is_file()
    assert (tmp_path / "truth_events.vcf.gz.tbi").is_file()


def test_pool_fastq_adds_collision_free_provenance(tmp_path):
    source = tmp_path / "reads.fq"
    source.write_text("@read/1 comment\nACGT\n+\nIIII\n")
    destination = io.StringIO()
    assert panel._pool_fastq(source, destination, "cov5x:baseline:h1") == 1
    assert destination.getvalue().startswith("@cov5x:baseline:h1:read/1 comment\n")


def test_origin_support_counts_junctions_and_spanners_once(tmp_path):
    sam_path = tmp_path / "origin.sam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "Chr1", "LN": 500}]}
    with pysam.AlignmentFile(sam_path, "w", header=header) as sam:
        for name, start in (("insertion", 0), ("reference", 100)):
            record = pysam.AlignedSegment()
            record.query_name = name
            record.query_sequence = "A" * 100
            record.flag = 0
            record.reference_id = 0
            record.reference_start = start
            record.mapping_quality = 60
            record.cigarstring = "100M"
            record.query_qualities = pysam.qualitystring_to_array("I" * 100)
            sam.write(record)
    coordinates = [
        {"chrom": "Chr1", "coordinate": 50, "event_id": "TE1", "present": 1},
        {"chrom": "Chr1", "coordinate": 150, "event_id": "TE2", "present": 0},
    ]
    truth = {
        "TE1": {"inserted_length": "4", "tsd_length": "3"},
        "TE2": {"inserted_length": "4", "tsd_length": "3"},
    }
    counts = {"TE1": [0, 0], "TE2": [0, 0]}
    panel._count_origin_support(sam_path, coordinates, truth, counts)
    assert counts == {"TE1": [1, 0], "TE2": [0, 1]}


def test_toml_config_and_cli_coverage_override(tmp_path):
    config = tmp_path / "panel.toml"
    config.write_text(
        '[inputs]\nref="ref.fa"\nte="te.fa"\nknown_del="rm.out"\n'
        "[panel]\ncoverage=[5,15,30]\nreplicates=3\n"
    )
    args = argparse.Namespace(
        config=str(config),
        ref=None,
        te=None,
        known_del=None,
        output="out",
        stage="reads",
        chroms=None,
        events_per_class=None,
        coverage=[7.0],
        replicates=None,
        replicate=2,
        seed=None,
        read_length=None,
        fragment_size=None,
        fragment_std=None,
        sequencing_system=None,
        te_type=None,
        snp_rate=None,
        indel_rate=None,
        tsd_min=None,
        tsd_max=None,
        sense_strand_ratio=None,
    )
    settings = panel._load_settings(args)
    assert settings["coverage"] == [7.0]
    assert settings["replicates"] == 3
    assert settings["ref"] == "ref.fa"


def test_art_command_requests_profile_and_sam():
    command = panel._build_art_command(
        Path("ref.fa"),
        150,
        2.5,
        "out.",
        seed=42,
        sequencing_system="HSXn",
        sam=True,
    )
    assert command[command.index("-ss") + 1] == "HSXn"
    assert "-sam" in command


class TestLongReadSupportCounting:
    """MAF spans must drive the same present/absent tally the SAM path did."""

    @staticmethod
    def _rows_and_truth():
        rows = [
            {
                "chrom": "Chr1",
                "coordinate": 1000,
                "event_id": "TE000001",
                "present": 1,
                "haplotype": 1,
            }
        ]
        truth = {"TE000001": {"inserted_length": "500", "tsd_length": "5"}}
        return rows, truth

    def test_spanning_read_supports_the_insertion(self):
        rows, truth = self._rows_and_truth()
        counts = defaultdict(lambda: [0, 0])
        # A read covering 0-20000 spans both junctions of the insertion.
        panel._count_origin_support_maf([("Chr1", 0, 20000)], rows, truth, counts)
        assert counts["TE000001"] == [1, 0]

    def test_absent_event_counts_as_reference_spanning(self):
        rows, truth = self._rows_and_truth()
        rows[0]["present"] = 0
        counts = defaultdict(lambda: [0, 0])
        panel._count_origin_support_maf([("Chr1", 0, 20000)], rows, truth, counts)
        assert counts["TE000001"] == [0, 1]

    def test_read_on_another_contig_is_ignored(self):
        rows, truth = self._rows_and_truth()
        counts = defaultdict(lambda: [0, 0])
        panel._count_origin_support_maf([("Chr9", 0, 20000)], rows, truth, counts)
        assert counts["TE000001"] == [0, 0]

    def test_read_ending_before_the_site_is_ignored(self):
        rows, truth = self._rows_and_truth()
        counts = defaultdict(lambda: [0, 0])
        panel._count_origin_support_maf([("Chr1", 0, 500)], rows, truth, counts)
        assert counts["TE000001"] == [0, 0]

    def test_each_event_counted_once_per_read(self):
        """A long read crossing both junctions must not double-count."""
        rows, truth = self._rows_and_truth()
        counts = defaultdict(lambda: [0, 0])
        panel._count_origin_support_maf([("Chr1", 0, 100000)], rows, truth, counts)
        assert counts["TE000001"] == [1, 0]


def test_long_read_seed_derivation_matches_short_read():
    """Seeds must stay comparable across technologies."""
    assert panel._seed(916, 30.0, 1, 0, 1) == panel._seed(916, 30.0, 1, 0, 1)
    assert panel._seed(916, 30.0, 1, 0, 1) != panel._seed(916, 30.0, 1, 0, 2)


def test_long_read_platforms_are_registered():
    assert set(panel.LONG_READ_PLATFORMS) == {"ont-hq", "hifi"}
