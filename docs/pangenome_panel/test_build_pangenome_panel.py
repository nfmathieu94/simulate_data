"""Focused tests for the multi-genome / pangenome TE panel builder."""

import collections
import importlib.util
import random
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent


def _load(module_name: str, filename: str):
    spec = importlib.util.spec_from_file_location(module_name, HERE / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


prepare = _load("prepare_chromosomes", "prepare_chromosomes.py")


ASSEMBLY_REPORT = (
    "# Assembly name: AzucenaRS1\n"
    "# Organism name: Oryza sativa\n"
    "# Sequence-Name\tSequence-Role\tAssigned-Molecule\t"
    "Assigned-Molecule-Location/Type\tGenBank-Accn\tRelationship\t"
    "RefSeq-Accn\tAssembly-Unit\tSequence-Length\tUCSC-style-name\n"
    "1\tassembled-molecule\t1\tChromosome\tCM020633.1\t<>\tna\t"
    "Primary Assembly\t43000000\tna\n"
    "2\tassembled-molecule\t2\tChromosome\tCM020634.1\t<>\tna\t"
    "Primary Assembly\t36000000\tna\n"
)


class TestChromosomeResolution(unittest.TestCase):
    """Chromosome identity comes from the assembly report, not header text."""

    def _report(self, directory: Path) -> Path:
        path = directory / "report.txt"
        path.write_text(ASSEMBLY_REPORT)
        return path

    def test_resolves_accession_for_chromosome(self):
        with tempfile.TemporaryDirectory() as tmp:
            report = self._report(Path(tmp))
            self.assertEqual(
                prepare.accession_for_chromosome(report, "1"), "CM020633.1"
            )

    def test_resolves_second_chromosome(self):
        with tempfile.TemporaryDirectory() as tmp:
            report = self._report(Path(tmp))
            self.assertEqual(
                prepare.accession_for_chromosome(report, "2"), "CM020634.1"
            )

    def test_missing_chromosome_raises(self):
        with tempfile.TemporaryDirectory() as tmp:
            report = self._report(Path(tmp))
            with self.assertRaises(ValueError):
                prepare.accession_for_chromosome(report, "13")

    def test_unplaced_scaffolds_are_ignored(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "report.txt"
            path.write_text(
                ASSEMBLY_REPORT
                + "scaf1\tunplaced-scaffold\tna\tna\tJAAA01.1\t<>\tna\t"
                "Primary Assembly\t1000\tna\n"
            )
            with self.assertRaises(ValueError):
                prepare.accession_for_chromosome(path, "na")


class TestPanSN(unittest.TestCase):
    """Graph tools expect sample#haplotype#contig naming."""

    def test_pansn_name(self):
        self.assertEqual(prepare.pansn_name("Azucena", "Chr1"), "Azucena#1#Chr1")

    def test_pansn_name_custom_haplotype(self):
        self.assertEqual(
            prepare.pansn_name("N22", "Chr1", haplotype=2), "N22#2#Chr1"
        )

    def test_sample_name_may_not_contain_separator(self):
        with self.assertRaises(ValueError):
            prepare.pansn_name("bad#name", "Chr1")


class TestSharingPatterns(unittest.TestCase):
    """All 15 non-empty subsets of four genomes."""

    def setUp(self):
        self.builder = _load("build_pangenome_panel", "05_build_pangenome_panel.py")

    def test_fifteen_non_empty_patterns(self):
        patterns = self.builder.sharing_patterns(4)
        self.assertEqual(len(patterns), 15)
        self.assertNotIn((False, False, False, False), patterns)

    def test_patterns_cover_every_cardinality(self):
        counts = collections.Counter(
            sum(pattern) for pattern in self.builder.sharing_patterns(4)
        )
        self.assertEqual(counts, {1: 4, 2: 6, 3: 4, 4: 1})

    def test_each_genome_appears_in_eight_patterns(self):
        patterns = self.builder.sharing_patterns(4)
        for index in range(4):
            self.assertEqual(sum(1 for p in patterns if p[index]), 8)

    def test_patterns_are_unique(self):
        patterns = self.builder.sharing_patterns(4)
        self.assertEqual(len(set(patterns)), len(patterns))

    def test_pattern_label_is_a_bitmask(self):
        self.assertEqual(self.builder.pattern_label((True, False, True, False)), "1010")


class TestAnchorFiltering(unittest.TestCase):
    """An anchor is only usable if it lifts 1:1 into every other genome."""

    def setUp(self):
        self.builder = _load("build_pangenome_panel", "05_build_pangenome_panel.py")
        self.targets = ["Azucena", "IR64", "N22"]

    def test_anchor_accepted_when_all_lift(self):
        lifted = {"Azucena": [100], "IR64": [200], "N22": [300]}
        self.assertTrue(self.builder.anchor_is_usable(lifted, self.targets))

    def test_anchor_rejected_when_one_genome_missing(self):
        lifted = {"Azucena": [100], "IR64": [200]}
        self.assertFalse(self.builder.anchor_is_usable(lifted, self.targets))

    def test_anchor_rejected_when_unlifted(self):
        lifted = {"Azucena": [100], "IR64": [200], "N22": []}
        self.assertFalse(self.builder.anchor_is_usable(lifted, self.targets))

    def test_multimapping_anchor_rejected(self):
        """paftools.js emits several intervals for an ambiguous anchor."""
        lifted = {"Azucena": [100, 5000], "IR64": [200], "N22": [300]}
        self.assertFalse(self.builder.anchor_is_usable(lifted, self.targets))


class TestLiftoverParsing(unittest.TestCase):
    """paftools.js discards the input BED name and rewrites column 4.

    Real output looks like:
        Azucena#1#Chr1  21495821  21495822  Nipponbare#1#Chr1_20862164_20862165  0  +
    so anchors are recovered by original reference coordinate, not by the
    name we supplied. These fixtures use the real format.
    """

    def setUp(self):
        self.builder = _load("build_pangenome_panel", "05_build_pangenome_panel.py")

    def test_source_key_matches_paftools_column_four(self):
        self.assertEqual(
            self.builder.source_key("Nipponbare#1#Chr1", 20862164),
            "Nipponbare#1#Chr1_20862164_20862165",
        )

    def test_parses_bed_keyed_by_source_coordinate(self):
        with tempfile.TemporaryDirectory() as tmp:
            bed = Path(tmp) / "lifted.bed"
            bed.write_text(
                "Azucena#1#Chr1\t100\t101\tNipponbare#1#Chr1_50_51\t0\t+\n"
                "Azucena#1#Chr1\t250\t251\tNipponbare#1#Chr1_70_71\t0\t+\n"
            )
            lifted = self.builder.parse_liftover_bed(bed)
            self.assertEqual(lifted["Nipponbare#1#Chr1_50_51"], [100])
            self.assertEqual(lifted["Nipponbare#1#Chr1_70_71"], [250])

    def test_collects_multiple_intervals_for_one_source(self):
        with tempfile.TemporaryDirectory() as tmp:
            bed = Path(tmp) / "lifted.bed"
            bed.write_text(
                "Azucena#1#Chr1\t100\t101\tNipponbare#1#Chr1_50_51\t0\t+\n"
                "Azucena#1#Chr1\t900\t901\tNipponbare#1#Chr1_50_51\t0\t+\n"
            )
            lifted = self.builder.parse_liftover_bed(bed)
            self.assertEqual(sorted(lifted["Nipponbare#1#Chr1_50_51"]), [100, 900])


class TestBackgroundTeFilter(unittest.TestCase):
    def setUp(self):
        self.builder = _load("build_pangenome_panel", "05_build_pangenome_panel.py")

    def test_site_near_same_family_te_is_rejected(self):
        background = {"Azucena": [("mPing", 10_000, 10_500)]}
        self.assertFalse(
            self.builder.site_clear_of_background(
                {"Azucena": 10_800}, "mPing", background, buffer_bp=1000
            )
        )

    def test_site_far_from_te_is_accepted(self):
        background = {"Azucena": [("mPing", 10_000, 10_500)]}
        self.assertTrue(
            self.builder.site_clear_of_background(
                {"Azucena": 50_000}, "mPing", background, buffer_bp=1000
            )
        )

    def test_different_family_does_not_block(self):
        background = {"Azucena": [("Tos17", 10_000, 10_500)]}
        self.assertTrue(
            self.builder.site_clear_of_background(
                {"Azucena": 10_800}, "mPing", background, buffer_bp=1000
            )
        )

    def test_blocked_by_any_genome(self):
        """A collision in any genome disqualifies the site everywhere."""
        background = {
            "Azucena": [],
            "IR64": [("mPing", 20_000, 20_400)],
        }
        self.assertFalse(
            self.builder.site_clear_of_background(
                {"Azucena": 50_000, "IR64": 20_500}, "mPing", background, 1000
            )
        )


class TestInsertion(unittest.TestCase):
    def setUp(self):
        self.builder = _load("build_pangenome_panel", "05_build_pangenome_panel.py")

    def test_descending_insertion_preserves_coordinates(self):
        """Later insertions must not shift the sites of earlier ones."""
        sequence = "A" * 1000
        events = [(100, "GGG", "TT"), (500, "CCC", "AA")]
        result = self.builder.insert_all(sequence, events)
        # Each insertion contributes its TE plus its TSD.
        self.assertEqual(len(result), 1000 + len("GGG") + 2 + len("CCC") + 2)
        self.assertEqual(result[100:105], "GGGTT")

    def test_insertion_is_idempotent_in_ordering(self):
        sequence = "A" * 1000
        forward = self.builder.insert_all(sequence, [(100, "GG", "T"), (500, "CC", "A")])
        reverse = self.builder.insert_all(sequence, [(500, "CC", "A"), (100, "GG", "T")])
        self.assertEqual(forward, reverse)


if __name__ == "__main__":
    unittest.main()


class TestSampleComposition(unittest.TestCase):
    """The sample must give ground truth for BOTH halves of the report."""

    def setUp(self):
        self.sample = _load("build_sample_genome", "06_build_sample_genome.py")

    @staticmethod
    def _panel_events(per_pattern=20):
        events = []
        patterns = ["{:04b}".format(i) for i in range(1, 16)]
        for pattern in patterns:
            for slot in range(per_pattern):
                events.append(
                    {
                        "event_id": f"PG{pattern}{slot:03d}",
                        "pattern": pattern,
                        "genomes": "Nipponbare",
                        "te_id": "Os0001#DNA/Foo",
                        "te_family": "Os0001",
                        "te_group": "MULE",
                        "strand": "+",
                        "tsd_length": "3",
                        "inserted_length": "500",
                        "anchor_id": f"anchor_{pattern}{slot:03d}",
                        "Nipponbare_pos": "1000",
                    }
                )
        return events

    def test_split_is_stratified_across_every_pattern(self):
        rng = random.Random(1)
        carried, absent = self.sample.choose_reference_events(
            self._panel_events(), 150, rng
        )
        self.assertEqual(len(carried), 150)
        self.assertEqual(len(absent), 300 - 150)
        carried_patterns = {e["pattern"] for e in carried}
        absent_patterns = {e["pattern"] for e in absent}
        # Every sharing pattern must contribute to BOTH sides, or attribution
        # goes untested for whichever patterns vanish.
        self.assertEqual(len(carried_patterns), 15)
        self.assertEqual(len(absent_patterns), 15)

    def test_split_is_disjoint_and_complete(self):
        rng = random.Random(2)
        events = self._panel_events()
        carried, absent = self.sample.choose_reference_events(events, 150, rng)
        ids_carried = {e["event_id"] for e in carried}
        ids_absent = {e["event_id"] for e in absent}
        self.assertEqual(ids_carried & ids_absent, set())
        self.assertEqual(len(ids_carried | ids_absent), len(events))

    def test_uneven_split_still_totals_correctly(self):
        rng = random.Random(3)
        carried, absent = self.sample.choose_reference_events(
            self._panel_events(), 157, rng
        )
        self.assertEqual(len(carried), 157)
        self.assertEqual(len(absent), 143)

    def test_used_anchor_ids(self):
        events = self._panel_events(per_pattern=2)
        self.assertEqual(len(self.sample.used_anchor_ids(events)), 30)
