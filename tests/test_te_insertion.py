"""Tests for the te_insertion module."""

import argparse
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from simulate_data.modules import te_insertion

DATA_DIR = Path(__file__).parent / "data"
MINI_GENOME = DATA_DIR / "mini_genome.fa"
MINI_TE = DATA_DIR / "mini_te.fa"
MINI_KNOWN_DEL = DATA_DIR / "mini_known_del.out"


def _mock_tevarsim_output(cmd):
    """Create the FASTA output expected from mocked TEvarSim Simulate."""
    if "Simulate" in cmd:
        outprefix = Path(cmd[cmd.index("--outprefix") + 1])
        chrom = outprefix.name.removeprefix("Sim_")
        outprefix.parent.mkdir(parents=True, exist_ok=True)
        outprefix.with_suffix(".fa").write_text(f">{chrom}_0\nNNNN\n")
    return MagicMock(returncode=0, stdout="", stderr="")


class TestRegisterParser:
    """Tests for register_parser."""

    def test_register_parser_creates_valid_parser(self):
        parser = argparse.ArgumentParser()
        te_insertion.register_parser(parser)

        args = parser.parse_args(
            [
                "--ref",
                "ref.fa",
                "--te",
                "te.fa",
                "--known-del",
                "del.out",
                "--num",
                "100",
                "--output",
                "results/",
            ]
        )
        assert args.ref == "ref.fa"
        assert args.te == "te.fa"
        assert args.known_del == "del.out"
        assert args.num == 100
        assert args.output == "results/"
        assert args.seed is None
        assert args.bed is None
        assert args.chroms == "all"
        assert args.num_genomes == 1
        assert args.ins_ratio == 0.6
        assert args.te_type is None
        assert args.snp_rate == 0.02
        assert args.indel_rate == 0.005
        assert args.indel_ins == 0.4
        assert args.indel_geom_p == 0.7
        assert args.truncated_ratio == 0.3
        assert args.truncated_max_length == 0.5
        assert args.polyA_ratio == 0.8
        assert args.polyA_min == 5
        assert args.polyA_max == 20
        assert args.sense_strand_ratio == 0.5

    def test_register_parser_with_all_options(self):
        parser = argparse.ArgumentParser()
        te_insertion.register_parser(parser)

        args = parser.parse_args(
            [
                "--ref",
                "ref.fa",
                "--te",
                "te.fa",
                "--known-del",
                "del.out",
                "--num",
                "50",
                "--output",
                "results/",
                "--seed",
                "42",
                "--bed",
                "positions.bed",
                "--chroms",
                "Chr1,Chr2",
                "--num-genomes",
                "3",
                "--ins-ratio",
                "0.8",
                "--te-type",
                "Harbinger",
                "--te-type",
                "Tourist,MuDR",
                "--snp-rate",
                "0.01",
                "--indel-rate",
                "0.002",
                "--indel-ins",
                "0.25",
                "--indel-geom-p",
                "0.8",
                "--truncated-ratio",
                "0.15",
                "--truncated-max-length",
                "0.4",
                "--polyA-ratio",
                "0.2",
                "--polyA-min",
                "3",
                "--polyA-max",
                "12",
                "--sense-strand-ratio",
                "0.7",
            ]
        )
        assert args.seed == 42
        assert args.bed == "positions.bed"
        assert args.chroms == "Chr1,Chr2"
        assert args.num_genomes == 3
        assert args.ins_ratio == 0.8
        assert args.te_type == ["Harbinger", "Tourist,MuDR"]
        assert args.snp_rate == 0.01
        assert args.indel_rate == 0.002
        assert args.indel_ins == 0.25
        assert args.indel_geom_p == 0.8
        assert args.truncated_ratio == 0.15
        assert args.truncated_max_length == 0.4
        assert args.polyA_ratio == 0.2
        assert args.polyA_min == 3
        assert args.polyA_max == 12
        assert args.sense_strand_ratio == 0.7

    def test_register_parser_missing_required(self):
        parser = argparse.ArgumentParser()
        te_insertion.register_parser(parser)

        with pytest.raises(SystemExit):
            parser.parse_args(["--ref", "ref.fa"])


class TestMain:
    """Tests for main function."""

    @patch("simulate_data.modules.te_insertion.extract_chromosomes")
    @patch("simulate_data.modules.te_insertion.check_tool_installed")
    @patch("simulate_data.modules.te_insertion.run_command")
    def test_main_all_chroms(self, mock_run, mock_check, mock_extract):
        """Test main with all chromosomes selected."""
        mock_run.side_effect = _mock_tevarsim_output

        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            te=str(MINI_TE),
            known_del=str(MINI_KNOWN_DEL),
            num=10,
            output="results/te_test/",
            seed=42,
            bed=None,
            chroms="all",
            num_genomes=1,
            ins_ratio=0.6,
        )

        te_insertion.main(ns)

        assert mock_run.call_count >= 2
        all_cmds = [call.args[0] for call in mock_run.call_args_list]
        terandom_cmds = [c for c in all_cmds if "TErandom" in c]
        simulate_cmds = [c for c in all_cmds if "Simulate" in c]
        assert len(terandom_cmds) >= 1
        assert len(simulate_cmds) >= 1

    @patch("simulate_data.modules.te_insertion.extract_chromosomes")
    @patch("simulate_data.modules.te_insertion.check_tool_installed")
    @patch("simulate_data.modules.te_insertion.run_command")
    def test_main_specific_chroms(self, mock_run, mock_check, mock_extract):
        """Test main with specific chromosomes selected."""
        mock_run.side_effect = _mock_tevarsim_output

        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            te=str(MINI_TE),
            known_del=str(MINI_KNOWN_DEL),
            num=5,
            output="results/te_test/",
            seed=None,
            bed=None,
            chroms="Chr1,Chr3",
            num_genomes=1,
            ins_ratio=0.6,
            te_type=["Harbinger"],
        )

        te_insertion.main(ns)

        assert mock_run.call_count == 4
        all_cmds = [call.args[0] for call in mock_run.call_args_list]
        terandom_cmds = [c for c in all_cmds if "TErandom" in c]
        assert all("--TEtype" in c and "Harbinger" in c for c in terandom_cmds)

    @patch("simulate_data.modules.te_insertion.extract_chromosomes")
    @patch("simulate_data.modules.te_insertion.check_tool_installed")
    @patch("simulate_data.modules.te_insertion.run_command")
    def test_main_range_chroms(self, mock_run, mock_check, mock_extract):
        """Test main with chromosome range Chr2-5."""
        mock_run.side_effect = _mock_tevarsim_output

        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            te=str(MINI_TE),
            known_del=str(MINI_KNOWN_DEL),
            num=20,
            output="results/te_test/",
            seed=None,
            bed=None,
            chroms="Chr2-5",
            num_genomes=1,
            ins_ratio=0.6,
        )

        te_insertion.main(ns)

        assert mock_run.call_count == 8

    @patch("simulate_data.modules.te_insertion.check_tool_installed")
    @patch("simulate_data.modules.te_insertion.run_command")
    def test_main_with_bed_file(self, mock_run, mock_check, tmp_path):
        """Test main with pre-generated BED file (skips TErandom)."""
        mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")

        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            te=str(MINI_TE),
            known_del=str(MINI_KNOWN_DEL),
            num=10,
            output=str(tmp_path / "te_output"),
            seed=None,
            bed=str(DATA_DIR / "mini_te.bed"),
            chroms="all",
            num_genomes=1,
            ins_ratio=0.6,
        )

        te_insertion.main(ns)

        mock_run.assert_called_once()
        cmd = mock_run.call_args[0][0]
        assert cmd[0] == "tevarsim"
        assert "Simulate" in cmd
        assert "--bed" in cmd

    def test_main_invalid_num(self):
        """Test that num <= 0 raises ValueError."""
        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            te=str(MINI_TE),
            known_del=str(MINI_KNOWN_DEL),
            num=0,
            output="results/",
            seed=None,
            bed=None,
            chroms="all",
            num_genomes=1,
            ins_ratio=0.6,
        )
        with pytest.raises(ValueError, match="must be positive"):
            te_insertion.main(ns)

    def test_main_invalid_divergence_rate(self):
        """Test that rate-style TEvarSim options are constrained to 0-1."""
        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            te=str(MINI_TE),
            known_del=str(MINI_KNOWN_DEL),
            num=10,
            output="results/",
            seed=None,
            bed=None,
            chroms="all",
            num_genomes=1,
            ins_ratio=0.6,
            snp_rate=1.5,
        )
        with pytest.raises(ValueError, match="snp-rate"):
            te_insertion.main(ns)

    def test_main_missing_ref(self):
        """Test that missing reference file raises FileNotFoundError."""
        ns = argparse.Namespace(
            ref="/nonexistent/ref.fa",
            te=str(MINI_TE),
            known_del=str(MINI_KNOWN_DEL),
            num=10,
            output="results/",
            seed=None,
            bed=None,
            chroms="all",
            num_genomes=1,
            ins_ratio=0.6,
        )
        with pytest.raises(FileNotFoundError):
            te_insertion.main(ns)

    def test_main_invalid_chromosomes(self):
        """Test that invalid chromosome names raise ValueError."""
        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            te=str(MINI_TE),
            known_del=str(MINI_KNOWN_DEL),
            num=10,
            output="results/",
            seed=None,
            bed=None,
            chroms="Chr99,Chr100",
            num_genomes=1,
            ins_ratio=0.6,
        )
        with pytest.raises(ValueError, match="not found"):
            te_insertion.main(ns)

    @patch("simulate_data.modules.te_insertion.extract_chromosomes")
    @patch("simulate_data.modules.te_insertion.check_tool_installed")
    @patch("simulate_data.modules.te_insertion.run_command")
    def test_main_no_seed(self, mock_run, mock_check, mock_extract):
        """Test that seed is not passed when None."""
        mock_run.side_effect = _mock_tevarsim_output

        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            te=str(MINI_TE),
            known_del=str(MINI_KNOWN_DEL),
            num=10,
            output="results/te_test/",
            seed=None,
            bed=None,
            chroms="all",
            num_genomes=1,
            ins_ratio=0.6,
        )

        te_insertion.main(ns)

        all_cmds = [call.args[0] for call in mock_run.call_args_list]
        for cmd in all_cmds:
            assert "--seed" not in cmd

    def test_main_num_genomes_requires_single_genome_for_chrom_subset(self):
        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            te=str(MINI_TE),
            known_del=str(MINI_KNOWN_DEL),
            num=10,
            output="results/",
            seed=None,
            bed=None,
            chroms="Chr1,Chr2",
            num_genomes=2,
            ins_ratio=0.6,
            te_type=None,
        )
        with patch("simulate_data.modules.te_insertion.check_tool_installed"):
            with pytest.raises(ValueError, match="num-genomes"):
                te_insertion.main(ns)


class TestBuildCommands:
    """Tests for command builder helpers."""

    def test_build_terandom_command_basic(self):
        cmd = te_insertion._build_terandom_command(
            te_fasta=Path("te.fa"),
            known_del=Path("del.out"),
            chrom="Chr1",
            num_te=100,
            outprefix="/tmp/terandom",
        )
        assert cmd[0] == "tevarsim"
        assert "TErandom" in cmd
        assert "--consensus" in cmd
        assert "--knownDEL" in cmd
        assert "--CHR" in cmd
        assert "--nTE" in cmd
        assert "--outprefix" in cmd
        assert "--ins-ratio" in cmd
        assert "--snp-rate" in cmd
        assert "--indel-rate" in cmd
        assert "--truncated-ratio" in cmd
        assert "--polyA-ratio" in cmd
        assert "--seed" not in cmd

    def test_build_terandom_command_with_seed(self):
        cmd = te_insertion._build_terandom_command(
            te_fasta=Path("te.fa"),
            known_del=Path("del.out"),
            chrom="Chr1",
            num_te=100,
            outprefix="/tmp/terandom",
            seed=42,
        )
        assert "--seed" in cmd

    def test_build_terandom_command_with_te_types(self):
        cmd = te_insertion._build_terandom_command(
            te_fasta=Path("te.fa"),
            known_del=Path("del.out"),
            chrom="Chr1",
            num_te=100,
            outprefix="/tmp/terandom",
            te_types=["Harbinger", "Tourist"],
        )
        assert cmd.count("--TEtype") == 2
        assert "Harbinger" in cmd
        assert "Tourist" in cmd

    def test_build_terandom_command_with_divergence_options(self):
        cmd = te_insertion._build_terandom_command(
            te_fasta=Path("te.fa"),
            known_del=Path("del.out"),
            chrom="Chr1",
            num_te=100,
            outprefix="/tmp/terandom",
            snp_rate=0.01,
            indel_rate=0.002,
            indel_ins=0.25,
            indel_geom_p=0.8,
            truncated_ratio=0.15,
            truncated_max_length=0.4,
            polya_ratio=0.2,
            polya_min=3,
            polya_max=12,
        )
        assert cmd[cmd.index("--snp-rate") + 1] == "0.01"
        assert cmd[cmd.index("--indel-rate") + 1] == "0.002"
        assert cmd[cmd.index("--indel-ins") + 1] == "0.25"
        assert cmd[cmd.index("--indel-geom-p") + 1] == "0.8"
        assert cmd[cmd.index("--truncated-ratio") + 1] == "0.15"
        assert cmd[cmd.index("--truncated-max-length") + 1] == "0.4"
        assert cmd[cmd.index("--polyA-ratio") + 1] == "0.2"
        assert cmd[cmd.index("--polyA-min") + 1] == "3"
        assert cmd[cmd.index("--polyA-max") + 1] == "12"

    def test_build_simulate_command_basic(self):
        cmd = te_insertion._build_simulate_command(
            ref_fasta=Path("ref.fa"),
            te_pool=Path("pool.fa"),
            bed_file=Path("pos.bed"),
            num_genomes=1,
            outprefix="/tmp/Sim",
        )
        assert cmd[0] == "tevarsim"
        assert "Simulate" in cmd
        assert "--ref" in cmd
        assert "--pool" in cmd
        assert "--bed" in cmd
        assert "--num" in cmd
        assert "--outprefix" in cmd
        assert "--sense-strand-ratio" in cmd
        assert "--seed" not in cmd

    def test_build_simulate_command_with_seed(self):
        cmd = te_insertion._build_simulate_command(
            ref_fasta=Path("ref.fa"),
            te_pool=Path("pool.fa"),
            bed_file=Path("pos.bed"),
            num_genomes=1,
            outprefix="/tmp/Sim",
            seed=42,
        )
        assert "--seed" in cmd

    def test_build_simulate_command_with_sense_strand_ratio(self):
        cmd = te_insertion._build_simulate_command(
            ref_fasta=Path("ref.fa"),
            te_pool=Path("pool.fa"),
            bed_file=Path("pos.bed"),
            num_genomes=1,
            outprefix="/tmp/Sim",
            sense_strand_ratio=0.75,
        )
        assert cmd[cmd.index("--sense-strand-ratio") + 1] == "0.75"


class TestRepeatMaskerGffConversion:
    """Tests for RepeatMasker GFF3 compatibility helpers."""

    def test_repeatmasker_gff_to_out(self, tmp_path):
        gff = tmp_path / "rm.gff"
        gff.write_text(
            "##gff-version 3\n"
            "Chr1\tRepeatMasker\tTransposon\t100\t530\t4057\t+\t.\t"
            "ID=TE1;Target=MPING 1 430;Class=DNA/Harbinger;"
            "PercDiv=0.0;PercDel=0.0;PercIns=0.0;\n"
        )
        out = tmp_path / "rm.out"

        te_insertion._repeatmasker_gff_to_out(gff, out)

        fields = out.read_text().split()
        assert fields[4] == "Chr1"
        assert fields[5] == "100"
        assert fields[6] == "530"
        assert fields[9] == "MPING"
        assert fields[10] == "DNA/Harbinger"

    def test_split_te_types(self):
        assert te_insertion._split_te_types(["Harbinger,Tourist", "MuDR"]) == [
            "Harbinger",
            "Tourist",
            "MuDR",
        ]
