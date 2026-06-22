"""Tests for the te_insertion module."""

import argparse
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from simulate_data.modules import te_insertion

DATA_DIR = Path(__file__).parent / "data"
MINI_GENOME = DATA_DIR / "mini_genome.fa"
MINI_TE = DATA_DIR / "mini_te.fa"


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
                "--num",
                "100",
                "--output",
                "results/",
            ]
        )
        assert args.ref == "ref.fa"
        assert args.te == "te.fa"
        assert args.num == 100
        assert args.output == "results/"
        assert args.seed is None
        assert args.bed is None
        assert args.chroms == "all"

    def test_register_parser_with_all_options(self):
        parser = argparse.ArgumentParser()
        te_insertion.register_parser(parser)

        args = parser.parse_args(
            [
                "--ref",
                "ref.fa",
                "--te",
                "te.fa",
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
            ]
        )
        assert args.seed == 42
        assert args.bed == "positions.bed"
        assert args.chroms == "Chr1,Chr2"

    def test_register_parser_missing_required(self):
        parser = argparse.ArgumentParser()
        te_insertion.register_parser(parser)

        with pytest.raises(SystemExit):
            parser.parse_args(["--ref", "ref.fa"])


class TestMain:
    """Tests for main function."""

    @patch("simulate_data.modules.te_insertion.check_tool_installed")
    @patch("simulate_data.modules.te_insertion.run_command")
    def test_main_all_chroms(self, mock_run, mock_check):
        """Test main with all chromosomes selected."""
        mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")

        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            te=str(MINI_TE),
            num=10,
            output="results/te_test/",
            seed=42,
            bed=None,
            chroms="all",
        )

        te_insertion.main(ns)

        mock_run.assert_called_once()
        cmd = mock_run.call_args[0][0]
        assert cmd[0] == "tevarsim"
        assert "Simulate" in cmd
        assert "--ref" in cmd
        assert "--seed" in cmd

    @patch("simulate_data.modules.te_insertion.check_tool_installed")
    @patch("simulate_data.modules.te_insertion.run_command")
    def test_main_specific_chroms(self, mock_run, mock_check, tmp_path):
        """Test main with specific chromosomes selected."""
        mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")

        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            te=str(MINI_TE),
            num=5,
            output=str(tmp_path / "te_output"),
            seed=None,
            bed=None,
            chroms="Chr1,Chr3",
        )

        te_insertion.main(ns)

        mock_run.assert_called_once()
        cmd = mock_run.call_args[0][0]
        assert cmd[0] == "tevarsim"

    @patch("simulate_data.modules.te_insertion.check_tool_installed")
    @patch("simulate_data.modules.te_insertion.run_command")
    def test_main_range_chroms(self, mock_run, mock_check, tmp_path):
        """Test main with chromosome range Chr2-5."""
        mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")

        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            te=str(MINI_TE),
            num=20,
            output=str(tmp_path / "te_output"),
            seed=None,
            bed=None,
            chroms="Chr2-5",
        )

        te_insertion.main(ns)

        mock_run.assert_called_once()

    @patch("simulate_data.modules.te_insertion.check_tool_installed")
    @patch("simulate_data.modules.te_insertion.run_command")
    def test_main_with_bed_file(self, mock_run, mock_check, tmp_path):
        """Test main with pre-generated BED file."""
        mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")

        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            te=str(MINI_TE),
            num=10,
            output=str(tmp_path / "te_output"),
            seed=None,
            bed=str(DATA_DIR / "mini_te.bed"),
            chroms="all",
        )

        te_insertion.main(ns)

        mock_run.assert_called_once()
        cmd = mock_run.call_args[0][0]
        assert "--bed" in cmd

    def test_main_invalid_num(self):
        """Test that num <= 0 raises ValueError."""
        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            te=str(MINI_TE),
            num=0,
            output="results/",
            seed=None,
            bed=None,
            chroms="all",
        )
        with pytest.raises(ValueError, match="must be positive"):
            te_insertion.main(ns)

    def test_main_missing_ref(self):
        """Test that missing reference file raises FileNotFoundError."""
        ns = argparse.Namespace(
            ref="/nonexistent/ref.fa",
            te=str(MINI_TE),
            num=10,
            output="results/",
            seed=None,
            bed=None,
            chroms="all",
        )
        with pytest.raises(FileNotFoundError):
            te_insertion.main(ns)

    def test_main_invalid_chromosomes(self):
        """Test that invalid chromosome names raise ValueError."""
        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            te=str(MINI_TE),
            num=10,
            output="results/",
            seed=None,
            bed=None,
            chroms="Chr99,Chr100",
        )
        with pytest.raises(ValueError, match="not found"):
            te_insertion.main(ns)

    @patch("simulate_data.modules.te_insertion.check_tool_installed")
    @patch("simulate_data.modules.te_insertion.run_command")
    def test_main_no_seed(self, mock_run, mock_check, tmp_path):
        """Test that seed is not passed when None."""
        mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")

        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            te=str(MINI_TE),
            num=10,
            output=str(tmp_path / "te_output"),
            seed=None,
            bed=None,
            chroms="all",
        )

        te_insertion.main(ns)

        cmd = mock_run.call_args[0][0]
        assert "--seed" not in cmd
