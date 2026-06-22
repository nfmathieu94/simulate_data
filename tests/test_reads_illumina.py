"""Tests for the reads_illumina module."""

import argparse
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from simulate_data.modules import reads_illumina

DATA_DIR = Path(__file__).parent / "data"
MINI_GENOME = DATA_DIR / "mini_genome.fa"


class TestRegisterParser:
    """Tests for register_parser."""

    def test_register_parser_defaults(self):
        parser = argparse.ArgumentParser()
        reads_illumina.register_parser(parser)

        args = parser.parse_args(
            [
                "--ref",
                "ref.fa",
                "--output",
                "results/",
            ]
        )
        assert args.read_length == 150
        assert args.coverage == 20.0
        assert args.paired is True
        assert args.fragment_size == 300
        assert args.fragment_std == 30
        assert args.profile is None
        assert args.seed is None

    def test_register_parser_custom_values(self):
        parser = argparse.ArgumentParser()
        reads_illumina.register_parser(parser)

        args = parser.parse_args(
            [
                "--ref",
                "ref.fa",
                "--output",
                "results/",
                "--read-length",
                "100",
                "--coverage",
                "30",
                "--seed",
                "42",
                "--fragment-size",
                "500",
                "--fragment-std",
                "50",
            ]
        )
        assert args.read_length == 100
        assert args.coverage == 30.0
        assert args.seed == 42
        assert args.fragment_size == 500
        assert args.fragment_std == 50

    def test_register_parser_missing_required(self):
        parser = argparse.ArgumentParser()
        reads_illumina.register_parser(parser)

        with pytest.raises(SystemExit):
            parser.parse_args(["--ref", "ref.fa"])


class TestMain:
    """Tests for main function."""

    @patch("simulate_data.modules.reads_illumina.check_tool_installed")
    @patch("simulate_data.modules.reads_illumina.run_command")
    def test_main_success(self, mock_run, mock_check, tmp_path):
        mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")

        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            read_length=150,
            coverage=20.0,
            output=str(tmp_path / "illumina_output"),
            seed=42,
            paired=True,
            fragment_size=300,
            fragment_std=30,
            profile=None,
        )

        reads_illumina.main(ns)

        mock_run.assert_called_once()
        cmd = mock_run.call_args[0][0]
        assert cmd[0] == "art_illumina"
        assert "-i" in cmd
        assert "-l" in cmd
        assert "-f" in cmd
        assert "-rs" in cmd  # seed flag

    @patch("simulate_data.modules.reads_illumina.check_tool_installed")
    @patch("simulate_data.modules.reads_illumina.run_command")
    def test_main_no_seed(self, mock_run, mock_check, tmp_path):
        mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")

        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            read_length=150,
            coverage=20.0,
            output=str(tmp_path / "illumina_output"),
            seed=None,
            paired=True,
            fragment_size=300,
            fragment_std=30,
            profile=None,
        )

        reads_illumina.main(ns)

        cmd = mock_run.call_args[0][0]
        assert "-rs" not in cmd

    @patch("simulate_data.modules.reads_illumina.check_tool_installed")
    @patch("simulate_data.modules.reads_illumina.run_command")
    def test_main_with_profile(self, mock_run, mock_check, tmp_path):
        mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")

        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            read_length=150,
            coverage=20.0,
            output=str(tmp_path / "illumina_output"),
            seed=None,
            paired=True,
            fragment_size=300,
            fragment_std=30,
            profile="HiSeq",
        )

        reads_illumina.main(ns)

        cmd = mock_run.call_args[0][0]
        assert "-sp" in cmd
        assert "HiSeq" in cmd

    def test_main_invalid_read_length(self):
        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            read_length=0,
            coverage=20.0,
            output="results/",
            seed=None,
            paired=True,
            fragment_size=300,
            fragment_std=30,
            profile=None,
        )
        with pytest.raises(ValueError, match="read-length must be positive"):
            reads_illumina.main(ns)

    def test_main_invalid_coverage(self):
        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            read_length=150,
            coverage=0,
            output="results/",
            seed=None,
            paired=True,
            fragment_size=300,
            fragment_std=30,
            profile=None,
        )
        with pytest.raises(ValueError, match="coverage must be positive"):
            reads_illumina.main(ns)

    def test_main_invalid_fragment_size(self):
        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            read_length=150,
            coverage=20.0,
            output="results/",
            seed=None,
            paired=True,
            fragment_size=0,
            fragment_std=30,
            profile=None,
        )
        with pytest.raises(ValueError, match="fragment-size must be positive"):
            reads_illumina.main(ns)

    def test_main_missing_ref(self):
        ns = argparse.Namespace(
            ref="/nonexistent/ref.fa",
            read_length=150,
            coverage=20.0,
            output="results/",
            seed=None,
            paired=True,
            fragment_size=300,
            fragment_std=30,
            profile=None,
        )
        with pytest.raises(FileNotFoundError):
            reads_illumina.main(ns)
