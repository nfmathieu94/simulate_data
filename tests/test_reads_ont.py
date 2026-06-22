"""Tests for the reads_ont module."""

import argparse
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from simulate_data.modules import reads_ont

DATA_DIR = Path(__file__).parent / "data"
MINI_GENOME = DATA_DIR / "mini_genome.fa"


class TestRegisterParser:
    """Tests for register_parser."""

    def test_register_parser_defaults(self):
        parser = argparse.ArgumentParser()
        reads_ont.register_parser(parser)

        args = parser.parse_args(
            [
                "--ref",
                "ref.fa",
                "--output",
                "results/",
            ]
        )
        assert args.coverage == 20.0
        assert args.read_length == 9000
        assert args.read_std == 7000
        assert args.error_model == "QSHMM-ONT"
        assert args.qscore_model is None
        assert args.seed is None

    def test_register_parser_custom_values(self):
        parser = argparse.ArgumentParser()
        reads_ont.register_parser(parser)

        args = parser.parse_args(
            [
                "--ref",
                "ref.fa",
                "--output",
                "results/",
                "--coverage",
                "30",
                "--read-length",
                "5000",
                "--read-std",
                "3000",
                "--seed",
                "42",
                "--error-model",
                "custom_model",
                "--qscore-model",
                "custom_qscore",
            ]
        )
        assert args.coverage == 30.0
        assert args.read_length == 5000
        assert args.read_std == 3000
        assert args.seed == 42
        assert args.error_model == "custom_model"
        assert args.qscore_model == "custom_qscore"

    def test_register_parser_missing_required(self):
        parser = argparse.ArgumentParser()
        reads_ont.register_parser(parser)

        with pytest.raises(SystemExit):
            parser.parse_args(["--ref", "ref.fa"])


class TestMain:
    """Tests for main function."""

    @patch("simulate_data.modules.reads_ont.check_tool_installed")
    @patch("simulate_data.modules.reads_ont.run_command")
    def test_main_success(self, mock_run, mock_check, tmp_path):
        mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")

        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            coverage=20.0,
            output=str(tmp_path / "ont_output"),
            seed=42,
            read_length=9000,
            read_std=7000,
            error_model="QSHMM-ONT",
            qscore_model=None,
        )

        reads_ont.main(ns)

        mock_run.assert_called_once()
        cmd = mock_run.call_args[0][0]
        assert cmd[0] == "pbsim"
        assert "--strategy" in cmd
        assert "wgs" in cmd
        assert "--depth" in cmd
        assert "--seed" in cmd

    @patch("simulate_data.modules.reads_ont.check_tool_installed")
    @patch("simulate_data.modules.reads_ont.run_command")
    def test_main_no_seed(self, mock_run, mock_check, tmp_path):
        mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")

        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            coverage=20.0,
            output=str(tmp_path / "ont_output"),
            seed=None,
            read_length=9000,
            read_std=7000,
            error_model="QSHMM-ONT",
            qscore_model=None,
        )

        reads_ont.main(ns)

        cmd = mock_run.call_args[0][0]
        assert "--seed" not in cmd

    @patch("simulate_data.modules.reads_ont.check_tool_installed")
    @patch("simulate_data.modules.reads_ont.run_command")
    def test_main_with_qscore_model(self, mock_run, mock_check, tmp_path):
        mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")

        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            coverage=20.0,
            output=str(tmp_path / "ont_output"),
            seed=None,
            read_length=9000,
            read_std=7000,
            error_model="QSHMM-ONT",
            qscore_model="custom_qscore_model",
        )

        reads_ont.main(ns)

        cmd = mock_run.call_args[0][0]
        assert "--qscore-model" in cmd
        assert "custom_qscore_model" in cmd

    def test_main_invalid_coverage(self):
        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            coverage=0,
            output="results/",
            seed=None,
            read_length=9000,
            read_std=7000,
            error_model="QSHMM-ONT",
            qscore_model=None,
        )
        with pytest.raises(ValueError, match="coverage must be positive"):
            reads_ont.main(ns)

    def test_main_invalid_read_length(self):
        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            coverage=20.0,
            output="results/",
            seed=None,
            read_length=0,
            read_std=7000,
            error_model="QSHMM-ONT",
            qscore_model=None,
        )
        with pytest.raises(ValueError, match="read-length must be positive"):
            reads_ont.main(ns)

    def test_main_negative_read_std(self):
        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            coverage=20.0,
            output="results/",
            seed=None,
            read_length=9000,
            read_std=-1,
            error_model="QSHMM-ONT",
            qscore_model=None,
        )
        with pytest.raises(ValueError, match="read-std must be non-negative"):
            reads_ont.main(ns)

    def test_main_missing_ref(self):
        ns = argparse.Namespace(
            ref="/nonexistent/ref.fa",
            coverage=20.0,
            output="results/",
            seed=None,
            read_length=9000,
            read_std=7000,
            error_model="QSHMM-ONT",
            qscore_model=None,
        )
        with pytest.raises(FileNotFoundError):
            reads_ont.main(ns)
