"""Tests for the reads_pacbio module."""

import argparse
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from simulate_data.modules import reads_pacbio

DATA_DIR = Path(__file__).parent / "data"
MINI_GENOME = DATA_DIR / "mini_genome.fa"


class TestRegisterParser:
    """Tests for register_parser."""

    def test_register_parser_defaults(self):
        parser = argparse.ArgumentParser()
        reads_pacbio.register_parser(parser)

        args = parser.parse_args(
            [
                "--ref",
                "ref.fa",
                "--output",
                "results/",
            ]
        )
        assert args.coverage == 20.0
        assert args.read_type == "CLR"
        assert args.read_length == 15000
        assert args.read_std == 13000
        assert args.pass_num == 10
        assert args.error_model == "QSHMM-RSII"
        assert args.qscore_model is None
        assert args.seed is None

    def test_register_parser_hifi_mode(self):
        parser = argparse.ArgumentParser()
        reads_pacbio.register_parser(parser)

        args = parser.parse_args(
            [
                "--ref",
                "ref.fa",
                "--output",
                "results/",
                "--read-type",
                "HiFi",
                "--pass-num",
                "20",
            ]
        )
        assert args.read_type == "HiFi"
        assert args.pass_num == 20

    def test_register_parser_invalid_read_type(self):
        parser = argparse.ArgumentParser()
        reads_pacbio.register_parser(parser)

        with pytest.raises(SystemExit):
            parser.parse_args(
                [
                    "--ref",
                    "ref.fa",
                    "--output",
                    "results/",
                    "--read-type",
                    "InvalidType",
                ]
            )

    def test_register_parser_missing_required(self):
        parser = argparse.ArgumentParser()
        reads_pacbio.register_parser(parser)

        with pytest.raises(SystemExit):
            parser.parse_args(["--ref", "ref.fa"])


class TestMain:
    """Tests for main function."""

    @patch("simulate_data.modules.reads_pacbio.check_tool_installed")
    @patch("simulate_data.modules.reads_pacbio.run_command")
    def test_main_clr_mode(self, mock_run, mock_check, tmp_path):
        mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")

        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            coverage=20.0,
            output=str(tmp_path / "pacbio_output"),
            seed=42,
            read_type="CLR",
            read_length=15000,
            read_std=13000,
            pass_num=10,
            error_model="QSHMM-RSII",
            qscore_model=None,
        )

        reads_pacbio.main(ns)

        mock_run.assert_called_once()
        cmd = mock_run.call_args[0][0]
        assert cmd[0] == "pbsim"
        assert "--strategy" in cmd
        assert "wgs" in cmd
        assert "--seed" in cmd
        # CLR mode should not include --pass-num
        assert "--pass-num" not in cmd

    @patch("simulate_data.modules.reads_pacbio.check_tool_installed")
    @patch("simulate_data.modules.reads_pacbio.run_command")
    def test_main_hifi_mode(self, mock_run, mock_check, tmp_path):
        mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")

        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            coverage=20.0,
            output=str(tmp_path / "pacbio_output"),
            seed=None,
            read_type="HiFi",
            read_length=15000,
            read_std=13000,
            pass_num=15,
            error_model="QSHMM-RSII",
            qscore_model=None,
        )

        reads_pacbio.main(ns)

        mock_run.assert_called_once()
        cmd = mock_run.call_args[0][0]
        # HiFi mode should include --pass-num
        assert "--pass-num" in cmd
        assert "15" in cmd

    @patch("simulate_data.modules.reads_pacbio.check_tool_installed")
    @patch("simulate_data.modules.reads_pacbio.run_command")
    def test_main_no_seed(self, mock_run, mock_check, tmp_path):
        mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")

        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            coverage=20.0,
            output=str(tmp_path / "pacbio_output"),
            seed=None,
            read_type="CLR",
            read_length=15000,
            read_std=13000,
            pass_num=10,
            error_model="QSHMM-RSII",
            qscore_model=None,
        )

        reads_pacbio.main(ns)

        cmd = mock_run.call_args[0][0]
        assert "--seed" not in cmd

    def test_main_invalid_coverage(self):
        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            coverage=0,
            output="results/",
            seed=None,
            read_type="CLR",
            read_length=15000,
            read_std=13000,
            pass_num=10,
            error_model="QSHMM-RSII",
            qscore_model=None,
        )
        with pytest.raises(ValueError, match="coverage must be positive"):
            reads_pacbio.main(ns)

    def test_main_invalid_read_length(self):
        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            coverage=20.0,
            output="results/",
            seed=None,
            read_type="CLR",
            read_length=0,
            read_std=13000,
            pass_num=10,
            error_model="QSHMM-RSII",
            qscore_model=None,
        )
        with pytest.raises(ValueError, match="read-length must be positive"):
            reads_pacbio.main(ns)

    def test_main_negative_read_std(self):
        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            coverage=20.0,
            output="results/",
            seed=None,
            read_type="CLR",
            read_length=15000,
            read_std=-1,
            pass_num=10,
            error_model="QSHMM-RSII",
            qscore_model=None,
        )
        with pytest.raises(ValueError, match="read-std must be non-negative"):
            reads_pacbio.main(ns)

    def test_main_hifi_invalid_pass_num(self):
        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            coverage=20.0,
            output="results/",
            seed=None,
            read_type="HiFi",
            read_length=15000,
            read_std=13000,
            pass_num=0,
            error_model="QSHMM-RSII",
            qscore_model=None,
        )
        with pytest.raises(ValueError, match="pass-num must be positive"):
            reads_pacbio.main(ns)

    def test_main_missing_ref(self):
        ns = argparse.Namespace(
            ref="/nonexistent/ref.fa",
            coverage=20.0,
            output="results/",
            seed=None,
            read_type="CLR",
            read_length=15000,
            read_std=13000,
            pass_num=10,
            error_model="QSHMM-RSII",
            qscore_model=None,
        )
        with pytest.raises(FileNotFoundError):
            reads_pacbio.main(ns)
