"""Tests for the sv_placement module."""

import argparse
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from simulate_data.modules import sv_placement

DATA_DIR = Path(__file__).parent / "data"
MINI_GENOME = DATA_DIR / "mini_genome.fa"


class TestRegisterParser:
    """Tests for register_parser."""

    def test_register_parser_creates_valid_parser(self):
        parser = argparse.ArgumentParser()
        sv_placement.register_parser(parser)

        args = parser.parse_args(
            [
                "--ref",
                "ref.fa",
                "--num-sv",
                "100",
                "--output",
                "results/",
            ]
        )
        assert args.ref == "ref.fa"
        assert args.num_sv == 100
        assert args.output == "results/"
        assert args.sv_types == "DEL,DUP,INV,TRA"
        assert args.seed is None
        assert args.config is None
        assert args.chroms == "all"

    def test_register_parser_with_options(self):
        parser = argparse.ArgumentParser()
        sv_placement.register_parser(parser)

        args = parser.parse_args(
            [
                "--ref",
                "ref.fa",
                "--num-sv",
                "50",
                "--output",
                "results/",
                "--sv-types",
                "DEL,INV",
                "--seed",
                "42",
                "--chroms",
                "Chr1-3",
            ]
        )
        assert args.sv_types == "DEL,INV"
        assert args.seed == 42
        assert args.chroms == "Chr1-3"

    def test_register_parser_missing_required(self):
        parser = argparse.ArgumentParser()
        sv_placement.register_parser(parser)

        with pytest.raises(SystemExit):
            parser.parse_args(["--ref", "ref.fa"])


class TestMain:
    """Tests for main function."""

    @patch("simulate_data.modules.sv_placement.check_tool_installed")
    @patch("simulate_data.modules.sv_placement.run_command")
    def test_main_all_chroms(self, mock_run, mock_check, tmp_path):
        """Test main with all chromosomes selected."""
        mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")

        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            num_sv=10,
            output=str(tmp_path / "sv_output"),
            sv_types="DEL,DUP,INV,TRA",
            seed=42,
            config=None,
            chroms="all",
        )

        sv_placement.main(ns)

        mock_run.assert_called_once()
        cmd = mock_run.call_args[0][0]
        assert cmd[0] == "SURVIVOR"
        assert "simSV" in cmd

    @patch("simulate_data.modules.sv_placement.check_tool_installed")
    @patch("simulate_data.modules.sv_placement.run_command")
    def test_main_specific_chroms(self, mock_run, mock_check, tmp_path):
        """Test main with specific chromosomes selected."""
        mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")

        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            num_sv=5,
            output=str(tmp_path / "sv_output"),
            sv_types="DEL,DUP,INV,TRA",
            seed=None,
            config=None,
            chroms="Chr1,Chr3",
        )

        sv_placement.main(ns)

        mock_run.assert_called_once()
        cmd = mock_run.call_args[0][0]
        assert cmd[0] == "SURVIVOR"

    @patch("simulate_data.modules.sv_placement.check_tool_installed")
    @patch("simulate_data.modules.sv_placement.run_command")
    def test_main_range_chroms(self, mock_run, mock_check, tmp_path):
        """Test main with chromosome range Chr2-5."""
        mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")

        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            num_sv=20,
            output=str(tmp_path / "sv_output"),
            sv_types="DEL,DUP,INV,TRA",
            seed=None,
            config=None,
            chroms="Chr2-5",
        )

        sv_placement.main(ns)

        mock_run.assert_called_once()

    @patch("simulate_data.modules.sv_placement.check_tool_installed")
    @patch("simulate_data.modules.sv_placement.run_command")
    def test_main_custom_sv_types(self, mock_run, mock_check, tmp_path):
        """Test main with custom SV types."""
        mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")

        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            num_sv=10,
            output=str(tmp_path / "sv_output"),
            sv_types="DEL,INV",
            seed=None,
            config=None,
            chroms="all",
        )

        sv_placement.main(ns)

        mock_run.assert_called_once()

    @patch("simulate_data.modules.sv_placement.check_tool_installed")
    @patch("simulate_data.modules.sv_placement.run_command")
    def test_main_with_config_file(self, mock_run, mock_check, tmp_path):
        """Test main with pre-generated config file."""
        config_path = tmp_path / "survivor_config.txt"
        config_path.write_text("# SURVIVOR config\nSV_TYPES=DEL\nNUM_SV=10\n")

        mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")

        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            num_sv=10,
            output=str(tmp_path / "sv_output"),
            sv_types="DEL,DUP,INV,TRA",
            seed=None,
            config=str(config_path),
            chroms="all",
        )

        sv_placement.main(ns)

        mock_run.assert_called_once()

    def test_main_invalid_num_sv(self):
        """Test that num_sv <= 0 raises ValueError."""
        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            num_sv=0,
            output="results/",
            sv_types="DEL,DUP,INV,TRA",
            seed=None,
            config=None,
            chroms="all",
        )
        with pytest.raises(ValueError, match="must be positive"):
            sv_placement.main(ns)

    def test_main_missing_ref(self):
        """Test that missing reference file raises FileNotFoundError."""
        ns = argparse.Namespace(
            ref="/nonexistent/ref.fa",
            num_sv=10,
            output="results/",
            sv_types="DEL,DUP,INV,TRA",
            seed=None,
            config=None,
            chroms="all",
        )
        with pytest.raises(FileNotFoundError):
            sv_placement.main(ns)

    def test_main_invalid_sv_type(self):
        """Test that invalid SV types raise ValueError."""
        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            num_sv=10,
            output="results/",
            sv_types="DEL,INVALID",
            seed=None,
            config=None,
            chroms="all",
        )
        with pytest.raises(ValueError, match="Invalid SV types"):
            sv_placement.main(ns)

    def test_main_invalid_chromosomes(self):
        """Test that invalid chromosome names raise ValueError."""
        ns = argparse.Namespace(
            ref=str(MINI_GENOME),
            num_sv=10,
            output="results/",
            sv_types="DEL,DUP,INV,TRA",
            seed=None,
            config=None,
            chroms="Chr99,Chr100",
        )
        with pytest.raises(ValueError, match="not found"):
            sv_placement.main(ns)


class TestGenerateSurvivorConfig:
    """Tests for _generate_survivor_config."""

    def test_generate_config_creates_file(self, tmp_path):
        config_path = tmp_path / "config.txt"
        result = sv_placement._generate_survivor_config(
            config_path=config_path,
            num_sv=50,
            sv_types=["DEL", "INV"],
            seed=123,
        )
        assert result == config_path
        assert config_path.exists()
        content = config_path.read_text()
        assert "NUM_SV=50" in content
        assert "DEL,INV" in content
        assert "SEED=123" in content

    def test_generate_config_no_seed(self, tmp_path):
        config_path = tmp_path / "config.txt"
        sv_placement._generate_survivor_config(
            config_path=config_path,
            num_sv=10,
            sv_types=["DEL"],
            seed=None,
        )
        content = config_path.read_text()
        assert "SEED" not in content
