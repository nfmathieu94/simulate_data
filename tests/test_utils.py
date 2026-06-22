"""Tests for simulate_data.utils."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from simulate_data.utils import (
    check_tool_installed,
    ensure_output_dir,
    extract_chromosomes,
    get_chromosomes_from_fasta,
    merge_fasta,
    parse_chromosome_spec,
    run_command,
    validate_fasta,
    validate_file_exists,
)

DATA_DIR = Path(__file__).parent / "data"
MINI_GENOME = DATA_DIR / "mini_genome.fa"
MINI_TE = DATA_DIR / "mini_te.fa"


class TestRunCommand:
    """Tests for run_command."""

    @patch("simulate_data.utils.subprocess.run")
    def test_run_command_success(self, mock_run):
        mock_run.return_value = MagicMock(returncode=0, stdout="ok", stderr="")
        result = run_command(["echo", "hello"])
        assert result.returncode == 0
        mock_run.assert_called_once()

    @patch("simulate_data.utils.subprocess.run")
    def test_run_command_failure_with_check(self, mock_run):
        mock_run.return_value = MagicMock(returncode=1, stdout="", stderr="error")
        with pytest.raises(RuntimeError, match="Command failed"):
            run_command(["false"])

    @patch("simulate_data.utils.subprocess.run")
    def test_run_command_failure_without_check(self, mock_run):
        mock_run.return_value = MagicMock(returncode=1, stdout="", stderr="error")
        result = run_command(["false"], check=False)
        assert result.returncode == 1


class TestValidateFileExists:
    """Tests for validate_file_exists."""

    def test_validate_file_exists_true(self):
        validate_file_exists(MINI_GENOME)

    def test_validate_file_exists_missing(self):
        with pytest.raises(FileNotFoundError, match="File not found"):
            validate_file_exists(Path("/nonexistent/file.fa"))

    def test_validate_file_exists_directory(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="not a file"):
            validate_file_exists(tmp_path)


class TestValidateFasta:
    """Tests for validate_fasta."""

    def test_validate_fasta_valid(self):
        validate_fasta(MINI_GENOME)

    def test_validate_fasta_missing(self):
        with pytest.raises(FileNotFoundError):
            validate_fasta(Path("/nonexistent/file.fa"))

    def test_validate_fasta_no_header(self, tmp_path):
        fasta = tmp_path / "no_header.fa"
        fasta.write_text("ACGTACGT\n")
        with pytest.raises(ValueError, match="must start with '>'"):
            validate_fasta(fasta)

    def test_validate_fasta_empty_sequence(self, tmp_path):
        fasta = tmp_path / "empty_seq.fa"
        fasta.write_text(">Chr1\n\n>Chr2\n\n")
        with pytest.raises(ValueError, match="no sequence data"):
            validate_fasta(fasta)


class TestCheckToolInstalled:
    """Tests for check_tool_installed."""

    @patch("simulate_data.utils.shutil.which")
    def test_check_tool_installed_found(self, mock_which):
        mock_which.return_value = "/usr/bin/python3"
        result = check_tool_installed("python3")
        assert result == "/usr/bin/python3"

    @patch("simulate_data.utils.shutil.which")
    def test_check_tool_installed_not_found(self, mock_which):
        mock_which.return_value = None
        with pytest.raises(FileNotFoundError, match="not found on PATH"):
            check_tool_installed("nonexistent_tool_xyz")


class TestEnsureOutputDir:
    """Tests for ensure_output_dir."""

    def test_ensure_output_dir_creates_new(self, tmp_path):
        new_dir = tmp_path / "output" / "subdir"
        result = ensure_output_dir(new_dir)
        assert result == new_dir
        assert new_dir.exists()
        assert new_dir.is_dir()

    def test_ensure_output_dir_existing(self, tmp_path):
        existing_dir = tmp_path / "existing"
        existing_dir.mkdir()
        result = ensure_output_dir(existing_dir)
        assert result == existing_dir
        assert existing_dir.exists()


class TestGetChromosomesFromFasta:
    """Tests for get_chromosomes_from_fasta."""

    def test_get_chromosomes_from_mini_genome(self):
        chroms = get_chromosomes_from_fasta(MINI_GENOME)
        assert chroms == [
            "Chr1",
            "Chr2",
            "Chr3",
            "Chr4",
            "Chr5",
            "Chr6",
            "Chr7",
            "Chr8",
            "Chr9",
            "Chr10",
            "Chr11",
            "Chr12",
        ]

    def test_get_chromosomes_missing_file(self):
        with pytest.raises(FileNotFoundError):
            get_chromosomes_from_fasta(Path("/nonexistent/file.fa"))

    def test_get_chromosomes_no_headers(self, tmp_path):
        fasta = tmp_path / "no_headers.fa"
        fasta.write_text("ACGTACGT\n")
        with pytest.raises(ValueError, match="No chromosome headers"):
            get_chromosomes_from_fasta(fasta)


class TestParseChromosomeSpec:
    """Tests for parse_chromosome_spec."""

    AVAILABLE = [
        "Chr1",
        "Chr2",
        "Chr3",
        "Chr4",
        "Chr5",
        "Chr6",
        "Chr7",
        "Chr8",
        "Chr9",
        "Chr10",
        "Chr11",
        "Chr12",
    ]

    def test_parse_all(self):
        result = parse_chromosome_spec("all", self.AVAILABLE)
        assert result == self.AVAILABLE

    def test_parse_all_uppercase(self):
        result = parse_chromosome_spec("ALL", self.AVAILABLE)
        assert result == self.AVAILABLE

    def test_parse_single_chromosome(self):
        result = parse_chromosome_spec("Chr1", self.AVAILABLE)
        assert result == ["Chr1"]

    def test_parse_comma_separated(self):
        result = parse_chromosome_spec("Chr1,Chr3,Chr5", self.AVAILABLE)
        assert result == ["Chr1", "Chr3", "Chr5"]

    def test_parse_range_short(self):
        """Chr2-5 should expand to Chr2, Chr3, Chr4, Chr5."""
        result = parse_chromosome_spec("Chr2-5", self.AVAILABLE)
        assert result == ["Chr2", "Chr3", "Chr4", "Chr5"]

    def test_parse_range_full(self):
        """Chr2-Chr5 should expand to Chr2, Chr3, Chr4, Chr5."""
        result = parse_chromosome_spec("Chr2-Chr5", self.AVAILABLE)
        assert result == ["Chr2", "Chr3", "Chr4", "Chr5"]

    def test_parse_range_with_single(self):
        """Chr1-3,Chr8 should expand correctly."""
        result = parse_chromosome_spec("Chr1-3,Chr8", self.AVAILABLE)
        assert result == ["Chr1", "Chr2", "Chr3", "Chr8"]

    def test_parse_range_double_digit(self):
        """Chr10-12 should expand to Chr10, Chr11, Chr12."""
        result = parse_chromosome_spec("Chr10-12", self.AVAILABLE)
        assert result == ["Chr10", "Chr11", "Chr12"]

    def test_parse_invalid_chromosome(self):
        with pytest.raises(ValueError, match="not found in reference"):
            parse_chromosome_spec("Chr99", self.AVAILABLE)

    def test_parse_empty_spec(self):
        with pytest.raises(ValueError, match="No chromosomes specified"):
            parse_chromosome_spec("", self.AVAILABLE)

    def test_parse_whitespace_handling(self):
        result = parse_chromosome_spec(" Chr1 , Chr2 ", self.AVAILABLE)
        assert result == ["Chr1", "Chr2"]

    def test_parse_duplicates_removed(self):
        result = parse_chromosome_spec("Chr1,Chr1,Chr2", self.AVAILABLE)
        assert result == ["Chr1", "Chr2"]

    def test_parse_reverse_range_error(self):
        with pytest.raises(ValueError, match="start.* > end"):
            parse_chromosome_spec("Chr5-Chr2", self.AVAILABLE)


class TestExtractChromosomes:
    """Tests for extract_chromosomes."""

    def test_extract_single_chromosome(self, tmp_path):
        output = tmp_path / "subset.fa"
        extract_chromosomes(MINI_GENOME, ["Chr1"], output)
        assert output.exists()
        chroms = get_chromosomes_from_fasta(output)
        assert chroms == ["Chr1"]

    def test_extract_multiple_chromosomes(self, tmp_path):
        output = tmp_path / "subset.fa"
        extract_chromosomes(MINI_GENOME, ["Chr1", "Chr3", "Chr5"], output)
        assert output.exists()
        chroms = get_chromosomes_from_fasta(output)
        assert chroms == ["Chr1", "Chr3", "Chr5"]

    def test_extract_all_chromosomes(self, tmp_path):
        output = tmp_path / "subset.fa"
        all_chroms = get_chromosomes_from_fasta(MINI_GENOME)
        extract_chromosomes(MINI_GENOME, all_chroms, output)
        assert output.exists()
        chroms = get_chromosomes_from_fasta(output)
        assert chroms == all_chroms

    def test_extract_missing_chromosome(self, tmp_path):
        output = tmp_path / "subset.fa"
        with pytest.raises(ValueError, match="not found in FASTA"):
            extract_chromosomes(MINI_GENOME, ["Chr99"], output)

    def test_extract_preserves_sequence_content(self, tmp_path):
        output = tmp_path / "subset.fa"
        extract_chromosomes(MINI_GENOME, ["Chr2"], output)

        original_content = ""
        with open(MINI_GENOME, "r") as f:
            in_chr2 = False
            for line in f:
                if line.startswith(">Chr2"):
                    in_chr2 = True
                    original_content += line
                elif line.startswith(">") and in_chr2:
                    break
                elif in_chr2:
                    original_content += line

        extracted_content = output.read_text()
        assert extracted_content.strip() == original_content.strip()


class TestMergeFasta:
    """Tests for merge_fasta."""

    def test_merge_replaces_modified_chromosomes(self, tmp_path):
        """Modified chromosomes should replace base chromosomes in output."""
        base = tmp_path / "base.fa"
        base.write_text(">Chr1\nACGTACGT\n>Chr2\nTTTTAAAA\n>Chr3\nGGGGCCCC\n")

        modified = tmp_path / "modified.fa"
        modified.write_text(">Chr1\nNNNNNNNN\n>Chr2\nNNNNNNNN\n")

        output = tmp_path / "merged.fa"
        merge_fasta(
            base_fasta=base,
            modified_fasta=modified,
            output_path=output,
            modified_chroms={"Chr1", "Chr2"},
        )

        assert output.exists()
        content = output.read_text()

        # Chr1 and Chr2 should be from modified (all N's)
        assert ">Chr1\nNNNNNNNN" in content
        assert ">Chr2\nNNNNNNNN" in content

        # Chr3 should be from base (unchanged)
        assert ">Chr3\nGGGGCCCC" in content

    def test_merge_keeps_unmodified_chromosomes(self, tmp_path):
        base = tmp_path / "base.fa"
        base.write_text(">Chr1\nACGT\n>Chr2\nTTTT\n")

        modified = tmp_path / "modified.fa"
        modified.write_text(">Chr1\nNNNN\n")

        output = tmp_path / "merged.fa"
        merge_fasta(
            base_fasta=base,
            modified_fasta=modified,
            output_path=output,
            modified_chroms={"Chr1"},
        )

        content = output.read_text()
        assert ">Chr1\nNNNN" in content
        assert ">Chr2\nTTTT" in content

    def test_merge_with_real_mini_genome(self, tmp_path):
        """Test merge using the actual mini_genome.fa test fixture."""
        # Create a modified FASTA that replaces Chr1 with N's
        modified = tmp_path / "modified.fa"
        modified.write_text(">Chr1\nNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n")

        output = tmp_path / "merged.fa"
        merge_fasta(
            base_fasta=MINI_GENOME,
            modified_fasta=modified,
            output_path=output,
            modified_chroms={"Chr1"},
        )

        assert output.exists()
        chroms = get_chromosomes_from_fasta(output)
        # All 12 chromosomes should be present
        assert len(chroms) == 12
        assert chroms[0] == "Chr1"
        assert chroms[-1] == "Chr12"
