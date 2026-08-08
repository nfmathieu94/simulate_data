"""Tests for the shared long-read simulation helpers."""

import gzip

import pytest

from simulate_data.longread import (
    GZIP_LEVEL,
    PLATFORMS,
    build_pbsim_command,
    contig_name_for,
    iter_maf_alignments,
    pool_fastq_gz,
    resolve_model_path,
)


@pytest.fixture
def model_dir(tmp_path, monkeypatch):
    """A fake CONDA_PREFIX/data holding every shipped PBSIM3 model."""
    data = tmp_path / "data"
    data.mkdir()
    for name in (
        "QSHMM-ONT",
        "QSHMM-ONT-HQ",
        "QSHMM-RSII",
        "ERRHMM-ONT",
        "ERRHMM-ONT-HQ",
        "ERRHMM-RSII",
        "ERRHMM-SEQUEL",
    ):
        (data / f"{name}.model").write_text("x")
    monkeypatch.setenv("CONDA_PREFIX", str(tmp_path))
    return data


class TestResolveModelPath:
    """PBSIM3 requires a path to a .model file, not a bare model name."""

    def test_resolves_bare_name(self, tmp_path, monkeypatch):
        data = tmp_path / "data"
        data.mkdir()
        (data / "QSHMM-ONT-HQ.model").write_text("x")
        monkeypatch.setenv("CONDA_PREFIX", str(tmp_path))
        assert resolve_model_path("QSHMM-ONT-HQ") == data / "QSHMM-ONT-HQ.model"

    def test_resolves_name_with_suffix(self, tmp_path, monkeypatch):
        data = tmp_path / "data"
        data.mkdir()
        (data / "ERRHMM-SEQUEL.model").write_text("x")
        monkeypatch.setenv("CONDA_PREFIX", str(tmp_path))
        assert resolve_model_path("ERRHMM-SEQUEL.model") == data / "ERRHMM-SEQUEL.model"

    def test_explicit_existing_path_passes_through(self, tmp_path):
        model = tmp_path / "custom.model"
        model.write_text("x")
        assert resolve_model_path(str(model)) == model

    def test_missing_model_lists_available(self, tmp_path, monkeypatch):
        data = tmp_path / "data"
        data.mkdir()
        (data / "QSHMM-ONT.model").write_text("x")
        monkeypatch.setenv("CONDA_PREFIX", str(tmp_path))
        with pytest.raises(FileNotFoundError, match="QSHMM-ONT"):
            resolve_model_path("NOPE")

    def test_unset_conda_prefix_is_reported(self, tmp_path, monkeypatch):
        monkeypatch.delenv("CONDA_PREFIX", raising=False)
        with pytest.raises(RuntimeError, match="CONDA_PREFIX"):
            resolve_model_path("QSHMM-ONT")


class TestBuildPbsimCommand:
    """Platform presets must select the right method, model, and error profile."""

    def test_ont_hq_uses_qshmm_and_ont_difference_ratio(self, tmp_path, model_dir):
        cmd = build_pbsim_command(
            PLATFORMS["ont-hq"], tmp_path / "g.fa", 10.0, tmp_path / "out", seed=7
        )
        assert cmd[0] == "pbsim"
        assert cmd[cmd.index("--method") + 1] == "qshmm"
        assert cmd[cmd.index("--qshmm") + 1] == str(model_dir / "QSHMM-ONT-HQ.model")
        assert cmd[cmd.index("--difference-ratio") + 1] == "39:24:36"
        assert cmd[cmd.index("--seed") + 1] == "7"
        assert "--pass-num" not in cmd

    def test_hifi_uses_errhmm_and_multipass(self, tmp_path, model_dir):
        cmd = build_pbsim_command(
            PLATFORMS["hifi"], tmp_path / "g.fa", 10.0, tmp_path / "out", seed=7
        )
        assert cmd[cmd.index("--method") + 1] == "errhmm"
        assert cmd[cmd.index("--errhmm") + 1] == str(model_dir / "ERRHMM-SEQUEL.model")
        assert cmd[cmd.index("--difference-ratio") + 1] == "22:45:33"
        assert cmd[cmd.index("--pass-num") + 1] == "10"

    def test_depth_is_formatted_without_exponent(self, tmp_path, model_dir):
        """PBSIM3 cannot parse scientific notation."""
        cmd = build_pbsim_command(
            PLATFORMS["ont-hq"], tmp_path / "g.fa", 0.0000005, tmp_path / "out", seed=1
        )
        assert "e-" not in cmd[cmd.index("--depth") + 1]

    def test_hifi_is_flagged_as_needing_ccs(self):
        assert PLATFORMS["hifi"].needs_ccs is True
        assert PLATFORMS["ont-hq"].needs_ccs is False


def _write_maf(path, blocks):
    """Write a PBSIM3-shaped MAF: the reference 's' line is named 'ref'."""
    with gzip.open(path, "wt") as handle:
        for start, size, strand, src_size in blocks:
            handle.write("a\n")
            handle.write(f"s ref {start} {size} {strand} {src_size} ACGT\n")
            handle.write(f"s S1_1 0 {size} + {size} ACGT\n\n")


class TestMafParsing:
    """PBSIM3 names the MAF source 'ref', so contig identity is external."""

    def test_contig_name_comes_from_ref_file(self, tmp_path):
        (tmp_path / "p_0001.ref").write_text(">Chr1 some description\nACGT\n")
        assert contig_name_for(tmp_path / "p_0001.ref") == "Chr1"

    def test_non_fasta_ref_file_is_rejected(self, tmp_path):
        (tmp_path / "p_0001.ref").write_text("ACGT\n")
        with pytest.raises(ValueError, match="FASTA header"):
            contig_name_for(tmp_path / "p_0001.ref")

    def test_iter_maf_yields_forward_spans(self, tmp_path):
        maf = tmp_path / "p_0001.maf.gz"
        _write_maf(maf, [(100, 50, "+", 1000)])
        assert list(iter_maf_alignments(maf)) == [(100, 150)]

    def test_iter_maf_converts_reverse_strand_coordinates(self, tmp_path):
        maf = tmp_path / "p_0001.maf.gz"
        # On the minus strand, start is measured from the reverse strand:
        # true start = 1000 - 100 - 50.
        _write_maf(maf, [(100, 50, "-", 1000)])
        assert list(iter_maf_alignments(maf)) == [(850, 900)]

    def test_iter_maf_ignores_read_lines(self, tmp_path):
        maf = tmp_path / "p_0001.maf.gz"
        _write_maf(maf, [(0, 10, "+", 100), (20, 10, "+", 100)])
        assert list(iter_maf_alignments(maf)) == [(0, 10), (20, 30)]


class TestPooling:
    """Pooling is gzip-to-gzip: no uncompressed FASTQ is ever written."""

    def test_pooling_prefixes_names_and_counts_reads(self, tmp_path):
        source = tmp_path / "in.fq.gz"
        with gzip.open(source, "wt") as handle:
            handle.write("@S1_1\nACGT\n+\nIIII\n@S1_2 note\nTTTT\n+\nIIII\n")
        out = tmp_path / "out.fastq.gz"
        with gzip.open(out, "wt", compresslevel=GZIP_LEVEL) as dest:
            assert pool_fastq_gz(source, dest, "sampleA:baseline:h1:s42") == 2
        with gzip.open(out, "rt") as handle:
            text = handle.read()
        assert text.startswith("@sampleA:baseline:h1:s42:S1_1\n")
        assert "@sampleA:baseline:h1:s42:S1_2 note\n" in text
        assert "ACGT\n+\nIIII\n" in text

    def test_pooling_rejects_truncated_fastq(self, tmp_path):
        source = tmp_path / "in.fq.gz"
        with gzip.open(source, "wt") as handle:
            handle.write("@S1_1\nACGT\n+\n")
        out = tmp_path / "out.fastq.gz"
        with gzip.open(out, "wt", compresslevel=GZIP_LEVEL) as dest:
            with pytest.raises(ValueError, match="Truncated FASTQ"):
                pool_fastq_gz(source, dest, "p")

    def test_pooling_makes_colliding_names_unique(self, tmp_path):
        """Every component emits S1_1; the prefix is what disambiguates."""
        names = []
        for component in ("baseline", "clone40"):
            source = tmp_path / f"{component}.fq.gz"
            with gzip.open(source, "wt") as handle:
                handle.write("@S1_1\nACGT\n+\nIIII\n")
            out = tmp_path / f"{component}.out.gz"
            with gzip.open(out, "wt", compresslevel=GZIP_LEVEL) as dest:
                pool_fastq_gz(source, dest, f"cov5x_rep1:{component}:h1:s1")
            with gzip.open(out, "rt") as handle:
                names.append(handle.readline().strip())
        assert len(set(names)) == 2
