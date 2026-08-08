"""Tests for the shared long-read simulation helpers."""

import gzip

import pytest

from simulate_data.longread import resolve_model_path


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
