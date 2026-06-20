"""Tests for the CLI and module discovery."""

import argparse

from simulate_data import cli
from simulate_data.modules import fastq


def test_discover_modules_includes_fastq():
    modules = cli._discover_modules()
    assert "fastq" in modules


def test_build_parser_has_fastq_subcommand():
    parser = cli.build_parser()
    subparser_actions = [
        a for a in parser._actions if isinstance(a, argparse._SubParsersAction)
    ]
    assert subparser_actions, "No subparser action found"
    choices = subparser_actions[0].choices
    assert "fastq" in choices


def test_fastq_module_writes_expected_records(tmp_path):
    out_file = tmp_path / "output.fastq"
    ns = argparse.Namespace(
        num_reads=5,
        read_length=50,
        output=open(out_file, "w"),
        seed=42,
    )
    fastq.main(ns)

    content = out_file.read_text()
    assert content.count("@simulated_read_") == 5
    # Each FASTQ record is 4 lines
    assert len(content.strip().splitlines()) == 20


def test_fastq_reproducible_with_seed(tmp_path):
    out1 = tmp_path / "run1.fastq"
    out2 = tmp_path / "run2.fastq"

    for path in (out1, out2):
        ns = argparse.Namespace(
            num_reads=3,
            read_length=30,
            output=open(path, "w"),
            seed=99,
        )
        fastq.main(ns)

    assert out1.read_text() == out2.read_text()
