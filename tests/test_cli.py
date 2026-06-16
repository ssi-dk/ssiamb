from __future__ import annotations

from pathlib import Path

from ssiamb.cli import main
from ssiamb.pipeline import output_paths, safe_sample_name


def test_help(capsys) -> None:
    try:
        main(["--help"])
    except SystemExit as exc:
        assert exc.code == 0
    out = capsys.readouterr().out
    assert "ssiamb" in out
    assert "self" in out


def test_safe_sample_name() -> None:
    assert safe_sample_name("sample 1/a") == "sample_1_a"


def test_output_paths() -> None:
    paths = output_paths(Path("out"), "sample")
    assert paths.summary == Path("out/sample.summary.tsv")
    assert paths.vcf == Path("out/sample.variants.vcf")


def test_output_paths_keep_dotted_sample_name() -> None:
    paths = output_paths(Path("out"), "sample.v1")
    assert paths.summary == Path("out/sample.v1.summary.tsv")
