"""Private comparisons must preserve source models and stay off by default."""
import argparse
import hashlib
import sys
from pathlib import Path

import pytest

from efesto._dev import add_dev_arguments, load_fegenie_catalog
from efesto import cli


def model(root, category, stem, *, date="one", body="same weights"):
    path = root / category / (stem + ".hmm")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"HMMER3/f\nNAME  {stem}\nDATE  {date}\nNSEQ  10\n{body}\n//\n")
    return path


@pytest.fixture
def library(tmp_path):
    (tmp_path / "HMM-bitcutoffs.txt").write_text("HMM\tcutoff\nA\t31.9\nB\t88\n")
    (tmp_path / "FeGenie-map.txt").write_text("A\tgeneA\nB\tgeneB\n")
    model(tmp_path, "iron_aquisition-iron_transport", "A")
    model(tmp_path, "iron_reduction", "B.txt-nr")
    model(tmp_path, "iron_reduction", "NoMetadata")
    return tmp_path


def test_original_models_cutoff_formatting_and_first_dot_lookup(library):
    before = {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in library.rglob("*") if p.is_file()}
    cats, cuts, rows, manifest = load_fegenie_catalog(library)
    assert cuts == {"A": 31, "B.txt-nr": 88}
    assert cats["IronAcquisition-Transport"][0][0] == "A"
    assert cats["Iron-Reduction"][0][1] == library / "iron_reduction/B.txt-nr.hmm"
    assert manifest["skipped_missing_metadata"] == ["iron_reduction/NoMetadata.hmm"]
    assert manifest["models"][0]["declared_cutoff"] == 31.9
    assert rows[0]["source"] == "fegenie"
    assert before == {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                      for p in library.rglob("*") if p.is_file()}


def test_date_only_duplicates_search_once_and_record_both_paths(library):
    model(library, "iron_oxidation", "A", date="two")
    cats, _, _, manifest = load_fegenie_catalog(library)
    assert sum(map(len, cats.values())) == 2
    assert manifest["metadata_matched_paths"] == 3
    assert len(manifest["models"][0]["source_files"]) == 2


def test_different_profiles_with_same_stem_fail_instead_of_colliding(library):
    model(library, "iron_oxidation", "A", body="different weights")
    with pytest.raises(ValueError, match="share stem A"):
        load_fegenie_catalog(library)


@pytest.mark.parametrize("value", ["nan", "inf", "0", "-20"])
def test_invalid_metadata_does_not_silently_fall_back(library, value):
    (library / "HMM-bitcutoffs.txt").write_text("A\t" + value + "\n")
    with pytest.raises(ValueError, match="Invalid positive"):
        load_fegenie_catalog(library)


def test_empty_or_wrong_library_fails(tmp_path):
    with pytest.raises(ValueError, match="original FeGenie"):
        load_fegenie_catalog(tmp_path)


def test_flag_hidden_even_in_developer_help():
    parser = argparse.ArgumentParser()
    add_dev_arguments(parser)
    assert "dev-fegenie" not in parser.format_help()
    assert parser.parse_args(["--dev-fegenie-hmms", "/tmp/db"])._dev_fegenie_hmms == Path("/tmp/db")


def test_normal_cli_rejects_developer_flag(monkeypatch, capsys):
    monkeypatch.delenv("EFESTO_DEV_MODE", raising=False)
    monkeypatch.setattr(sys, "argv", ["efesto", "--faa_dir", "/unused",
                                    "--dev-fegenie-hmms", "/unused"])
    with pytest.raises(SystemExit) as exc:
        cli.main()
    assert exc.value.code == 2
    assert "unrecognized arguments" in capsys.readouterr().err


@pytest.mark.parametrize("extra", [["--normalize_hmms"], ["--annotate", "Fe"],
                                  ["--hmm_dir", "/unused"]])
def test_dev_cli_rejects_overrides_and_source_mutation(monkeypatch, capsys, extra):
    monkeypatch.setenv("EFESTO_DEV_MODE", "1")
    monkeypatch.setattr(sys, "argv", ["efesto", "--faa_dir", "/unused",
                                    "--dev-fegenie-hmms", "/unused", *extra])
    with pytest.raises(SystemExit) as exc:
        cli.main()
    assert exc.value.code == 2
    assert "cannot combine" in capsys.readouterr().err


def test_dev_cli_rejects_old_output_cache(tmp_path, monkeypatch, capsys):
    (tmp_path / "old.tblout").write_text("old data")
    monkeypatch.setenv("EFESTO_DEV_MODE", "1")
    monkeypatch.setattr(sys, "argv", ["efesto", "--faa_dir", "/unused",
                                    "--dev-fegenie-hmms", "/unused",
                                    "--out", str(tmp_path)])
    with pytest.raises(SystemExit) as exc:
        cli.main()
    assert exc.value.code == 2
    assert "fresh or empty" in capsys.readouterr().err
