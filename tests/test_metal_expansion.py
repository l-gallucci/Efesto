"""Prevent unsupported metal assignments and preserve official-profile provenance."""
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess

import pytest

from efesto.io import filter_categories, load_registry, read_cutoffs
from efesto.operon import (filter_cluster_fegenie, filter_cluster_json,
                           _DEFAULT_OPERON_RULES, _REPORT_ALL_PATTERNS)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/validation/metal_expansion"


def hit(stem, cat="Metallophore-Synthesis"):
    return {"hmm_stem": stem, "cat": cat, "genome": "test", "orf": stem}


@pytest.mark.parametrize("core", [[], ["K01780"], ["K23120"], ["K23122"],
                                  ["K01780", "K23120"], ["K01780", "K23122"]])
def test_partial_cnt_does_not_become_siderophore_or_complete_metallophore(core):
    other = hit("copper", "MetalResistance-Copper")
    rows = [hit(stem) for stem in core] + [other]
    for result in (filter_cluster_fegenie(rows, _REPORT_ALL_PATTERNS),
                   filter_cluster_json(rows, _DEFAULT_OPERON_RULES, _REPORT_ALL_PATTERNS)):
        assert result == [other]


@pytest.mark.parametrize("optional", [[], ["K01780"]])
def test_cnt_core_accepts_pseudopaline_without_requiring_histidine_racemase(optional):
    rows = [hit(stem) for stem in ["K23120", "K23122"] + optional]
    assert filter_cluster_fegenie(rows, _REPORT_ALL_PATTERNS) == rows
    assert filter_cluster_json(rows, _DEFAULT_OPERON_RULES, _REPORT_ALL_PATTERNS) == rows


def test_cnt_core_does_not_bypass_unrelated_incomplete_fleet():
    rows = [hit("K23120"), hit("K23122"), hit("EetA", "Iron-Oxidation")]
    for result in (filter_cluster_fegenie(rows, _REPORT_ALL_PATTERNS),
                   filter_cluster_json(rows, _DEFAULT_OPERON_RULES, _REPORT_ALL_PATTERNS)):
        assert {r["hmm_stem"] for r in result} == {"K23120", "K23122"}


def test_catalog_and_all_results_keep_individual_family_evidence():
    rows = [hit("K23120")]
    assert filter_cluster_fegenie(rows, _REPORT_ALL_PATTERNS, catalog_mode=True) == rows
    assert filter_cluster_fegenie(rows, _REPORT_ALL_PATTERNS, all_results=True) == rows


def test_ln_and_w_selectors_do_not_imply_iron_or_resistance():
    cats = {cat: [("m", Path("m.hmm"))] for cat in (
        "Metallophore-Synthesis", "MetalHomeostasis-Lanthanides",
        "MetalAcquisition-Tungsten", "MetalAcquisition-MolybdenumTungsten",
        "Iron-Reduction", "Siderophore-Synthesis", "MetalResistance-Copper")}
    assert set(filter_categories(cats, ["Ln"])) == {"MetalHomeostasis-Lanthanides"}
    assert set(filter_categories(cats, ["W"])) == {"MetalAcquisition-Tungsten", "MetalAcquisition-MolybdenumTungsten"}
    assert "MetalAcquisition-MolybdenumTungsten" in filter_categories(cats, ["Mo"])
    assert "Metallophore-Synthesis" not in filter_categories(cats, ["Fe"])


def test_imoa_exclusion_and_accession_categories_survive_import():
    rows = {r["stem"]: r for r in load_registry(ROOT / "hmm_library")}
    assert rows["ImoA"]["status"].startswith("deprecated")
    assert rows["ImoA"]["category"] != "Iron-Oxidation"
    assert rows["CymA"]["status"] == "experimental"
    for rules in (_DEFAULT_OPERON_RULES, json.loads((ROOT / "operon_rules.json").read_text())["rules"]):
        assert "ImoA" not in next(r for r in rules if r["name"] == "MtrMto")["genes"]
    spec = importlib.util.spec_from_file_location("metal_importer", ROOT / "scripts/curate_hmm_library.py")
    importer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(importer)
    assert importer.infer_category("K23122", "cntM", "")[0] == "Metallophore-Synthesis"
    assert importer.infer_category("NF053806", "LanM", "")[0] == "MetalHomeostasis-Lanthanides"
    assert importer.infer_category("wtpA_TIGR03730.1", "WtpA", "")[0] == "MetalAcquisition-MolybdenumTungsten"
    model = dict(stem="ImoA", category="Iron-Oxidation", name="ImoA")
    assert importer._apply_category_curations([model])[0]["status"].startswith("deprecated")


def test_official_model_bytes_and_ga_are_preserved():
    registry = {r["stem"]: r for r in load_registry(ROOT / "hmm_library")}
    cutoffs = read_cutoffs(ROOT / "hmm_library/HMM-bitcutoffs.txt")
    for record in json.loads((DATA / "profile_sources.json").read_text())["profiles"]:
        row = registry[record["stem"]]
        raw = (ROOT / "hmm_library" / row["hmm_file"]).read_bytes()
        assert hashlib.sha256(raw).hexdigest() == record["sha256"]
        assert float(row["cutoff"]) == record["sequence_cutoff"] == cutoffs[row["stem"]]
        ga = next(line.split()[1] for line in raw.decode().splitlines() if line.startswith("GA "))
        assert float(ga) == record["sequence_cutoff"]


def test_reference_proteins_and_vup_confusers(tmp_path):
    executable = shutil.which("hmmsearch")
    if not executable:
        pytest.skip("HMMER is required for the real-sequence integration check")
    registry = {r["stem"]: r for r in load_registry(ROOT / "hmm_library")}
    checks = json.loads((DATA / "reference_expectations.json").read_text())
    for stem, expected in checks.items():
        row = registry[stem]
        output = tmp_path / (stem + ".tbl")
        subprocess.run([executable, "--cpu", "1", "--noali", "--tblout", str(output),
                        str(ROOT / "hmm_library" / row["hmm_file"]), str(DATA / "references.faa")],
                       check=True, stdout=subprocess.DEVNULL)
        passed = {line.split()[0] for line in output.read_text().splitlines()
                  if not line.startswith("#") and float(line.split()[5]) >= float(row["cutoff"])}
        assert set(expected["positive"]) <= passed
        assert not set(expected["negative"]) & passed


def test_rebuilding_library_does_not_round_official_ga(tmp_path):
    spec = importlib.util.spec_from_file_location("precise_importer", ROOT / "scripts/curate_hmm_library.py")
    importer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(importer)
    row = next(r for r in load_registry(ROOT / "hmm_library") if r["stem"] == "wtpA_TIGR03730")
    model = dict(row, hmm_file=str(ROOT / "hmm_library" / row["hmm_file"]))
    importer.write_library(tmp_path, [model], [])
    assert read_cutoffs(tmp_path / "HMM-bitcutoffs.txt")[row["stem"]] == 242.75
