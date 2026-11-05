"""Category migration must preserve profiles and keep filtering functional."""
import csv
import importlib.util
import json
from pathlib import Path

import pytest

from efesto.categories import canonical_category
from efesto.io import filter_categories, load_registry
from efesto.operon import (filter_cluster_fegenie, filter_cluster_json,
                          load_operon_rules, second_pass, _REPORT_ALL_PATTERNS)

ROOT = Path(__file__).resolve().parents[1]


def load_script(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def hit(stem, category):
    return {"hmm_stem": stem, "cat": category, "genome": "g", "orf": stem}


def test_legacy_library_names_and_canonical_selectors():
    catalog = {"iron_aquisition-iron_transport": [("a", Path("a.hmm"))],
               "MetalAcquisition-Nickel": [("b", Path("b.hmm"))],
               "MetalResistance-Nickel": [("c", Path("c.hmm"))]}
    assert set(filter_categories(catalog, ["Fe-acquisition"])) == {"IronAcquisition-Transport"}
    assert set(filter_categories(catalog, ["Ni"])) == {"MetalAcquisition-Nickel", "MetalResistance-Nickel"}
    assert set(filter_categories(catalog, ["Ni-resistance"])) == {"MetalResistance-Nickel"}
    assert set(filter_categories(catalog, ["MetalAcquisition-Nickel"])) == {"MetalAcquisition-Nickel"}
    assert canonical_category("Custom-Function") == "Custom-Function"


def test_new_categories_do_not_disable_context_filters():
    singleton = [hit("biosynthesis", "Siderophore-Synthesis")]
    assert filter_cluster_fegenie(singleton, _REPORT_ALL_PATTERNS) == []
    complete = [hit(stem, "Siderophore-Synthesis") for stem in ("a", "b", "c")]
    assert filter_cluster_fegenie(complete, _REPORT_ALL_PATTERNS) == complete
    acquisition = [hit("ModC", "MetalAcquisition-Molybdenum")]
    assert filter_cluster_fegenie(acquisition, _REPORT_ALL_PATTERNS) == acquisition
    mixed_transport = [hit("receptor", "IronAcquisition-ComplexTransport")]
    assert filter_cluster_fegenie(mixed_transport, _REPORT_ALL_PATTERNS) == []


def test_original_json_category_names_still_apply(tmp_path):
    data = {"report_all_categories": ["metal_resistance-*"], "rules": [{
        "name": "test", "categories": ["iron_acquisition-siderophore_synthesis"],
        "genes": [], "rule": "require_n_cat", "min_genes": 3, "on_fail": "drop"}]}
    (tmp_path / "operon_rules.json").write_text(json.dumps(data))
    rules, patterns, present = load_operon_rules(tmp_path)
    assert present and patterns == ["MetalResistance-*"]
    assert filter_cluster_json([hit("single", "Siderophore-Synthesis")], rules, patterns) == []


def test_redox_disambiguation_and_regulator_second_pass():
    rows = [hit("MtoA", "Iron-Oxidation"), hit("MtrB_TIGR03509", "Iron-Reduction")]
    result = filter_cluster_fegenie(rows, _REPORT_ALL_PATTERNS)
    assert {r["cat"] for r in result} == {"Iron-Oxidation"}
    regulator = [hit("Fur", "Iron-Regulation")]
    assert second_pass(regulator, {"Fur": "Iron-Regulation"}, {}) == regulator


def test_library_paths_curations_and_exclusions_are_consistent():
    library = ROOT / "hmm_library"
    rows = load_registry(library)
    assert rows
    assert len({r["stem"] for r in rows}) == len(rows)
    for row in rows:
        expected = Path(row["category"]) / (row["stem"] + ".hmm")
        assert Path(row["hmm_file"]) == expected
        assert (library / expected).is_file()
        assert canonical_category(row["category"]) == row["category"]
        assert row["category_review_status"]
    by_stem = {r["stem"]: r for r in rows}
    assert by_stem["K04759"]["category"] == "IronAcquisition-Transport"
    assert by_stem["K02361"]["category"] == "Siderophore-Synthesis"
    assert by_stem["NF040986"]["category"] == "Magnetosome-Formation"
    assert by_stem["modC_TIGR02142"]["category"] == "MetalAcquisition-Molybdenum"
    for stem in ("NF041018", "NF041019", "NF041020"):
        assert stem not in by_stem
        assert not list(library.glob("*/" + stem + ".hmm"))
    assert by_stem["NF040982"]["status"] == "deprecated_out_of_scope"
    assert by_stem["Sid_FpvI_regulator_PA2387_Paeruginosa_PAO1_180620"]["status"] == "deprecated_misassigned_duplicate"
    assert by_stem["Sid_FpvH_iron_release_PA2404_Paeruginosa_PAO1_180602"]["status"] == "active"
    assert by_stem["NF040984"]["status"].startswith("deprecated")
    assert by_stem["As_lyase_1"]["category_review_status"] == "requires_review"


def test_importer_does_not_reintroduce_old_assignments_or_reactivate_exclusions(tmp_path):
    importer = load_script("curate_hmm_library")
    assert importer.infer_category("K04759", "feoB", "")[0] == "IronAcquisition-Transport"
    assert importer.infer_category("K02361", "entC", "")[0] == "Siderophore-Synthesis"
    assert importer.infer_category("NF040986", "MamH", "")[0] == "Magnetosome-Formation"
    model = dict(stem="NF040982", category="siderophore_transport", name="ComGD", acc="NF040982",
                 nseq=10, cutoff=99, source="interpro", hmm_file=str(tmp_path / "original.hmm"))
    original = b"HMMER3/f\nNAME ComGD\nNSEQ 10\n//\n"
    Path(model["hmm_file"]).write_bytes(original)
    importer.write_library(tmp_path / "out", importer._apply_category_curations([model]), [])
    with (tmp_path / "out/hmm_registry.tsv").open() as fh:
        row = next(csv.DictReader(fh, delimiter="\t"))
    assert row["status"] == "deprecated_out_of_scope"
    assert row["hmm_file"] == "Competence-Pilus/NF040982.hmm"
    assert (tmp_path / "out" / row["hmm_file"]).read_bytes() == original


def test_benchmark_normalizes_labels_without_recategorizing_biology(tmp_path):
    benchmark = load_script("benchmark_fegenie")
    (tmp_path / "Efesto-results-long.tsv").write_text(
        "genome\tcategory\tgene\torf\tbitscore\n"
        "g.faa\tiron_reduction\tMtrA\tc_1\t100\n"
        "g.faa\tiron_storage\tTbpB\tc_2\t200\n")
    hits = benchmark.parse_efesto(tmp_path)
    assert hits[0]["category"] == "Iron-Reduction"
    assert hits[1]["category"] == "Iron-Storage"  # Historical assignment remains visible.


def test_migration_preflights_collisions_preserves_bytes_and_is_idempotent(tmp_path):
    migration = load_script("migrate_categories")
    old = tmp_path / "iron_storage/model.hmm"
    old.parent.mkdir();old.write_bytes(b"original HMM bytes\n")
    (tmp_path / "hmm_registry.tsv").write_text(
        "stem\tcategory\thmm_file\tstatus\nmodel\tiron_storage\tiron_storage/model.hmm\tactive\n")
    (tmp_path / "category_curation.tsv").write_text(
        "stem\tcategory\tstatus_override\treview_status\tnote\tevidence\n")
    target = tmp_path / "Iron-Storage/model.hmm"
    target.parent.mkdir();target.write_text("collision")
    with pytest.raises(ValueError, match="Destination already exists"):
        migration.migrate(tmp_path, apply=True)
    assert old.read_bytes() == b"original HMM bytes\n"
    target.unlink()
    assert migration.migrate(tmp_path)["file_moves"] == 1
    assert old.exists()
    migration.migrate(tmp_path, apply=True)
    assert target.read_bytes() == b"original HMM bytes\n"
    assert migration.migrate(tmp_path)["file_moves"] == 0


def test_removed_carbon_profiles_cannot_be_reimported_even_from_reviewed_report(tmp_path):
    importer = load_script("curate_hmm_library")
    models = [dict(stem=stem, category="CarbonMetabolism-Dehydrogenase", name=stem,
                   acc=stem, nseq=15, cutoff=100, source="interpro", hmm_file="missing.hmm")
              for stem in ("NF041018", "NF041019", "NF041020")]
    assert importer._apply_category_curations(models) == []
    importer.write_library(tmp_path / "out", models, [])
    with (tmp_path / "out/hmm_registry.tsv").open() as handle:
        assert list(csv.DictReader(handle, delimiter="\t")) == []
    assert not list((tmp_path / "out").glob("*/*.hmm"))
