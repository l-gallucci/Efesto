"""Tests for src/efesto/io.py — build_confirmation_map, build_topology_expectation_map."""

import csv

from efesto.io import build_confirmation_map, build_topology_expectation_map


def _row(stem, needs_confirmation=None, kegg_ko=None, gene_name="gene"):
    return {"stem": stem, "gene_name": gene_name,
            "needs_confirmation": needs_confirmation, "kegg_ko": kegg_ko}


class TestBuildConfirmationMap:
    def test_flags_yes_rows(self):
        registry = [_row("sufA_x", needs_confirmation="yes", kegg_ko="K05997")]
        result = build_confirmation_map(registry)
        assert "sufA_x" in result
        assert result["sufA_x"]["kegg_ko"] == ["K05997"]

    def test_ignores_non_flagged_rows(self):
        registry = [_row("sufB_x", needs_confirmation=None, kegg_ko="K09014")]
        assert build_confirmation_map(registry) == {}

    def test_ignores_empty_string_flag(self):
        registry = [_row("sufB_x", needs_confirmation="", kegg_ko="K09014")]
        assert build_confirmation_map(registry) == {}

    def test_handles_none_kegg_ko(self):
        # a flagged model with no KO recorded should still get an entry,
        # just with an empty ko list — the escalation can still run, it
        # just won't have a KO-based signal for the comparison
        registry = [_row("sufA_x", needs_confirmation="yes", kegg_ko=None)]
        result = build_confirmation_map(registry)
        assert result["sufA_x"]["kegg_ko"] == []

    def test_comma_separated_multiple_ko(self):
        registry = [_row("foo", needs_confirmation="yes", kegg_ko="K00001,K00002")]
        result = build_confirmation_map(registry)
        assert result["foo"]["kegg_ko"] == ["K00001", "K00002"]

    def test_case_insensitive_yes(self):
        registry = [_row("foo", needs_confirmation="Yes", kegg_ko="K00001")]
        assert "foo" in build_confirmation_map(registry)

    def test_real_registry_tsv_ragged_rows_dont_crash(self, tmp_path):
        """
        Regression test: csv.DictReader fills missing trailing columns with
        None (restval=None), not "". A short row (from before needs_confirmation/
        kegg_ko columns existed) must not crash build_confirmation_map, and
        must not be treated as flagged.
        """
        p = tmp_path / "hmm_registry.tsv"
        p.write_text(
            "stem\tname\tacc\tcategory\tgene_name\tsource\thmm_file\tnseq\t"
            "cutoff\tadded_date\tstatus\treference\tvalidated_in\t"
            "needs_confirmation\tkegg_ko\n"
            # short row: only 12 fields, no needs_confirmation/kegg_ko at all
            "old_stem\told\t\tcat\told\tfegenie\told.hmm\t10\t100.0\t"
            "2026-01-01\tactive\tdoi:x\n"
            # full row with the flag set
            "sufA_x\tsufA\t\tIronSulfur-Assembly\tsufA\tncbifam\t"
            "sufA_x.hmm\t4\t136.2\t2026-07-19\tactive\tdoi:y\t\tyes\tK05997\n"
        )
        with open(p, newline="") as fh:
            registry = list(csv.DictReader(fh, delimiter="\t"))
        result = build_confirmation_map(registry)  # must not raise
        assert "old_stem" not in result
        assert "sufA_x" in result
        assert result["sufA_x"]["kegg_ko"] == ["K05997"]


def _topo_row(stem, expected_heme_c=None, expected_topology=None):
    return {"stem": stem, "expected_heme_c": expected_heme_c,
            "expected_topology": expected_topology}


class TestBuildTopologyExpectationMap:
    def test_flags_rows_with_heme_expectation(self):
        registry = [_topo_row("MtrA", expected_heme_c="10")]
        result = build_topology_expectation_map(registry)
        assert result["MtrA"] == {"expected_heme_c": 10, "expected_topology": None}

    def test_flags_rows_with_topology_expectation(self):
        registry = [_topo_row("MtrA", expected_topology="signal_peptide")]
        result = build_topology_expectation_map(registry)
        assert result["MtrA"] == {"expected_heme_c": None,
                                   "expected_topology": "signal_peptide"}

    def test_both_fields_set(self):
        registry = [_topo_row("MtrA", expected_heme_c="10",
                              expected_topology="signal_peptide")]
        result = build_topology_expectation_map(registry)
        assert result["MtrA"] == {"expected_heme_c": 10,
                                   "expected_topology": "signal_peptide"}

    def test_ignores_rows_with_neither_field(self):
        registry = [_topo_row("MtrB")]
        assert build_topology_expectation_map(registry) == {}

    def test_ignores_empty_string_fields(self):
        registry = [_topo_row("MtrB", expected_heme_c="", expected_topology="")]
        assert build_topology_expectation_map(registry) == {}

    def test_invalid_heme_c_value_becomes_none(self):
        registry = [_topo_row("MtrA", expected_heme_c="not_a_number",
                              expected_topology="signal_peptide")]
        result = build_topology_expectation_map(registry)
        assert result["MtrA"]["expected_heme_c"] is None
        assert result["MtrA"]["expected_topology"] == "signal_peptide"

    def test_real_registry_tsv_ragged_rows_dont_crash(self, tmp_path):
        """
        Regression test, same shape as
        TestBuildConfirmationMap.test_real_registry_tsv_ragged_rows_dont_crash:
        a short row from before expected_heme_c/expected_topology existed
        must not crash the loader, and must not be treated as having an
        expectation set.
        """
        p = tmp_path / "hmm_registry.tsv"
        p.write_text(
            "stem\tname\tacc\tcategory\tgene_name\tsource\thmm_file\tnseq\t"
            "cutoff\tadded_date\tstatus\treference\tvalidated_in\t"
            "needs_confirmation\tkegg_ko\texpected_heme_c\texpected_topology\n"
            # short row: only 12 fields, none of the trailing columns at all
            "old_stem\told\t\tcat\told\tfegenie\told.hmm\t10\t100.0\t"
            "2026-01-01\tactive\tdoi:x\n"
            # full row with both new fields set
            "MtrA\tMtrA\t\tIron-Reduction\tMtrA\tcurated\tMtrA.hmm\t15\t580.0\t"
            "2026-05-23\tactive\tdoi:y\t\t\t\t10\tsignal_peptide\n"
        )
        with open(p, newline="") as fh:
            registry = list(csv.DictReader(fh, delimiter="\t"))
        result = build_topology_expectation_map(registry)  # must not raise
        assert "old_stem" not in result
        assert "MtrA" in result
        assert result["MtrA"] == {"expected_heme_c": 10,
                                   "expected_topology": "signal_peptide"}
