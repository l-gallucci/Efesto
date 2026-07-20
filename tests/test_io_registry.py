"""Tests for src/efesto/io.py — build_confirmation_map."""

import csv

from efesto.io import build_confirmation_map


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
            "sufA_x\tsufA\t\tiron_sulfur_assembly\tsufA\tncbifam\t"
            "sufA_x.hmm\t4\t136.2\t2026-07-19\tactive\tdoi:y\t\tyes\tK05997\n"
        )
        with open(p, newline="") as fh:
            registry = list(csv.DictReader(fh, delimiter="\t"))
        result = build_confirmation_map(registry)  # must not raise
        assert "old_stem" not in result
        assert "sufA_x" in result
        assert result["sufA_x"]["kegg_ko"] == ["K05997"]
