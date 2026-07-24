"""Tests for baktfold.py — Foldseek results parsing, confirm_struct_hit verdicts."""

from efesto.baktfold import _parse_foldseek_results, confirm_struct_hit


def _write_m8(tmp_path, data_lines):
    p = tmp_path / "test.m8"
    p.write_text("\n".join(data_lines) + "\n")
    return str(p)


def _row(query="orf_1", target="MtrA_ref_1zzz", fident="0.9", alnlen="300",
        mismatch="10", gapopen="1", qstart="1", qend="300", tstart="1",
        tend="300", evalue="1e-50", bits="450.0"):
    return "\t".join([query, target, fident, alnlen, mismatch, gapopen,
                      qstart, qend, tstart, tend, evalue, bits])


class TestParseFoldseekResults:
    def test_parses_basic_row(self, tmp_path):
        hits = _parse_foldseek_results(_write_m8(tmp_path, [_row()]))
        assert hits["orf_1"]["target"] == "MtrA_ref_1zzz"
        assert hits["orf_1"]["bits"] == 450.0
        assert hits["orf_1"]["evalue"] == 1e-50

    def test_keeps_best_bitscore_per_query(self, tmp_path):
        lines = [_row(target="MtrA_ref_1", bits="200.0"),
                 _row(target="MtrA_ref_2", bits="500.0"),
                 _row(target="MtrA_ref_3", bits="350.0")]
        hits = _parse_foldseek_results(_write_m8(tmp_path, lines))
        assert hits["orf_1"]["target"] == "MtrA_ref_2"
        assert hits["orf_1"]["bits"] == 500.0

    def test_multiple_queries(self, tmp_path):
        lines = [_row(query="orf_1", bits="300.0"),
                 _row(query="orf_2", target="MtoA_ref_1", bits="600.0")]
        hits = _parse_foldseek_results(_write_m8(tmp_path, lines))
        assert set(hits) == {"orf_1", "orf_2"}
        assert hits["orf_2"]["target"] == "MtoA_ref_1"

    def test_missing_file_returns_empty_dict(self):
        assert _parse_foldseek_results("/nonexistent/path.m8") == {}

    def test_empty_file_returns_empty_dict(self, tmp_path):
        p = tmp_path / "empty.m8"
        p.write_text("")
        assert _parse_foldseek_results(str(p)) == {}

    def test_short_line_skipped(self, tmp_path):
        p = tmp_path / "test.m8"
        p.write_text("orf_1\tMtrA_ref\n")
        assert _parse_foldseek_results(str(p)) == {}

    def test_non_numeric_score_line_skipped(self, tmp_path):
        p = tmp_path / "test.m8"
        line = _row(bits="not_a_number")
        p.write_text(line + "\n")
        assert _parse_foldseek_results(str(p)) == {}


class TestConfirmStructHit:
    def test_matching_target_name_confirmed(self):
        row = {"target": "MtrA_ref_1zzz", "bits": 450.0, "evalue": 1e-50}
        assert confirm_struct_hit("MtrA", row) == "confirmed"

    def test_case_insensitive_match(self):
        row = {"target": "mtra_ref_1zzz", "bits": 450.0, "evalue": 1e-50}
        assert confirm_struct_hit("MtrA", row) == "confirmed"

    def test_different_target_name_contradicted(self):
        row = {"target": "MtoA_ref_2yyy", "bits": 400.0, "evalue": 1e-40}
        assert confirm_struct_hit("MtrA", row) == "contradicted"

    def test_no_foldseek_row_is_neutral(self):
        assert confirm_struct_hit("MtrA", None) == "neutral"

    def test_empty_expected_gene_name_is_neutral(self):
        row = {"target": "MtrA_ref_1zzz", "bits": 450.0, "evalue": 1e-50}
        assert confirm_struct_hit("", row) == "neutral"

    def test_empty_target_is_neutral(self):
        row = {"target": "", "bits": 0.0, "evalue": 1.0}
        assert confirm_struct_hit("MtrA", row) == "neutral"
