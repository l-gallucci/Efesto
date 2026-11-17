# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Efesto contributors.

"""Tests for baktfold.py — Foldseek results parsing, confirm_struct_hit verdicts."""

from efesto.baktfold import _parse_foldseek_results


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
