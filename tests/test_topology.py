"""Tests for topology.py — Phobius output parsing and verdict logic."""

from efesto.topology import (_parse_phobius_output, heme_verdict,
                             topology_verdict)


def _write_phobius(tmp_path, data_lines):
    p = tmp_path / "test.phobius"
    lines = ["SEQENCE ID  TM SP PREDICTION"] + data_lines
    p.write_text("\n".join(lines) + "\n")
    return str(p)


class TestParsePhobiusOutput:
    def test_parses_basic_row(self, tmp_path):
        line = "orf_1  0  Y  n8-18c19/20o"
        hits = _parse_phobius_output(_write_phobius(tmp_path, [line]))
        assert hits["orf_1"] == {"signal_peptide": True, "tm_helices": 0}

    def test_parses_multi_pass_tm_row(self, tmp_path):
        line = "orf_2  6  0  i7-29o48-70i"
        hits = _parse_phobius_output(_write_phobius(tmp_path, [line]))
        assert hits["orf_2"] == {"signal_peptide": False, "tm_helices": 6}

    def test_skips_header_line(self, tmp_path):
        hits = _parse_phobius_output(_write_phobius(tmp_path, ["orf_1  0  0  o"]))
        assert len(hits) == 1

    def test_missing_file_returns_empty_dict(self):
        assert _parse_phobius_output("/nonexistent/path.phobius") == {}

    def test_empty_file_returns_empty_dict(self, tmp_path):
        p = tmp_path / "empty.phobius"
        p.write_text("")
        assert _parse_phobius_output(str(p)) == {}

    def test_short_line_skipped(self, tmp_path):
        hits = _parse_phobius_output(_write_phobius(tmp_path, ["orf_1  0"]))
        assert hits == {}


class TestTopologyVerdict:
    def test_no_expectation_is_neutral(self):
        row = {"signal_peptide": True, "tm_helices": 0}
        assert topology_verdict(None, row) == "neutral"

    def test_no_phobius_row_is_neutral(self):
        assert topology_verdict("signal_peptide", None) == "neutral"

    def test_signal_peptide_confirmed(self):
        row = {"signal_peptide": True, "tm_helices": 0}
        assert topology_verdict("signal_peptide", row) == "confirmed"

    def test_signal_peptide_contradicted_by_multi_pass_tm(self):
        row = {"signal_peptide": False, "tm_helices": 6}
        assert topology_verdict("signal_peptide", row) == "contradicted"

    def test_signal_peptide_borderline_single_tm_is_neutral(self):
        row = {"signal_peptide": False, "tm_helices": 1}
        assert topology_verdict("signal_peptide", row) == "neutral"

    def test_multi_pass_tm_confirmed(self):
        row = {"signal_peptide": False, "tm_helices": 4}
        assert topology_verdict("multi_pass_tm", row) == "confirmed"

    def test_multi_pass_tm_contradicted_by_soluble(self):
        row = {"signal_peptide": False, "tm_helices": 0}
        assert topology_verdict("multi_pass_tm", row) == "contradicted"

    def test_cytoplasmic_confirmed(self):
        row = {"signal_peptide": False, "tm_helices": 0}
        assert topology_verdict("cytoplasmic", row) == "confirmed"

    def test_cytoplasmic_contradicted_by_signal_peptide(self):
        row = {"signal_peptide": True, "tm_helices": 0}
        assert topology_verdict("cytoplasmic", row) == "contradicted"

    def test_unknown_expected_topology_is_neutral(self):
        row = {"signal_peptide": True, "tm_helices": 0}
        assert topology_verdict("some_future_category", row) == "neutral"


class TestHemeVerdict:
    def test_no_expectation_is_neutral(self):
        assert heme_verdict(None, 10) == "neutral"

    def test_no_observed_value_is_neutral(self):
        assert heme_verdict(10, None) == "neutral"

    def test_matching_count_confirmed(self):
        assert heme_verdict(10, 10) == "confirmed"

    def test_within_tolerance_confirmed(self):
        assert heme_verdict(10, 8) == "confirmed"  # 8 >= 0.7*10

    def test_zero_observed_against_low_expectation_is_neutral(self):
        # expected < 4: a single missed motif shouldn't read as a strong signal
        assert heme_verdict(2, 0) == "neutral"

    def test_zero_observed_against_high_expectation_contradicted(self):
        assert heme_verdict(10, 0) == "contradicted"

    def test_below_tolerance_but_nonzero_is_neutral(self):
        assert heme_verdict(10, 3) == "neutral"
