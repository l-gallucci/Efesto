"""Tests for eggnog.py — annotations parsing, confirm_hit verdicts, FASTA export."""

import pytest

from efesto.eggnog import (_parse_eggnog_annotations, confirm_hit,
                            write_flagged_faa)

_HEADER = ("#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\t"
           "COG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\t"
           "KEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\t"
           "KEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs\tannotation_confidence")


def _write_annotations(tmp_path, data_lines, header=_HEADER):
    p = tmp_path / "test.emapper.annotations"
    lines = ["## Fri Jul 20 12:00:00 2026", "## emapper.py test", header]
    lines += data_lines
    lines.append("## N queries scanned")
    p.write_text("\n".join(lines) + "\n")
    return str(p)


class TestParseEggnogAnnotations:
    def test_parses_basic_row(self, tmp_path):
        line = ("orf_1\t511145.b1684\t1e-80\t280.0\tCOG0316@1|root\t1|root\tP\t"
                "Fe-S cluster assembly protein\tsufA\t-\t-\tko:K05997\t-\t-\t-\t"
                "-\t-\t-\t-\t-\tFe-S_biosyn\t1.0")
        hits = _parse_eggnog_annotations(_write_annotations(tmp_path, [line]))
        assert "orf_1" in hits
        assert hits["orf_1"]["preferred_name"] == "sufA"
        assert hits["orf_1"]["kegg_ko"] == ["K05997"]

    def test_skips_comment_lines(self, tmp_path):
        line = ("orf_1\t-\t-\t-\t-\t-\t-\t-\tsufA\t-\t-\tko:K05997\t-\t-\t-\t"
                "-\t-\t-\t-\t-\t-\t1.0")
        hits = _parse_eggnog_annotations(_write_annotations(tmp_path, [line]))
        assert len(hits) == 1

    def test_dash_placeholder_normalized_to_empty(self, tmp_path):
        line = "orf_1\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-"
        hits = _parse_eggnog_annotations(_write_annotations(tmp_path, [line]))
        assert hits["orf_1"]["preferred_name"] == ""
        assert hits["orf_1"]["kegg_ko"] == []

    def test_multiple_ko_comma_separated(self, tmp_path):
        line = ("orf_1\t-\t-\t-\t-\t-\t-\t-\tfoo\t-\t-\tko:K00001,ko:K00002\t"
                "-\t-\t-\t-\t-\t-\t-\t-\t-\t1.0")
        hits = _parse_eggnog_annotations(_write_annotations(tmp_path, [line]))
        assert hits["orf_1"]["kegg_ko"] == ["K00001", "K00002"]

    def test_missing_file_returns_empty_dict(self):
        assert _parse_eggnog_annotations("/nonexistent/path.annotations") == {}

    def test_empty_file_returns_empty_dict(self, tmp_path):
        p = tmp_path / "empty.annotations"
        p.write_text("")
        assert _parse_eggnog_annotations(str(p)) == {}


class TestConfirmHit:
    def test_matching_ko_is_confirmed(self):
        row = {"preferred_name": "sufA", "kegg_ko": ["K05997"], "description": ""}
        assert confirm_hit("sufA", ["K05997"], row) == "confirmed"

    def test_different_ko_is_contradicted(self):
        row = {"preferred_name": "csdE", "kegg_ko": ["K05978"], "description": ""}
        assert confirm_hit("sufE", ["K02426"], row) == "contradicted"

    def test_no_eggnog_row_is_neutral(self):
        assert confirm_hit("sufA", ["K05997"], None) == "neutral"

    def test_uninformative_row_is_neutral(self):
        row = {"preferred_name": "", "kegg_ko": [], "description": ""}
        assert confirm_hit("sufA", ["K05997"], row) == "neutral"

    def test_no_expected_ko_falls_back_to_name(self):
        row = {"preferred_name": "sufa", "kegg_ko": [], "description": ""}
        assert confirm_hit("sufA", [], row) == "confirmed"

    def test_no_expected_ko_name_mismatch_contradicted(self):
        row = {"preferred_name": "csde", "kegg_ko": [], "description": ""}
        assert confirm_hit("sufA", [], row) == "contradicted"

    def test_overlapping_ko_lists_confirmed(self):
        # HMM-side expects one of several possible KOs; any overlap confirms.
        row = {"preferred_name": "", "kegg_ko": ["K00002", "K05997"], "description": ""}
        assert confirm_hit("sufA", ["K05997"], row) == "confirmed"


class TestWriteFlaggedFaa:
    def _write_faa(self, tmp_path):
        p = tmp_path / "genome.faa"
        p.write_text(">orf_1 desc one\nMSEQAAAA\n>orf_2 desc two\nMSEQBBBB\n"
                     ">orf_3 desc three\nMSEQCCCC\n")
        return str(p)

    def test_writes_only_flagged_orfs(self, tmp_path):
        faa = self._write_faa(tmp_path)
        out = str(tmp_path / "flagged.faa")
        n = write_flagged_faa({"orf_1", "orf_3"}, faa, out)
        assert n == 2
        content = open(out).read()
        assert ">orf_1" in content
        assert ">orf_2" not in content
        assert ">orf_3" in content

    def test_empty_flagged_set_writes_empty_file(self, tmp_path):
        faa = self._write_faa(tmp_path)
        out = str(tmp_path / "flagged.faa")
        n = write_flagged_faa(set(), faa, out)
        assert n == 0
        assert open(out).read() == ""
