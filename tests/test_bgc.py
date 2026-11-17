# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Efesto contributors.

"""Tests for src/efesto/bgc.py"""

from pathlib import Path

import pytest

from efesto.bgc import (
    _extract_bgc_type,
    _is_siderophore_type,
    parse_antismash_gff,
)


def _make_gff(tmp_path, stem, regions):
    """Write a minimal antiSMASH-style GFF3 with given regions list."""
    gff = tmp_path / f"{stem}.gff3"
    lines = ["##gff-version 3\n"]
    for contig, start, end, bgc_type in regions:
        lines.append(
            f"{contig}\tantismash\tregion\t{start}\t{end}\t.\t.\t.\t"
            f"ID={contig}.region001;product={bgc_type}\n"
        )
    gff.write_text("".join(lines))
    return gff


class FakeFaa:
    """Minimal stand-in for a Path object with .stem and .name."""
    def __init__(self, stem):
        self.stem = stem
        self.name = stem + ".faa"


class TestIsSiderophoreType:
    def test_siderophore_literal(self):
        assert _is_siderophore_type("siderophore")

    def test_ni_siderophore(self):
        assert _is_siderophore_type("NI-siderophore")

    def test_metallophore(self):
        assert _is_siderophore_type("NRP-metallophore")

    def test_nrps_not_siderophore(self):
        assert not _is_siderophore_type("NRPS")

    def test_terpene_not_siderophore(self):
        assert not _is_siderophore_type("terpene")

    def test_case_insensitive(self):
        assert _is_siderophore_type("SIDEROPHORE")


class TestExtractBgcType:
    def test_product_attr(self):
        assert _extract_bgc_type("ID=r1;product=siderophore;contig_edge=False") == "siderophore"

    def test_type_attr_fallback(self):
        assert _extract_bgc_type("ID=r1;type=NRPS") == "NRPS"

    def test_product_wins_over_type(self):
        assert _extract_bgc_type("product=siderophore;type=NRPS") == "siderophore"

    def test_no_attr_returns_empty(self):
        assert _extract_bgc_type("ID=r1;score=.") == ""


class TestParseAntismashGff:
    def test_returns_empty_when_bgc_dir_none(self):
        assert parse_antismash_gff(None, []) == {}

    def test_finds_siderophore_region(self, tmp_path):
        stem = "genome1"
        _make_gff(tmp_path, stem, [("contig1", 1000, 50000, "siderophore")])
        result = parse_antismash_gff(str(tmp_path), [FakeFaa(stem)])
        assert "genome1.faa" in result
        assert result["genome1.faa"][0][3] == "siderophore"

    def test_ignores_non_siderophore(self, tmp_path):
        stem = "genome1"
        _make_gff(tmp_path, stem, [("contig1", 1000, 50000, "NRPS")])
        result = parse_antismash_gff(str(tmp_path), [FakeFaa(stem)])
        assert "genome1.faa" not in result

    def test_mixed_regions_keeps_only_siderophore(self, tmp_path):
        stem = "genome1"
        _make_gff(tmp_path, stem, [
            ("contig1", 1000, 20000, "terpene"),
            ("contig2", 5000, 30000, "siderophore"),
        ])
        result = parse_antismash_gff(str(tmp_path), [FakeFaa(stem)])
        assert len(result["genome1.faa"]) == 1
        assert result["genome1.faa"][0][0] == "contig2"

    def test_no_gff_file_no_entry(self, tmp_path):
        result = parse_antismash_gff(str(tmp_path), [FakeFaa("missing_genome")])
        assert result == {}

    def test_nested_directory_lookup(self, tmp_path):
        stem = "genome2"
        nested = tmp_path / stem
        nested.mkdir()
        _make_gff(nested, stem, [("c1", 100, 5000, "NI-siderophore")])
        result = parse_antismash_gff(str(tmp_path), [FakeFaa(stem)])
        assert "genome2.faa" in result


def test_status_distinguishes_missing_input_from_no_qualifying_regions(tmp_path):
    _make_gff(tmp_path, "present", [("ctg", 1, 1000, "NRPS")])
    regions, statuses = parse_antismash_gff(tmp_path, [FakeFaa("present"), FakeFaa("missing")], return_status=True)
    assert regions == {}
    assert statuses == {"present.faa": "parsed", "missing.faa": "input_unavailable"}
