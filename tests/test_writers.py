"""Tests for write_gff3 and write_summary_stats in writers.py."""

import csv
from pathlib import Path

import pytest

from metalgenie_evo.writers import write_gff3, write_summary_stats


def _make_row(orf="orf1", genome="g1.faa", contig="c1", cat="iron_reduction",
              gene_name="mtrA", hmm_stem="MtrA", cluster_id=0,
              cluster_confidence=0.85, confidence="calibrated",
              bitscore=200.0, evalue=1e-50, heme_motifs=0,
              cutoff=140.0, contig_len=50000):
    return {
        "orf":                orf,
        "genome":             genome,
        "contig":             contig,
        "cat":                cat,
        "gene_name":          gene_name,
        "hmm_stem":           hmm_stem,
        "cluster_id":         cluster_id,
        "cluster_confidence": cluster_confidence,
        "confidence":         confidence,
        "bitscore":           bitscore,
        "evalue":             evalue,
        "heme_motifs":        heme_motifs,
        "cutoff":             cutoff,
        "contig_len":         contig_len,
        "sequence":           "",
    }


def _make_coords(orf="orf1", start=5000, end=5300, strand="+", contig="c1"):
    return {"start": start, "end": end, "strand": strand, "contig": contig}


class TestWriteGff3:
    def test_creates_file(self, tmp_path):
        rows   = [_make_row()]
        coords = {"g1.faa": {"orf1": _make_coords()}}
        out    = str(tmp_path / "out.gff3")
        write_gff3(out, rows, coords)
        assert Path(out).exists()

    def test_gff3_header(self, tmp_path):
        rows   = [_make_row()]
        coords = {"g1.faa": {"orf1": _make_coords()}}
        out    = str(tmp_path / "out.gff3")
        write_gff3(out, rows, coords)
        header = Path(out).read_text().splitlines()[0]
        assert header == "##gff-version 3"

    def test_feature_columns(self, tmp_path):
        rows   = [_make_row()]
        coords = {"g1.faa": {"orf1": _make_coords(start=100, end=400, strand="+")}}
        out    = str(tmp_path / "out.gff3")
        write_gff3(out, rows, coords)
        lines = [l for l in Path(out).read_text().splitlines() if not l.startswith("#")]
        assert len(lines) == 1
        parts = lines[0].split("\t")
        assert parts[0] == "c1"
        assert parts[2] == "CDS"
        assert parts[3] == "100"
        assert parts[4] == "400"
        assert parts[6] == "+"

    def test_attributes_contain_gene(self, tmp_path):
        rows   = [_make_row(gene_name="mtrA")]
        coords = {"g1.faa": {"orf1": _make_coords()}}
        out    = str(tmp_path / "out.gff3")
        write_gff3(out, rows, coords)
        content = Path(out).read_text()
        assert "gene=mtrA" in content

    def test_cluster_confidence_in_attributes(self, tmp_path):
        rows   = [_make_row(cluster_confidence=0.750)]
        coords = {"g1.faa": {"orf1": _make_coords()}}
        out    = str(tmp_path / "out.gff3")
        write_gff3(out, rows, coords)
        assert "cluster_confidence=0.750" in Path(out).read_text()

    def test_row_without_coords_skipped(self, tmp_path):
        rows   = [_make_row(orf="orf_no_coord")]
        coords = {"g1.faa": {}}
        out    = str(tmp_path / "out.gff3")
        n = write_gff3(out, rows, coords)
        assert n == 0
        data_lines = [l for l in Path(out).read_text().splitlines()
                      if l and not l.startswith("#")]
        assert len(data_lines) == 0

    def test_returns_count(self, tmp_path):
        rows = [_make_row("orf1"), _make_row("orf2")]
        coords = {
            "g1.faa": {
                "orf1": _make_coords("orf1"),
                "orf2": _make_coords("orf2", start=6000, end=6300),
            }
        }
        out = str(tmp_path / "out.gff3")
        n = write_gff3(out, rows, coords)
        assert n == 2

    def test_empty_rows(self, tmp_path):
        out = str(tmp_path / "out.gff3")
        n = write_gff3(out, [], {})
        assert n == 0


class TestWriteSummaryStats:
    def _rows(self, n=3):
        return [_make_row(orf=f"orf{i}", cluster_id=i,
                          cluster_confidence=0.9 - i * 0.2) for i in range(n)]

    def test_creates_file(self, tmp_path):
        out = str(tmp_path / "stats.tsv")
        write_summary_stats(out, self._rows(), ["g1.faa"])
        assert Path(out).exists()

    def test_has_tsv_header(self, tmp_path):
        out = str(tmp_path / "stats.tsv")
        write_summary_stats(out, self._rows(), ["g1.faa"])
        lines = [l for l in Path(out).read_text().splitlines()
                 if not l.startswith("#")]
        assert lines[0].startswith("section\tmetric\tvalue")

    def test_run_section_present(self, tmp_path):
        out = str(tmp_path / "stats.tsv")
        write_summary_stats(out, self._rows(), ["g1.faa"])
        content = Path(out).read_text()
        assert "RUN\ttotal_orf_hits" in content

    def test_confidence_section(self, tmp_path):
        out = str(tmp_path / "stats.tsv")
        write_summary_stats(out, self._rows(), ["g1.faa"])
        content = Path(out).read_text()
        assert "CONFIDENCE" in content
        assert "high_confidence_clusters" in content

    def test_category_section(self, tmp_path):
        out = str(tmp_path / "stats.tsv")
        write_summary_stats(out, self._rows(), ["g1.faa"])
        assert "CATEGORY" in Path(out).read_text()

    def test_genome_section_all_genomes(self, tmp_path):
        out = str(tmp_path / "stats.tsv")
        write_summary_stats(out, self._rows(), ["g1.faa", "g2.faa"])
        content = Path(out).read_text()
        assert "g2.faa" in content

    def test_runtime_written_when_provided(self, tmp_path):
        out = str(tmp_path / "stats.tsv")
        write_summary_stats(out, self._rows(), ["g1.faa"], runtime_s=42.5)
        assert "runtime_seconds" in Path(out).read_text()

    def test_empty_rows(self, tmp_path):
        out = str(tmp_path / "stats.tsv")
        write_summary_stats(out, [], ["g1.faa"])
        assert Path(out).exists()
