"""Tests for write_gff3, write_summary_stats, write_hit_faa, write_hit_fna."""

import csv
from pathlib import Path

import pytest

from efesto.writers import (
    write_gff3, write_hit_faa, write_hit_fna, write_summary_stats,
)


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


class TestWriteHitFaa:
    def _row(self, orf="orf1", seq="MAST", cluster_confidence=0.9):
        r = _make_row(orf=orf, cluster_confidence=cluster_confidence)
        r["sequence"] = seq
        r["bakta_id"] = orf
        return r

    def test_creates_file(self, tmp_path):
        out = str(tmp_path / "hits.faa")
        write_hit_faa(out, [self._row()])
        assert Path(out).exists()

    def test_fasta_header_fields(self, tmp_path):
        out = str(tmp_path / "hits.faa")
        write_hit_faa(out, [self._row(orf="orf1", seq="MAST")])
        lines = Path(out).read_text().splitlines()
        assert lines[0].startswith(">orf1")
        assert "gene=mtrA" in lines[0]
        assert "category=iron_reduction" in lines[0]
        assert "bitscore=200.0" in lines[0]
        assert "cluster_confidence=0.900" in lines[0]
        assert "confidence=calibrated" in lines[0]

    def test_sequence_written(self, tmp_path):
        out = str(tmp_path / "hits.faa")
        write_hit_faa(out, [self._row(seq="MKVL")])
        lines = Path(out).read_text().splitlines()
        assert lines[1] == "MKVL"

    def test_rows_without_sequence_skipped(self, tmp_path):
        out = str(tmp_path / "hits.faa")
        row = self._row(seq="")
        n = write_hit_faa(out, [row])
        assert n == 0
        assert Path(out).read_text() == ""

    def test_returns_count(self, tmp_path):
        out = str(tmp_path / "hits.faa")
        rows = [self._row("o1", "MAST"), self._row("o2", "MKLL")]
        n = write_hit_faa(out, rows)
        assert n == 2

    def test_multiple_rows(self, tmp_path):
        out = str(tmp_path / "hits.faa")
        rows = [self._row("o1", "MAST"), self._row("o2", "MKLL")]
        write_hit_faa(out, rows)
        lines = Path(out).read_text().splitlines()
        assert lines[0].startswith(">o1")
        assert lines[2].startswith(">o2")


class TestWriteHitFna:
    def _make_fna(self, tmp_path, name="g1", seq="A" * 100):
        fna = tmp_path / f"{name}.fna"
        fna.write_text(f">contig1\n{seq}\n")
        return tmp_path

    def _row_with_bakta(self, orf="orf1"):
        r = _make_row(orf=orf, genome="g1.faa", contig="contig1")
        r["sequence"] = "MAST"
        r["bakta_id"] = orf
        return r

    def test_creates_file(self, tmp_path):
        fna_dir = self._make_fna(tmp_path)
        coords  = {"g1.faa": {"orf1": _make_coords(start=1, end=9, strand="+")}}
        out     = str(tmp_path / "hits.fna")
        write_hit_fna(out, [self._row_with_bakta()], coords, str(fna_dir))
        assert Path(out).exists()

    def test_sequence_extracted(self, tmp_path):
        fna_dir = self._make_fna(tmp_path, seq="ATGCGTAAATTT" + "N" * 88)
        coords  = {"g1.faa": {"orf1": _make_coords(start=1, end=9, strand="+", contig="contig1")}}
        out     = str(tmp_path / "hits.fna")
        write_hit_fna(out, [self._row_with_bakta()], coords, str(fna_dir))
        lines = Path(out).read_text().splitlines()
        assert lines[1] == "ATGCGTAAA"

    def test_reverse_complement(self, tmp_path):
        fna_dir = self._make_fna(tmp_path, seq="ATGCGT" + "N" * 94)
        coords  = {"g1.faa": {"orf1": _make_coords(start=1, end=6, strand="-", contig="contig1")}}
        out     = str(tmp_path / "hits.fna")
        write_hit_fna(out, [self._row_with_bakta()], coords, str(fna_dir))
        lines = Path(out).read_text().splitlines()
        assert lines[1] == "ACGCAT"

    def test_header_fields(self, tmp_path):
        fna_dir = self._make_fna(tmp_path)
        coords  = {"g1.faa": {"orf1": _make_coords(start=1, end=9, strand="+", contig="contig1")}}
        out     = str(tmp_path / "hits.fna")
        write_hit_fna(out, [self._row_with_bakta()], coords, str(fna_dir))
        header = Path(out).read_text().splitlines()[0]
        assert "gene=mtrA" in header
        assert "start=1" in header
        assert "strand=+" in header
        assert "cluster_confidence=" in header

    def test_missing_fna_skipped(self, tmp_path):
        coords = {"g1.faa": {"orf1": _make_coords()}}
        out    = str(tmp_path / "hits.fna")
        n = write_hit_fna(out, [self._row_with_bakta()], coords, str(tmp_path))
        assert n == 0

    def test_orf_without_coords_skipped(self, tmp_path):
        fna_dir = self._make_fna(tmp_path)
        coords  = {"g1.faa": {}}
        out     = str(tmp_path / "hits.fna")
        n = write_hit_fna(out, [self._row_with_bakta()], coords, str(fna_dir))
        assert n == 0
