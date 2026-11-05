"""Private exports must preserve counts, grouping, and coordinate uncertainty."""
import csv
import importlib.util
import json
from pathlib import Path

import pytest

from efesto._dev_figures import (aggregate, attach_coordinates, export_figures,
                                 read_results, select_clusters)


def write_results(path, rows):
    fields=["genome","contig","orf","gene","category","cluster_id","bitscore_cutoff"]
    with path.open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter="\t")
        writer.writeheader();writer.writerows(rows)


def hit(genome="a.faa",orf="c_1",cluster="0",cutoff="50",category="iron_reduction"):
    return dict(genome=genome,contig="c",orf=orf,gene="MtrA",category=category,
                cluster_id=cluster,bitscore_cutoff=cutoff)


def test_counts_deduplicate_exact_records_and_reject_conflicting_assignments(tmp_path):
    path=tmp_path/"Efesto-results-long.tsv"
    write_results(path,[hit(),hit(),hit(orf="c_2",cutoff="0")])
    rows=read_results(path)
    counts,categories,genomes,cutoff=aggregate(rows)
    assert len(rows)==2
    assert counts[("a.faa","Iron-Reduction")]==2
    assert cutoff=={"Model-specific cutoff":1,"Generic cutoff":1}
    write_results(path,[hit(),hit(category="Iron-Oxidation")])
    with pytest.raises(ValueError,match="Conflicting assignments"):read_results(path)


def test_cluster_zero_is_valid_and_ids_are_scoped_to_genome_and_contig(tmp_path):
    path=tmp_path/"results.tsv"
    write_results(path,[hit(),hit(orf="c_2"),hit(genome="b.faa")])
    clusters,total=select_clusters(read_results(path),6)
    assert total==1 and clusters[0][0]==("a.faa","c","0")
    assert len(clusters[0][1])==2


def test_prodigal_id_alias_and_ambiguous_coordinates(tmp_path):
    source=tmp_path/"Efesto-results-long.tsv";write_results(source,[hit()])
    gff=tmp_path/"a.gff"
    gff.write_text("c\tProdigal_v2.6.3\tCDS\t10\t60\t.\t-\t0\tID=1_1;\n")
    rows=read_results(source);used,warnings=attach_coordinates(rows,tmp_path,tmp_path)
    assert rows[0]["_coordinates"]==(10,60,"-") and not warnings
    gff.write_text(gff.read_text()+"c\tProdigal_v2.6.3\tCDS\t100\t150\t.\t+\t0\tID=1_1;\n")
    rows=read_results(source);_,warnings=attach_coordinates(rows,tmp_path,tmp_path)
    assert "_coordinates" not in rows[0] and warnings


def test_invalid_nan_cutoff_is_rejected(tmp_path):
    source=tmp_path/"r.tsv";write_results(source,[hit(cutoff="nan")])
    with pytest.raises(ValueError,match="Invalid cutoff"):read_results(source)


@pytest.mark.skipif(importlib.util.find_spec("matplotlib") is None,reason="optional plotting dependency")
def test_export_real_coordinates_or_explicit_schematic_and_empty_results(tmp_path):
    source=tmp_path/"input";source.mkdir()
    write_results(source/"Efesto-results-long.tsv",[hit(),hit(orf="c_2")])
    out=tmp_path/"figures"
    result=export_figures(source,out,max_clusters=1,dpi=60)
    assert len(result["outputs"])==12
    assert result["cluster_selection"][0]["genomic_coordinates"] is False
    assert (out/"01_functional_overview.png").stat().st_size>1000
    assert json.loads((out/"figure_manifest.json").read_text())["orf_count"]==2
    write_results(source/"Efesto-results-long.tsv",[])
    (source/"Efesto-heatmap-data.csv").write_text("X,a.faa,b.faa\n")
    result=export_figures(source,tmp_path/"empty",dpi=60)
    assert result["orf_count"]==0 and len(result["outputs"])==9
    assert result["total_samples"]==2 and result["samples_with_hits"]==0
