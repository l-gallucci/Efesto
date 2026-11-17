# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Efesto contributors.

"""Repeated ORF and Bakta labels must remain distinct across genomes."""
import csv
import pytest
from efesto.identifiers import ProteinIdentifiers, scoped_id, export_ids, write_id_map
from efesto.writers import write_anvio_functions, write_anvio_misc_data, write_hit_faa, write_hit_fna, write_gff3
from efesto._dev_figures import attach_coordinates


def rows():
    return [dict(genome=g,orf='contig_1',contig='contig',bakta_id='BAKTA_1',
                 sequence=s,gene_name='MtrA',hmm_stem='MtrA',cat='Iron-Reduction',
                 cluster_id=i,bitscore=200,evalue=1e-40,cutoff=140,heme_motifs=0)
            for i,(g,s) in enumerate([('a.faa','MAST'),('b.faa','MKKL')])]


def test_scoped_queries_preserve_sequences_and_tool_assignments(tmp_path):
    rr=rows();ids=ProteinIdentifiers({r['genome']:{r['orf']:r['sequence']} for r in rr})
    path=tmp_path/'flagged.faa';assert ids.write_fasta(path,rr+rr)==2
    assert path.read_text()=='>a.faa::contig_1\nMAST\n>b.faa::contig_1\nMKKL\n'
    grouped=ids.group_hits({scoped_id('a.faa','contig_1'):'first',scoped_id('b.faa','contig_1'):'second'},rr)
    assert grouped=={'a.faa':{'contig_1':'first'},'b.faa':{'contig_1':'second'}}


def test_external_annotations_reject_ambiguous_ids_even_without_second_hit():
    rr=rows();ids=ProteinIdentifiers({r['genome']:{r['orf']:r['sequence']} for r in rr})
    with pytest.raises(ValueError,match='Ambiguous annotation query'):
        ids.group_hits({'contig_1':'KO'},rr[:1],allow_raw=True)
    assert ids.group_hits({'irrelevant':'KO'},rr,allow_raw=True)=={}
    assert ProteinIdentifiers({'a.faa':{'contig_1':'MAST'}}).group_hits({'contig_1':'KO'},rr[:1],True)=={'a.faa':{'contig_1':'KO'}}


def test_scoped_encoding_is_unambiguous():
    assert scoped_id('a::b','c')!=scoped_id('a','b::c')
    ids=ProteinIdentifiers({'a::b':{'c %':'MA'}})
    row={'genome':'a::b','orf':'c %'}
    assert ids.group_hits({scoped_id('a::b','c %'):1},[row])=={'a::b':{'c %':1}}


def test_exports_retain_both_genomes_and_crosswalk(tmp_path):
    rr=rows()
    for writer,name in [(write_anvio_functions,'functions'),(write_anvio_misc_data,'scores')]:
        path=tmp_path/name;writer(path,rr+rr)
        data=list(csv.DictReader(path.open(),delimiter='\t'))
        assert [r['gene_callers_id'] for r in data]==['a.faa::contig_1','b.faa::contig_1']
    path=tmp_path/'hits.faa';write_hit_faa(path,rr)
    assert [l.split()[0] for l in path.read_text().splitlines() if l.startswith('>')]==['>a.faa::contig_1','>b.faa::contig_1']
    coords={r['genome']:{r['orf']:dict(contig='contig',start=1+i*3,end=6+i*3,strand='+')} for i,r in enumerate(rr)}
    path=tmp_path/'Efesto-hits.gff3';write_gff3(path,rr,coords)
    lines=[l.split('\t') for l in path.read_text().splitlines() if not l.startswith('#')]
    assert [l[0] for l in lines]==['a.faa::contig','b.faa::contig']
    assert 'ID=a.faa::contig_1' in lines[0][8]
    # Long-format IDs remain original; coordinates must join by genome as well.
    report=[{'genome':r['genome'],'orf':r['orf'],'contig':'contig'} for r in rr]
    _,warnings=attach_coordinates(report,tmp_path)
    assert not warnings
    assert [r['_coordinates'] for r in report]==[(1,6,'+'),(4,9,'+')]
    mapped=[dict(r,orf='BAKTA_1') for r in report]
    _,warnings=attach_coordinates(mapped,tmp_path)
    assert not warnings
    assert [r['_coordinates'] for r in mapped]==[(1,6,'+'),(4,9,'+')]
    for g in ('a','b'):(tmp_path/(g+'.fna')).write_text('>contig\nATGCGTAAATTT\n')
    path=tmp_path/'hits.fna';assert write_hit_fna(path,rr,coords,tmp_path)==2
    assert path.read_text().count('::contig_1')==2
    path=tmp_path/'map';write_id_map(path,rr,rr)
    assert len(list(csv.DictReader(path.open(),delimiter='\t')))==2


def test_unique_export_labels_stay_compatible():
    rr=rows()[:1];assert export_ids(rr)=={('a.faa','contig_1'):'BAKTA_1'}
    # Different ORFs mapping to one locus also need distinct IDs.
    second=dict(rr[0],orf='other')
    assert len(set(export_ids(rr+[second]).values()))==2
