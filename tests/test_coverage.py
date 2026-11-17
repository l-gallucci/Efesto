# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Efesto contributors.

"""Depth formats and streaming subprocess errors must not silently alter depth."""
import io
import subprocess
from pathlib import Path
import pytest
from efesto import coverage


@pytest.mark.parametrize('text,expected', [
    ('c\t1\t10\nc\t2\t20\n', {'depth':15.,'length':2,'reads':0}),
    ('#CHROM\tPOS\tfile.bam\nc\t1\t0\nc\t2\t20\n', {'depth':10.,'length':2,'reads':0}),
    ('c\t15\n', {'depth':15.,'length':0,'reads':0}),
    ('#contig\tdepth\tlength\nc\t15\t200\n', {'depth':15.,'length':200,'reads':0}),
    ('contigName\tcontigLen\ttotalAvgDepth\tx.bam\tx.bam-var\nc\t200\t15\t15\t0\n', {'depth':15.,'length':200,'reads':0}),
    ('#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\nc\t1\t200\t4\t100\t50\t15\t30\t60\n', {'depth':15.,'length':200,'reads':4}),
    ('c\t1\t200\t4\t100\t50\t15\t30\t60\n', {'depth':15.,'length':200,'reads':4}),
    ('#ID\tAvg_fold\tLength\tRef_GC\nc\t15\t200\t0.4\n', {'depth':15.,'length':200,'reads':0}),
])
def test_supported_formats(tmp_path, text, expected):
    path=tmp_path/'depth'; path.write_text(text)
    assert coverage.load_depth_file(path)=={'c':expected}


@pytest.mark.parametrize('text', ['c\t1\t10\nc\t1\t20\n', 'c\t0\t10\n',
                                  'c\t1\tNaN\n', 'c\t1\t-1\n', 'c\t1\t10\t20\n',
                                  'c\t15\nc\t20\n', 'c\t1\tbad\n'])
def test_invalid_depth_rejected_with_location(tmp_path,text):
    path=tmp_path/'depth';path.write_text(text)
    with pytest.raises(ValueError,match=r'depth:\d+:'):
        coverage.load_depth_file(path)


def test_sparse_positions_mean_only_reported_bases(tmp_path):
    path=tmp_path/'depth';path.write_text('c\t2\t10\nc\t4\t20\n')
    assert coverage.load_depth_file(path)['c']=={'depth':15., 'length':4, 'reads':0}


def test_missing_file_is_an_error(tmp_path):
    with pytest.raises(FileNotFoundError):coverage.load_depth_file(tmp_path/'missing')


def test_depth_parser_consumes_iterator():
    class Stream:
        def __iter__(self):return iter(['c\t1\t10\n', 'c\t2\t20\n'])
        def read(self,*args):pytest.fail('Must stream, not read whole file')
        def readlines(self):pytest.fail('Must stream, not read all lines')
    assert coverage._parse_depth_lines(Stream(), 'stream')['c']['depth']==15


@pytest.mark.parametrize('code',[0,2])
def test_bam_fallback_streams_and_checks_status(tmp_path,monkeypatch,code):
    monkeypatch.setattr(coverage.subprocess, 'run',lambda *a,**k:subprocess.CompletedProcess(a,1,'','coverage unavailable'))
    calls=[]
    class Process:
        stdout=io.StringIO('c\t1\t10\nc\t2\t20\n')
        def wait(self):return code
        def kill(self):pass
    def launch(command,**kwargs):
        calls.append((command,kwargs));return Process()
    monkeypatch.setattr(coverage.subprocess,'Popen',launch)
    if code:
        with pytest.raises(RuntimeError,match='exit 2'):coverage.compute_coverage_from_bam('test.bam',tmp_path,'test')
    else:
        assert coverage.compute_coverage_from_bam('test.bam',tmp_path,'test')['c']['depth']==15
    assert calls[0][0]==['samtools','depth','-aa','test.bam']
    assert calls[0][1]['stdout']==subprocess.PIPE
    assert calls[0][1]['stderr']!=subprocess.PIPE



def test_equal_depth_is_equal_abundance_despite_contig_length():
    info={'short':{'depth':10.,'length':1000}, 'long':{'depth':10.,'length':2000}}
    result=coverage.compute_coverage_per_million(info)
    assert result=={'short':500000.,'long':500000.}
    assert coverage.compute_tpm(info,{'short':500,'long':9000})==result


def test_normalisation_zero_and_empty_inputs():
    assert coverage.compute_coverage_per_million({})=={}
    assert coverage.compute_coverage_per_million({'a':{'depth':0.}})=={'a':0.}
    assert sum(coverage.compute_coverage_per_million({'a':{'depth':10},'b':{'depth':20}}).values())==pytest.approx(1e6)


def test_heatmap_normalises_without_contig_lengths_once_per_genome(tmp_path,monkeypatch):
    from efesto import writers
    calls=[];original=writers.compute_coverage_per_million
    def count(info):calls.append(info);return original(info)
    monkeypatch.setattr(writers,'compute_coverage_per_million',count)
    rows=[{'genome':'g.faa','contig':'a','cat':'Fe'}, {'genome':'g.faa','contig':'b','cat':'Fe'}]
    path=tmp_path/'out.csv'
    writers.write_coverage_heatmap(path,rows,['g.faa'],{'g.faa':{'a':{'depth':10},'b':{'depth':30}}},True)
    assert len(calls)==1
    assert path.read_text()=='# coverage metric: coverage_per_million_sum\nX,g.faa\nFe,1000000.0000\n'
