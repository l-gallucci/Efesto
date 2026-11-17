# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Efesto contributors.
"""Prediction cache invalidation and failed-output isolation."""
import sys
from types import SimpleNamespace
import pytest
from efesto import gene_calling as gc


@pytest.fixture
def prediction(tmp_path, monkeypatch):
    class Gene:
        begin=1;strand=1;start_type='ATG';rbs_motif=None;partial_end=False
        def __init__(self, seq): self.end=len(seq);self.seq=seq
        def translate(self):
            if Finder.fail: raise ValueError('translation failed')
            return 'M'+'A'*(len(self.seq)//3)
    class Finder:
        calls=0;fail=False
        def __init__(self, **kwargs): pass
        def train(self,*sequences): pass
        def find_genes(self, seq):
            Finder.calls+=1
            return [Gene(seq)]
    monkeypatch.setitem(sys.modules,'pyrodigal',SimpleNamespace(GeneFinder=Finder))
    monkeypatch.setitem(sys.modules,'pyrodigal_gv',SimpleNamespace(ViralGeneFinder=Finder))
    monkeypatch.setattr(gc,'version',lambda _: 'test-1')
    fna=tmp_path/'sample.fna';fna.write_text('>ctg\n'+'ATG'*20+'\n')
    args=(str(fna),str(tmp_path/'sample.faa'),str(tmp_path/'sample.gff'),True,'pyrodigal')
    return Finder,fna,args


def test_cache_requires_content_configuration_and_complete_outputs(prediction,monkeypatch):
    finder,fna,args=prediction
    assert gc._pyrodigal_job(args)[1]
    assert gc._pyrodigal_job(args)[1] and finder.calls==1
    fna.write_text('>ctg\n'+'ATG'*30+'\n')
    assert gc._pyrodigal_job(args)[1] and finder.calls==2
    assert gc._pyrodigal_job((*args[:3],False,args[4]))[1] and finder.calls==3
    monkeypatch.setattr(gc,'version',lambda _: 'test-2')
    assert gc._pyrodigal_job(args)[1] and finder.calls==4
    assert gc._pyrodigal_job((*args[:4],'pyrodigal-gv'))[1] and finder.calls==5
    # Manual truncation of a completed file is also a cache miss.
    from pathlib import Path
    Path(args[1]).write_text('>partial\nM\n')
    assert gc._pyrodigal_job((*args[:4],'pyrodigal-gv'))[1] and finder.calls==6


def test_partial_failure_is_never_accepted_on_retry(prediction):
    finder,_,args=prediction
    finder.fail=True
    for _ in range(2):
        assert not gc._pyrodigal_job(args)[1]
    assert finder.calls==2
    from pathlib import Path
    assert not Path(args[1]).exists() and not Path(args[2]).exists()
    assert not Path(args[1]).with_suffix('.cache.json').exists()
    finder.fail=False
    assert gc._pyrodigal_job(args)[1] and finder.calls==3


def test_legacy_files_do_not_skip_prediction(prediction):
    from pathlib import Path
    finder,_,args=prediction
    Path(args[1]).write_text('>old\nMA\n');Path(args[2]).write_text('partial')
    assert gc._pyrodigal_job(args)[1] and finder.calls==1
    assert '>ctg_1' in Path(args[1]).read_text()


def test_failed_prediction_stops_the_run(prediction,tmp_path):
    finder,fna,_=prediction;finder.fail=True
    with pytest.raises(RuntimeError,match='annotation stopped'):
        gc.run_prodigal([fna],tmp_path/'out',meta_mode=True,gene_caller='pyrodigal')
