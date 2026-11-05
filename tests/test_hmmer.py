"""Search scheduling must preserve hits while bounding concurrency and memory."""
import subprocess
import threading
import time
from pathlib import Path

import pytest

from efesto import hmmer


@pytest.mark.parametrize('cpus,hmmer_workers', [(1,'1'), (2,'2'), (4,'4')])
def test_hmmsearch_cpu_mode_preserves_thresholds(tmp_path,monkeypatch,cpus,hmmer_workers):
    commands=[]
    def execute(command,**kwargs):
        commands.append(command)
        return subprocess.CompletedProcess(command,0,stdout='',stderr='')
    monkeypatch.setattr(hmmer.subprocess,'run',execute)
    target=tmp_path/'result.tblout'
    assert hmmer._hmmsearch_job(('model.hmm','proteins.faa',str(target),73.5,cpus,30))==(str(target),True,'')
    command=commands[0]
    assert command[command.index('--cpu')+1]==hmmer_workers
    assert command[command.index('-T')+1]=='73.5'
    assert command[-2:]==['model.hmm','proteins.faa']


def test_existing_result_does_not_relaunch_hmmer(tmp_path,monkeypatch):
    target=tmp_path/'result.tblout';target.write_text('# [ok]\n')
    def unexpected(*args,**kwargs):
        pytest.fail('Cached HMMER output should not be recomputed')
    monkeypatch.setattr(hmmer.subprocess,'run',unexpected)
    assert hmmer._hmmsearch_job(('model.hmm','proteins.faa',str(target),0,1,30))==(str(target),True,'')


@pytest.mark.parametrize('threads_total,expected_workers', [(1,1),(4,2),(16,2)])
def test_searches_complete_with_bounded_queue_and_cpu_budget(tmp_path,monkeypatch,capsys,threads_total,expected_workers):
    files=[Path(f'genome_{i}.faa') for i in range(3)]
    models={'iron':[(f'h{i}',Path(f'h{i}.hmm')) for i in range(7)]}
    seen=[];inflight=[];active=maximum=0;lock=threading.Lock()
    actual_wait=hmmer.wait
    def observe_wait(futures,**kwargs):
        inflight.append(len(futures))
        return actual_wait(futures,**kwargs)
    def search(job):
        nonlocal active,maximum
        with lock:
            active+=1;maximum=max(maximum,active);seen.append(job)
        time.sleep(.002)
        with lock:active-=1
        # Exercise unsuccessful tool results without losing remaining jobs.
        ok=job[0]!='h0.hmm'
        return job[2],ok,'' if ok else 'test search failure'
    monkeypatch.setattr(hmmer,'wait',observe_wait)
    monkeypatch.setattr(hmmer,'_hmmsearch_job',search)
    hmmer.run_all_hmmsearches(files,models,{'h1':62},tmp_path,threads_total=threads_total,hmm_threads=1)
    assert len(seen)==21
    assert len({(j[0],j[1]) for j in seen})==21
    assert maximum<=expected_workers
    assert max(inflight)<=2*expected_workers
    assert all(j[4:]==(1,30) for j in seen)
    assert {j[3] for j in seen if j[0]=='h1.hmm'}=={62}
    output=capsys.readouterr()
    assert '21/21  (3 errors)' in output.out
    assert output.err.count('test search failure')==3


def test_worker_exception_propagates(tmp_path,monkeypatch):
    def fail(job):raise OSError('search launcher failed')
    monkeypatch.setattr(hmmer,'_hmmsearch_job',fail)
    with pytest.raises(OSError,match='search launcher failed'):
        hmmer.run_all_hmmsearches([Path('g.faa')],{'iron':[('m',Path('m.hmm'))]}, {},tmp_path,1)


def test_hmm_worker_budget_limits_simultaneous_searches(tmp_path,monkeypatch):
    observed=[]
    actual_pool=hmmer.ThreadPoolExecutor
    def pool(*args,**kwargs):
        observed.append(kwargs['max_workers'])
        return actual_pool(*args,**kwargs)
    monkeypatch.setattr(hmmer,'ThreadPoolExecutor',pool)
    monkeypatch.setattr(hmmer,'_hmmsearch_job',lambda job:(job[2],True,''))
    hmmer.run_all_hmmsearches([Path('g.faa')],{'iron':[(str(i),Path(f'{i}.hmm')) for i in range(4)]},{},tmp_path,threads_total=4,hmm_threads=4)
    assert observed==[1]
