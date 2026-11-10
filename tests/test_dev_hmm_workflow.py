import csv
import importlib.util
import json
from pathlib import Path

import pytest

# Works both in the analysis directory and in Efesto/scripts after application.
ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/'scripts/dev_hmm_workflow.py'
if not SCRIPT.exists():SCRIPT=ROOT/'dev_hmm_workflow.py'
spec=importlib.util.spec_from_file_location('dev_hmm_workflow_under_test',SCRIPT)
W=importlib.util.module_from_spec(spec);spec.loader.exec_module(W)


def row(i,label='target',group=None,partition='development'):
    return dict(id=i,label=label,group=group or i,partition=partition,evidence='experimental_anchor',reference='doi:test')


def inputs(tmp_path, rows, seqs):
    faa=tmp_path/'in.faa';faa.write_text(''.join('>'+i+'\n'+s+'\n' for i,s in seqs.items()))
    manifest=tmp_path/'manifest.tsv'
    with manifest.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    return faa,manifest


def test_exact_alias_cannot_leak_across_groups(tmp_path):
    faa,m=inputs(tmp_path,[row('a'),row('b')],{'a':'ACDEF','b':'ACDEF'})
    with pytest.raises(ValueError,match='Identical sequences'):W.read_inputs(faa,m)


def test_homology_group_cannot_cross_external_boundary(tmp_path):
    faa,m=inputs(tmp_path,[row('a',group='same'),row('b',group='same',partition='external')],{'a':'ACDEF','b':'ACDEG'})
    with pytest.raises(ValueError,match='leaks'):W.read_inputs(faa,m)


def test_unknown_cannot_become_test_truth(tmp_path):
    faa,m=inputs(tmp_path,[row('a',label='unknown')],{'a':'ACDEF'})
    with pytest.raises(ValueError,match='Unknown labels'):W.read_inputs(faa,m)


def test_conflicting_alias_labels_are_rejected(tmp_path):
    faa,m=inputs(tmp_path,[row('a',group='same'),row('b',label='background',group='same')],{'a':'ACDEF','b':'ACDEF'})
    with pytest.raises(ValueError,match='Conflicting'):W.read_inputs(faa,m)


def test_aliases_do_not_inflate_support(tmp_path):
    faa,m=inputs(tmp_path,[row('a',group='same'),row('b',group='same')],{'a':'ACDEF','b':'ACDEF'})
    _,rows=W.read_inputs(faa,m)
    assert len(rows)==1


def test_folds_keep_all_members_of_a_group_together():
    rows=[row(str(i),group='g'+str(i//2)) for i in range(12)]
    rows += [row('n'+str(i),'background','n'+str(i//2)) for i in range(10)]
    folds=W.make_folds(rows,3)
    groups=[{r['group'] for r in f} for f in folds]
    assert all(not groups[a]&groups[b] for a in range(3) for b in range(a))
    assert all(any(r['label']=='target' for r in f) for f in folds)


def test_sparse_groups_are_not_padded_with_duplicate_sequences():
    with pytest.raises(ValueError,match='Insufficient'):W.make_folds([row('a',group='one'),row('b',group='one')],2)


def test_cutoff_tie_preserves_specificity_and_runtime_evalue():
    rows=[row('p'),row('n','background')]
    scores={'p':dict(score=30,evalue=1e-9),'n':dict(score=40,evalue=1)}
    cutoff,m=W.choose_cutoff(rows,scores,'target')
    assert cutoff==30 and m['tp']==1 and m['fp']==0
    assert not W.accepts(None,-100,1e-5)


def test_abstention_can_exclude_the_highest_scoring_negative():
    rows=[row('p'),row('n','background')]
    scores={'n':dict(score=100,evalue=1e-30)}
    cutoff,m=W.choose_cutoff(rows,scores,'target')
    assert cutoff>100 and m['fp']==0 and m['fn']==1


def test_competition_uses_predeclared_margin_and_preserves_raw_comparison():
    rows=[row('q')];scores={'A':{'q':dict(score=200,evalue=1e-30)},'B':{'q':dict(score=150,evalue=1e-20)}}
    assert W.predictions(rows,scores,{'A':190,'B':100},1e-5)[0]['predicted']=='B'
    assert W.predictions(rows,scores,{'A':190,'B':100},1e-5,'raw')[0]['predicted']=='A'


def test_nested_training_never_sees_outer_external_or_unknown_records(tmp_path,monkeypatch):
    rows=[row('p'+str(i)) for i in range(3)]+[row('n'+str(i),'background') for i in range(6)]
    rows += [row('external',partition='external'),row('candidate','unknown',partition='unlabelled')]
    aa='ACDEFGHIKLMNPQRSTVWY';seqs={r['id']:'ACDEF'+aa[i] for i,r in enumerate(rows)}
    faa,m=inputs(tmp_path,rows,seqs);calls=[]
    def build(train,evaluate,seqs,targets,out,tools):
        out.mkdir(parents=True)
        calls.append((out,train,evaluate))
        return {t:{r['id']:dict(score=100 if r['label']==t else 10,evalue=1e-20) for r in evaluate} for t in targets}
    monkeypatch.setattr(W,'build_and_score',build)
    result=W.workflow(faa,m,tmp_path/'out')
    assert all(not {'external','candidate'} & {r['id'] for r in train} for _,train,_ in calls)
    for path,train,evaluate in calls:
        if 'outer_' in str(path):
            assert not {r['group'] for r in train}&{r['group'] for r in evaluate}
            assert not {'external','candidate'}&{r['id'] for r in evaluate}
    assert result['external_test']['target']['tp']==1
    assert result['promotion']['eligible'] is False
    assert len(result['outer_predictions'])==9
    assert json.loads((tmp_path/'out/result.json').read_text())['final_cutoffs']['target']==100


def test_unusual_training_residue_needs_explicit_alignment_policy(tmp_path):
    faa,m=inputs(tmp_path,[row('a')],{'a':'ACDEU'})
    with pytest.raises(ValueError,match='alignment policy'):W.read_inputs(faa,m)
