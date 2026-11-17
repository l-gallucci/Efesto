# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Efesto contributors.

"""Actual rule decisions, final rescue status and genome-scoped candidate identity."""

import csv
import json
import shutil
import pytest

from efesto.decisions import FilterDecisionLog
from efesto.operon import (filter_cluster_fegenie, filter_cluster_json,
                           second_pass, _REPORT_ALL_PATTERNS)


def hit(stem, orf=None, category="Iron-Reduction", genome="g.faa"):
    return dict(cat=category, genome=genome, contig="ctg", orf=orf or stem,
                hmm_stem=stem, bitscore=109.6, evalue=1e-30, cutoff=40,
                search_bitscore_cutoff=40, threshold_type="model_specific",
                threshold_source="model_cutoff_table", model_source="test",
                cluster_id=0, contig_len=10000, retention_mode="context_filter")


def audit(tmp_path, rows, context, retained, events):
    log = FilterDecisionLog(directory=tmp_path)
    try:
        log.record(rows, events, context, retained, engine="test", version="hash")
        path=tmp_path/'decisions.tsv'
        log.write(path, retained, {}, {})
        with path.open() as h:return list(csv.DictReader(h,delimiter='\t'))
    finally:
        log.close()


def test_mtrd_like_candidate_is_traced_as_dfe_exclusion_without_removing_mtr_members(tmp_path):
    rows=[hit('MtrA'),hit('MtrB_TIGR03509'),hit('MtrC_TIGR03507'),hit('DFE_0461','MtrD_locus')]
    events=[]
    filtered=filter_cluster_fegenie(rows,_REPORT_ALL_PATTERNS,decisions=events)
    assert filtered==filter_cluster_fegenie(rows,_REPORT_ALL_PATTERNS)
    records={r['original_orf']:r for r in audit(tmp_path,rows,filtered,filtered,events)}
    mtrd=records['MtrD_locus']
    assert mtrd['final_decision']=='excluded' and mtrd['final_cluster_id']==''
    assert json.loads(mtrd['reason_codes'])==['DFE2:minimum_models_not_met']
    decision=json.loads(mtrd['rule_decisions'])[0]
    assert decision['observed_models']==['DFE_0461'] and decision['minimum_models']==3
    assert records['MtrA']['final_decision']=='retained'
    assert records['MtrA']['category_before_filter']==records['MtrA']['final_category']=='Iron-Reduction'


@pytest.mark.parametrize('on_fail,expected',[('drop',set()),('passthrough_non_members',{'other'}),('keep_all',{'DFE_0461','other'})])
def test_json_failure_records_actual_cluster_or_member_scope(tmp_path,on_fail,expected):
    rows=[hit('DFE_0461'),hit('other')]
    rule={'name':'custom_dfe','genes':['DFE_0461','DFE_0462'],'rule':'require_n_of','min_genes':2,'on_fail':on_fail}
    events=[]
    filtered=filter_cluster_json(rows,[rule],[],decisions=events)
    assert {r['orf'] for r in filtered}==expected
    assert filtered==filter_cluster_json(rows,[rule],[])
    records=audit(tmp_path,rows,filtered,filtered,events)
    for row in records:
        if row['final_decision']=='excluded':
            assert json.loads(row['reason_codes'])==['custom_dfe:minimum_models_not_met']
    if on_fail=='keep_all':
        assert events[0]['status']=='failed' and events[0]['excluded_orfs']==[]


@pytest.mark.parametrize('on_fail',['drop','passthrough_non_members'])
def test_json_category_drop_explains_affected_nonmembers(tmp_path,on_fail):
    rows=[hit('transport',category='Siderophore-Transport'),hit('unrelated')]
    rules=[{'name':'transport_minimum','categories':['Siderophore-Transport'],'rule':'require_n_cat','min_genes':2,'on_fail':on_fail}]
    events=[]
    assert filter_cluster_json(rows,rules,[],decisions=events)==[]
    records=audit(tmp_path,rows,[],[],events)
    assert all(json.loads(r['reason_codes'])==['transport_minimum:minimum_models_not_met'] for r in records)


def test_json_relaxed_threshold_and_existing_two_model_exception_are_explicit():
    rows=[hit('a'),hit('b')]
    rule={'name':'minimum','genes':['a','b','c','d'],'rule':'require_n_of','min_genes':4,'on_fail':'drop'}
    events=[]
    assert filter_cluster_json(rows,[rule],[],contig_len=1000,decisions=events)==rows
    assert events[0]['minimum_models']==2 and events[0]['relaxed'] is True
    rule=dict(name='category_min',categories=['Iron-Reduction'],rule='require_n_cat_or_lone_trusted',min_genes=5,on_fail='drop')
    events=[]
    assert filter_cluster_json(rows,[rule],[],decisions=events)==rows
    assert events[0]['reason_code']=='two_model_exception'


@pytest.mark.parametrize('options,reason',[({'all_results':True},'all_results'),({'catalog_mode':True},'catalog_context_bypass')])
def test_explicit_context_bypass_is_distinguished_from_satisfied_rule(options,reason):
    rows=[hit('DFE_0461')];events=[]
    assert filter_cluster_fegenie(rows,_REPORT_ALL_PATTERNS,decisions=events,**options)==rows
    assert events[0]['status']=='bypassed' and events[0]['reason_code']==reason


def test_default_anchor_failure_and_trusted_lone_are_traced():
    rows=[hit('FoxY',category='Iron-Oxidation'),hit('LbtU-LvtA-PiuA-PirA-RhtA',category='Siderophore-Transport')]
    events=[]
    retained=filter_cluster_fegenie(rows,_REPORT_ALL_PATTERNS,decisions=events)
    assert [r['hmm_stem'] for r in retained]==['LbtU-LvtA-PiuA-PirA-RhtA']
    assert any(e['rule_id']=='FOXEYZ' and e['required_anchors']==['FoxE'] and e['status']=='failed' for e in events)
    assert any(e['reason_code']=='trusted_lone_observed' for e in events)


@pytest.mark.parametrize('sequence,kept',[('M'+'A'*364,False),('MCXXCH'+'A'*359,True),('MCXXCH',False)])
def test_second_pass_reports_actual_cyc2_sequence_gate(tmp_path,sequence,kept):
    rows=[hit('Cyc2')];events=[]
    retained=second_pass(rows,{}, {'g.faa':{'Cyc2':sequence}},decisions=events)
    assert bool(retained)==kept
    record=audit(tmp_path,rows,rows,retained,events)[0]
    assert record['context_decision']=='retained'
    assert record['sequence_decision']==('retained' if kept else 'excluded')
    assert json.loads(record['rule_decisions'])[0]['protein_length']==len(sequence)


def test_spooled_rows_preserve_identity_and_initial_decision_after_rescue(tmp_path):
    a=hit('DFE_0461','same_id',genome='a.faa');b=hit('DFE_0461','same_id',genome='b.faa')
    events=[]
    filter_cluster_fegenie([a],_REPORT_ALL_PATTERNS,decisions=events)
    log=FilterDecisionLog(directory=tmp_path)
    try:
        log.record([a],events,[],[],engine='fegenie_derived',version='hash')
        log.record([b],events,[],[],engine='fegenie_derived',version='hash')
        # Later rescue mutates the original candidate, but cannot overwrite the spool.
        a.update(cluster_id=99,retention_mode='eggnog_rescue',cat='Iron-RedoxPotential')
        path=tmp_path/'decision.tsv'
        log.write(path,[a],{}, {'a.faa':{'same_id':'BAKTA1'}})
        with path.open() as h:records={r['genome']:r for r in csv.DictReader(h,delimiter='\t')}
        assert records['a.faa']['final_decision']=='rescued'
        assert records['a.faa']['orf']=='BAKTA1'
        assert records['a.faa']['raw_cluster_id']=='0' and records['a.faa']['final_cluster_id']=='99'
        assert records['a.faa']['category_before_filter']=='Iron-Reduction'
        assert records['a.faa']['final_category']=='Iron-RedoxPotential'
        assert records['b.faa']['final_decision']=='excluded'
    finally:
        log.close()


def test_empty_decision_table_has_header(tmp_path):
    log=FilterDecisionLog(directory=tmp_path)
    try:
        path=tmp_path/'empty.tsv';log.write(path,[],{}, {})
        with path.open() as h:
            reader=csv.DictReader(h,delimiter='\t')
            assert 'final_decision' in reader.fieldnames and list(reader)==[]
    finally:
        log.close()


@pytest.mark.skipif(not shutil.which('hmmsearch'),reason='HMMER required')
def test_cli_writes_unflagged_context_rejections_without_retaining_them(tmp_path):
    from tests.test_pipeline_regressions import setup_library, run_cli
    from efesto.io import read_fasta
    from pathlib import Path
    root=Path(__file__).resolve().parents[1]
    lib=setup_library(tmp_path,'MtrA')
    (lib/'operon_rules.json').write_text(json.dumps({'rules':[{'name':'minimum2','genes':['MtrA','other'],
        'rule':'require_n_of','min_genes':2,'on_fail':'drop'}]}))
    seq=next(iter(read_fasta(root/'data/seeds/curated/MtrA_seeds.faa').values()))
    (tmp_path/'faa/a.faa').write_text('>ctg_1\n'+seq+'\n')
    (tmp_path/'gff/a.gff').write_text('##sequence-region ctg 1 10000\nctg\tcaller\tCDS\t100\t1100\t.\t+\t0\tID=ctg_1\n')
    result=run_cli(tmp_path,lib)
    assert result.returncode==0,result.stderr
    with (tmp_path/'out/Efesto-results-long.tsv').open() as h:assert list(csv.DictReader(h,delimiter='\t'))==[]
    with (tmp_path/'out/Efesto-confirmation-candidates.tsv').open() as h:assert list(csv.DictReader(h,delimiter='\t'))==[]
    with (tmp_path/'out/Efesto-filter-decisions.tsv').open() as h:records=list(csv.DictReader(h,delimiter='\t'))
    assert len(records)==1 and records[0]['final_decision']=='excluded'
    assert records[0]['coordinate_status']=='available' and records[0]['contig_length_status']=='known'
    assert json.loads(records[0]['reason_codes'])==['minimum2:minimum_models_not_met']
