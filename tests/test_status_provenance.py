# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Efesto contributors.
"""Model status and provenance survive result export and library reconstruction."""
import csv
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from efesto.decisions import FilterDecisionLog
from efesto.io import load_registry
from efesto.writers import write_long_format

ROOT = Path(__file__).resolve().parents[1]
LIB = ROOT / 'hmm_library'
REG = {r['stem']:r for r in load_registry(LIB)}


def importer():
    spec = importlib.util.spec_from_file_location('provenance_importer', ROOT/'scripts/curate_hmm_library.py')
    module = importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def read_table(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def model(stem):
    r = dict(REG[stem]);r['hmm_file'] = str(LIB/r['hmm_file']);return r


def candidate(status):
    return dict(cat='Siderophore-FunctionUnresolved',genome='g.faa',contig='ctg',orf='orf1',
                gene_name='shared-domain',hmm_stem='profile',bitscore=100,evalue=1e-30,cutoff=25.2,
                threshold_type='model_specific',threshold_source='model_cutoff_table',
                search_bitscore_cutoff=25.2,model_source='fegenie',model_status=status,
                cluster_id=0,contig_len=10000,heme_motifs=0,retention_mode='context_filter')


@pytest.mark.parametrize('status',['active','experimental',None,''])
def test_long_format_exports_status_and_does_not_invent_it(tmp_path,status):
    row=candidate(status);out=tmp_path/'hits.tsv';write_long_format(out,[row])
    assert read_table(out)[0]['model_status'] == (status or 'unknown')


@pytest.mark.parametrize('retention',['retained','excluded','rescued'])
def test_status_survives_spooling_and_rescue(tmp_path,retention):
    row=candidate('experimental');log=FilterDecisionLog(directory=tmp_path)
    try:
        log.record([row],[],[] if retention!='retained' else [row],[],engine='test',version='hash')
        kept=[] if retention=='excluded' else [dict(row,retention_mode='eggnog_rescue' if retention=='rescued' else 'context_filter')]
        out=tmp_path/'decisions.tsv';log.write(out,kept,{},{});r=read_table(out)[0]
        assert r['model_status']=='experimental'
        assert r['final_decision']==retention
    finally:log.close()


def test_writer_preserves_complete_registry_metadata_and_auxiliary_tables(tmp_path):
    m=importer();pf=next(s for s in REG if s.startswith('PF06276-'))
    stems=['CbcL',pf,'IroC-salmochelin_transport-rep','rusticyanin_TIGR03095']
    models=[model(s) for s in stems]
    models[0]['future_provenance']='kept verbatim'
    out=tmp_path/'library';m.write_library(out,models,[])
    written={r['stem']:r for r in read_table(out/'hmm_registry.tsv')}
    for original in models:
        result=written[original['stem']]
        for field,value in original.items():
            if field!='hmm_file':assert result[field]==value,(field,original['stem'])
        assert (out/result['hmm_file']).read_bytes()==Path(original['hmm_file']).read_bytes()
    assert written['IroC-salmochelin_transport-rep']['status'].startswith('deprecated')
    assert written[pf]['status']=='experimental'
    scope=read_table(out/'siderophore_metadata.tsv')
    assert {r['stem'] for r in scope}=={pf,'IroC-salmochelin_transport-rep'}
    for r in scope:assert r['status']==written[r['stem']]['status']
    mapping=read_table(out/'siderophore_reclassification.tsv')
    assert {r['stem'] for r in mapping}=={r['stem'] for r in scope}
    for r in mapping:
        assert r['category_current']==written[r['stem']]['category']
        assert r['status_current']==written[r['stem']]['status']
    assert {r['stem'] for r in read_table(out/'model_licenses.tsv')}==set(stems)
    assert (out/'category_curation.tsv').exists() and (out/'category_aliases.tsv').exists()
    counts=read_table(out/'category_catalog.tsv')
    assert sum(int(r['total']) for r in counts)==4
    assert sum(int(r['experimental']) for r in counts)==1
    assert sum(int(r['deprecated']) for r in counts)==1
    # The registry-aware loader must preserve source, original date and exclusions.
    reloaded={r['stem']:r for r in m.load_fegenie(out)}
    for stem in stems:
        for key in ('source','reference','added_date','status','category_original'):
            assert reloaded[stem][key]==written[stem][key]


def test_report_route_retains_metadata_and_resolves_relative_hmm_paths(tmp_path,monkeypatch):
    m=importer();source=tmp_path/'source';source.mkdir();cat=source/'Custom';cat.mkdir()
    (cat/'custom.hmm').write_text('HMMER3/f\nNAME custom\nNSEQ 20\n//\n')
    report=source/'report.tsv'
    record=dict(stem='custom',hmm_file='Custom/custom.hmm',name='custom',acc='CUSTOM1',
                nseq='20',category='Custom',gene_name='custom',cutoff='50.5',source='custom',
                status='experimental',reference='https://example.org/source',validated_in='reference organism',
                needs_confirmation='yes',kegg_ko='K00001',expected_heme_c='2',expected_topology='periplasmic',
                added_date='2020-01-01',category_original='historical',future_provenance='keep me')
    with report.open('w',newline='') as h:
        w=csv.DictWriter(h,record.keys(),delimiter='\t');w.writeheader();w.writerow(record)
    out=tmp_path/'output';monkeypatch.setattr(sys,'argv',['curate_hmm_library.py','--from_report',str(report),'--out_dir',str(out)])
    monkeypatch.setattr(m,'deduplicate',lambda models,threshold:(models,[]))
    m.main();written=read_table(out/'hmm_registry.tsv')[0]
    for field,value in record.items():
        if field!='hmm_file':assert written[field]==value,field
    assert (out/written['hmm_file']).read_bytes()==(cat/'custom.hmm').read_bytes()
    # Known source notes must not be invented for an unrecognized profile.
    assert read_table(out/'siderophore_metadata.tsv')==[]
    assert read_table(out/'model_licenses.tsv')==[]


def test_source_notes_are_not_transferred_between_unrelated_import_channels(tmp_path):
    m=importer();original=model('K02014');original['source']='unrelated'
    out=tmp_path/'library';m.write_library(out,[original],[])
    assert read_table(out/'siderophore_metadata.tsv')==[]
    assert read_table(out/'model_licenses.tsv')==[]


@pytest.mark.skipif(not shutil.which('hmmsearch'),reason='HMMER required')
@pytest.mark.parametrize('mode',['retained','excluded','unregistered'])
def test_cli_propagates_experimental_and_unknown_status(tmp_path,mode):
    pf=next(s for s in REG if s.startswith('PF06276-'));m=importer();lib=tmp_path/'library'
    m.write_library(lib,[model(pf)],[])
    if mode=='unregistered':(lib/'hmm_registry.tsv').unlink()
    faa=tmp_path/'proteins';faa.mkdir()
    shutil.copyfile(ROOT/'data/validation/siderophores/references.faa',faa/'panel.faa')
    out=tmp_path/'out';env=os.environ.copy();env['PYTHONPATH']=str(ROOT/'src');env['PYTHONDONTWRITEBYTECODE']='1'
    cmd=[sys.executable,'-m','efesto.cli','--faa_dir',str(faa),'--hmm_dir',str(lib),
         '--out',str(out),'--threads','1']
    if mode!='excluded':cmd.append('--all_results')
    result=subprocess.run(cmd,env=env,capture_output=True,text=True)
    assert result.returncode==0,result.stdout+'\n'+result.stderr
    decisions=read_table(out/'Efesto-filter-decisions.tsv');hits=read_table(out/'Efesto-results-long.tsv')
    expected='unknown' if mode=='unregistered' else 'experimental'
    assert len(decisions)==3 and {r['model_status'] for r in decisions}=={expected}
    if mode=='excluded':
        assert hits==[] and {r['final_decision'] for r in decisions}=={'excluded'}
    else:
        assert len(hits)==3 and {r['model_status'] for r in hits}=={expected}


def test_library_loader_skips_local_staging_models(tmp_path):
    lib=tmp_path/'library';staging=lib/'_staging';staging.mkdir(parents=True)
    (staging/'unpublished.hmm').write_text('HMMER3/f\nNAME unpublished\nNSEQ 10\n//\n')
    assert importer().load_fegenie(lib)==[]
