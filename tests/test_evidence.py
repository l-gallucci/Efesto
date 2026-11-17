# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Efesto contributors.

"""Regressions for separate evidence, missing data, and mixed neighborhoods."""
import csv
import json

import pytest

from efesto.evidence import context_evidence, hit_evidence, coordinate_evidence
from efesto.operon import _DEFAULT_OPERON_RULES, _REPORT_ALL_PATTERNS
from efesto.writers import write_context_evidence, write_uniop_pair_evidence, write_long_format


def row(stem="MtrA", orf="a", category="Iron-Reduction", genome="g.faa"):
    return dict(hmm_stem=stem, orf=orf, cat=category, genome=genome,
                cluster_id=0, heme_motifs=10)


def context(rows, **kwargs):
    coords={r['orf']:dict(start=5000+i*1000,end=5500+i*1000,contig='ctg') for i,r in enumerate(rows)}
    return context_evidence(rows,_DEFAULT_OPERON_RULES,_REPORT_ALL_PATTERNS,
                            mode=kwargs.pop('mode','normal'),coords=kwargs.pop('coords',coords),**kwargs)


@pytest.mark.parametrize('stems',[['MtrA']*3,['MtrA','lanM_NF053806','lanM_NF053806']])
def test_copies_and_unrelated_genes_do_not_complete_mtr(stems):
    records=context([row(s,str(i)) for i,s in enumerate(stems)])
    mtr=next(r for r in records if r['rule_id']=='MtrMto')
    assert mtr['observed_rule_models']==['MtrA']
    assert mtr['distinct_rule_models']==1
    assert mtr['rule_status']=='category_disambiguation_only'
    assert 'MtrB_TIGR03509' in mtr['rule_models_not_observed']
    assert not any('completeness' in key for key in mtr)


def test_suf_uses_model_ids_not_translated_gene_labels():
    stems=['sufB_TIGR01980.1','sufC_TIGR01978.1','sufS_TIGR01979.1']
    rows=[dict(row(s,str(i),'IronSulfur-Assembly'),gene_name=s.split('_')[0]) for i,s in enumerate(stems)]
    suf=next(r for r in context(rows) if r['rule_id']=='SUF_OPERON')
    assert suf['distinct_rule_models']==3
    assert suf['minimum_models']==3
    assert suf['rule_status']=='satisfied_minimum'
    assert len(suf['rule_models_not_observed'])==3


def test_mixed_cluster_reports_unsatisfied_suf_even_if_retained_by_other_filter():
    rows=[row('EetA','a'),row('sufA_TIGR01997.1','b','IronSulfur-Assembly')]
    records=context(rows)
    assert {r['rule_id']:r['rule_status'] for r in records}=={'FLEET':'not_satisfied','SUF_OPERON':'not_satisfied'}
    assert all(r['assessment_scope']=='retained_hits_not_filter_trace' for r in records)


@pytest.mark.parametrize('mode',['all_results','catalog'])
def test_bypassed_rule_is_not_called_satisfied(mode):
    assert context([row('FoxA')],mode=mode)[0]['rule_status']=='bypassed'


def test_lanm_singleton_has_no_operon_requirement():
    assert context([row('lanM_NF053806',category='MetalHomeostasis-Lanthanides')])[0]['rule_status']=='not_applicable'


def test_unknown_coordinates_not_reported_as_failed_gene_count_rule():
    assert context([row('FoxA')],coords={})[0]['rule_status']=='not_evaluable_without_coordinates'


def test_unsearched_required_models_are_distinct_from_not_observed():
    record=context([row('FoxA')],available_models={'FoxA'})[0]
    assert record['rule_status']=='not_evaluable_models_not_searched'
    assert record['rule_models_not_searched']==['FoxB','FoxC']


def test_contig_distances_are_measurements_without_penalties():
    coords={'a':dict(contig='ctg',start=101,end=400)}
    assert coordinate_evidence('a',coords,{'ctg':1000})['distance_to_contig_edge_bp']==100
    assert coordinate_evidence('a',coords,{})['distance_to_contig_edge_bp']==''
    assert coordinate_evidence('b',coords,{})['coordinate_status']=='unavailable'


def evidence(**updates):
    data=dict(coords={},contig_lengths={},confirmation_map={},topology_map={},
              pending_keys=set(),eggnog_hits={},eggnog_status='not_requested',
              foldseek_hits={},foldseek_status='not_requested',topology_requested=False,
              phobius_hits={},phobius_status='not_requested',bgc_regions=[],bgc_status='not_requested',
              uniop_status='not_requested',operon_map={})
    data.update(updates)
    return hit_evidence(row(),**data)


def test_missing_optional_tools_are_not_support():
    result=evidence()
    assert result['uniop_status']=='not_requested'
    assert result['eggnog_status']=='not_applicable'
    assert result['bgc_overlap']=='not_evaluated'
    assert not any(key.endswith('_weight') for key in result)


def test_failed_phobius_does_not_hide_independent_motif_observation():
    result=evidence(topology_requested=True,phobius_status='failed',
                    topology_map={'MtrA':dict(expected_heme_c=10,expected_topology='signal_peptide')})
    assert result['phobius_status']=='failed'
    assert result['topology_verdict']=='not_evaluated'
    assert result['heme_verdict']=='confirmed'


def test_failed_structure_tool_and_empty_result_are_distinct():
    params=dict(confirmation_map={'MtrA':dict(gene_name='MtrA',kegg_ko=[])},pending_keys={('g.faa','a')})
    failed=evidence(**params,foldseek_status='failed')
    empty=evidence(**params,foldseek_status='completed')
    assert failed['structure_status']=='failed'
    assert failed['structure_label_match']=='not_evaluated'
    assert empty['structure_label_match']=='no_informative_result'


def test_structure_label_is_not_a_functional_confirmation():
    result=evidence(confirmation_map={'MtrA':dict(gene_name='MtrA',kegg_ko=[])},
                    pending_keys={('g.faa','a')},foldseek_status='completed',
                    foldseek_hits={'a':dict(target='MtrA_reference',bits=.1,evalue=10)})
    assert result['structure_label_match']=='matches_label'
    assert result['structure_evalue']==10
    assert 'structure_verdict' not in result


def test_bgc_overlap_applies_only_to_overlapping_orf():
    coords={'a':dict(contig='ctg',start=5000,end=5300)}
    params=dict(coords=coords,bgc_status='parsed',bgc_regions=[('ctg',1,1000,'siderophore')])
    result=evidence(**params)
    assert result['bgc_overlap']=='no_overlap'
    params['bgc_regions']=[('ctg',5000,7000,'siderophore')]
    result=evidence(**params)
    assert result['bgc_overlap']=='overlap'
    assert json.loads(result['bgc_overlap_regions'])[0]['product']=='siderophore'


def test_context_json_roundtrip_and_rule_hash(tmp_path):
    records=context([row(),row(orf='b')])
    out=tmp_path/'context.tsv';write_context_evidence(out,records)
    parsed=list(csv.DictReader(out.open(),delimiter='\t'))
    assert json.loads(parsed[0]['model_copy_counts'])=={'MtrA':2}
    assert len(parsed[0]['rule_version'])==64


def test_zero_uniop_probability_is_exported_for_the_correct_genome(tmp_path):
    rows=[row(orf=orf,genome=genome) for genome in ['one.faa','two.faa'] for orf in ['a','b']]
    out=tmp_path/'pairs.tsv'
    write_uniop_pair_evidence(out,rows,{'one.faa':{('a','b'):0.0},'two.faa':{('a','b'):.9}})
    pairs=list(csv.DictReader(out.open(),delimiter='\t'))
    assert {r['genome']:float(r['operon_probability']) for r in pairs}=={'one.faa':0.0,'two.faa':.9}
    assert pairs[0]['orf_a']!=pairs[1]['orf_a']


def test_sufficient_observation_does_not_require_searching_every_alternative():
    record=context([row('FoxA','a'),row('FoxB','b')],available_models={'FoxA','FoxB'})[0]
    assert record['rule_status']=='satisfied_minimum'
    assert record['rule_models_not_searched']==['FoxC']


def test_catalog_still_reports_mtr_category_disambiguation():
    assert context([row()],mode='catalog')[0]['rule_status']=='category_disambiguation_only'


def test_undefined_topology_is_not_a_neutral_topology_confirmation():
    result=evidence(topology_requested=True,phobius_status='completed',
                    topology_map={'MtrA':dict(expected_heme_c=10,expected_topology=None)},
                    phobius_hits={'a':dict(tm_helices=0,signal_peptide=True)})
    assert result['topology_verdict']=='not_applicable'
    assert result['heme_verdict']=='confirmed'


def test_failed_uniop_does_not_claim_the_gene_is_outside_an_operon():
    assert evidence(uniop_status='failed')['uniop_context']=='not_evaluated'


def test_last_cds_coordinate_is_not_a_known_contig_edge(tmp_path):
    from efesto.io import build_contig_length_dict
    faa=tmp_path/'g.faa';faa.write_text('>a\nMAAA\n')
    gff=tmp_path/'g.gff';gff.write_text('ctg\tcaller\tCDS\t100\t400\t.\t+\t0\tID=a\n')
    lengths,exact=build_contig_length_dict([faa],gff_dir=tmp_path,return_exact=True)
    assert lengths=={'g.faa':{'ctg':400}}
    assert exact=={'g.faa':{}}
    assert coordinate_evidence('a',{'a':dict(contig='ctg',start=100,end=400)},exact['g.faa'])['distance_to_contig_edge_bp']==''
    gff.write_text('##sequence-region ctg 1 1000\n'+gff.read_text())
    _,exact=build_contig_length_dict([faa],gff_dir=tmp_path,return_exact=True)
    assert exact=={'g.faa':{'ctg':1000}}
