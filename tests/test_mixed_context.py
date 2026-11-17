# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Efesto contributors.
"""Unrelated systems must neither bypass nor satisfy each other's rules."""
from itertools import permutations
import pytest
from efesto.operon import filter_cluster_fegenie, _REPORT_ALL_PATTERNS


def hit(stem, category='Iron-Reduction'):
    return {'hmm_stem':stem,'cat':category,'orf':stem,'genome':'g'}


def retained(rows, **options):
    return {r['hmm_stem'] for r in filter_cluster_fegenie(rows,_REPORT_ALL_PATTERNS,**options)}


def test_two_incomplete_systems_do_not_rescue_each_other():
    rows=[hit('EetA'),hit('sufA_TIGR01997.1','IronSulfur-Assembly')]
    for ordering in permutations(rows): assert retained(list(ordering))==set()


@pytest.mark.parametrize('other',['EetA','MamA','FoxA','FoxY','DFE_0448','DFE_0461','PF04964-Flp_Fap_pilin'])
def test_complete_suf_survives_while_an_unrelated_partial_system_is_removed(other):
    suf=['sufB_TIGR01980.1','sufC_TIGR01978.1','sufS_TIGR01979.1']
    rows=[hit(other)]+[hit(s,'IronSulfur-Assembly') for s in suf]
    assert retained(rows)==set(suf)
    assert retained(rows[::-1])==set(suf)


def test_complete_fleet_does_not_bypass_suf_or_siderophore_rules():
    fleet=['EetA','EetB','Ndh2','FmnA','FmnB']
    rows=[hit(s) for s in fleet]+[hit('sufA_TIGR01997.1','IronSulfur-Assembly'),
         hit('synthesis','Siderophore-Synthesis'),hit('copper','MetalResistance-Copper')]
    assert retained(rows)==set(fleet)|{'copper'}


def test_redox_assignment_does_not_bypass_unrelated_rules_or_mutate_input():
    rows=[hit('MtoA','Iron-Oxidation'),hit('MtrB_TIGR03509'),
          hit('receptor','Siderophore-Transport')]
    result=filter_cluster_fegenie(rows,_REPORT_ALL_PATTERNS)
    assert {r['hmm_stem'] for r in result}=={'MtoA','MtrB_TIGR03509'}
    assert {r['cat'] for r in result}=={'Iron-Oxidation'}
    assert rows[1]['cat']=='Iron-Reduction'


def test_category_groups_cannot_supply_each_others_minimum():
    rows=[hit('transport','Siderophore-Transport'),hit('synthesis','Siderophore-Synthesis')]
    assert retained(rows)==set()
    trusted=hit('LbtU-LvtA-PiuA-PirA-RhtA','Siderophore-Transport')
    assert retained([trusted])=={trusted['hmm_stem']}
    assert retained([trusted,rows[1]])=={trusted['hmm_stem']}


@pytest.mark.parametrize('mode',['catalog_mode','all_results'])
def test_explicit_bypasses_still_retain_partial_families(mode):
    rows=[hit('EetA'),hit('sufA_TIGR01997.1','IronSulfur-Assembly')]
    assert retained(rows,**{mode:True})=={r['hmm_stem'] for r in rows}


@pytest.mark.parametrize('stem,category',[('OmcF','Iron-Reduction'),
                                        ('Ni_protein_hupE_2','MetalAcquisition-CobaltNickel')])
def test_unsupported_synthesis_member_does_not_drop_an_unrelated_hit(stem,category):
    rows=[hit(stem,category),hit('PvsE-AcsE_vibrioferrin_biosynthesis-rep','Siderophore-Synthesis')]
    assert retained(rows)=={stem}


def test_partial_dfe_families_do_not_preserve_one_another():
    assert retained([hit('DFE_0451'),hit('DFE_0462')])==set()
