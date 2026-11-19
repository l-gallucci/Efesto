# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Efesto contributors.
"""Functional labels must not change HMMs or silently bypass context rules."""
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import random
import shutil
import subprocess

import pytest

from efesto.categories import context_category
from efesto.evidence import context_evidence
from efesto.io import filter_categories, load_registry, read_cutoffs, read_map
from efesto.operon import (filter_cluster_fegenie, filter_cluster_json,
                           _DEFAULT_OPERON_RULES, _REPORT_ALL_PATTERNS)

ROOT = Path(__file__).resolve().parents[1]
LIB = ROOT / 'hmm_library'
DATA = ROOT / 'data/validation/siderophores'
BEFORE = json.loads((DATA / 'profile_integrity.json').read_text())
REG = {r['stem']: r for r in load_registry(LIB)}
CHANGED = [s for s in BEFORE if BEFORE[s]['category'] != REG[s]['category']]


def hit(stem, category, i=0):
    return dict(hmm_stem=stem, cat=category, genome='g', orf=str(i), cluster_id=0)


def test_profile_bytes_thresholds_status_and_metadata_stay_consistent():
    cutoffs = read_cutoffs(LIB / 'HMM-bitcutoffs.txt')
    names = read_map(LIB / 'FeGenie-map.txt')
    metadata = list(csv.DictReader((LIB / 'siderophore_metadata.tsv').open(), delimiter='\t'))
    assert len(metadata) == len({r['stem'] for r in metadata}) == len(BEFORE) == 236
    assert len(CHANGED) == 15
    assert set(BEFORE) == {r['stem'] for r in metadata}
    for m in metadata:
        stem = m['stem']; row = REG[stem]; old = BEFORE[stem]
        assert hashlib.sha256((LIB / row['hmm_file']).read_bytes()).hexdigest() == old['sha256']
        assert float(row['cutoff']) == float(old['cutoff']) == cutoffs.get(stem, 0.0)
        expected_status = 'experimental' if stem == 'K14759' or stem.startswith('PF06276-') else old['status']
        assert row['status'] == expected_status == m['status']
        assert row['category'] == m['category']
        assert row['gene_name'] == names[stem] == m['gene_name']
        assert m['interpretation'] and m['functional_reference']
        assert context_category(hit(stem, row['category'])) == m['context_category']


def test_irob_models_are_distinct_but_share_the_functional_target():
    iro = next(r for s, r in REG.items() if s.startswith('Sid_IroB_'))
    ko = REG['K23725']
    assert iro['category'] == ko['category'] == 'Siderophore-Synthesis'
    assert iro['nseq'] != ko['nseq']
    assert (LIB / iro['hmm_file']).read_bytes() != (LIB / ko['hmm_file']).read_bytes()
    md = {r['stem']: r for r in csv.DictReader((LIB / 'siderophore_metadata.tsv').open(), delimiter='\t')}
    assert md[iro['stem']]['functional_overlap_group'] == md['K23725']['functional_overlap_group'] != ''


@pytest.mark.parametrize('selector', ['Fe', 'Fe-acquisition', 'Fe-metabolism'])
def test_fe_selectors_include_all_reclassified_roles(selector):
    cats = {REG[s]['category']: [(s, Path(REG[s]['hmm_file']))] for s in CHANGED}
    assert set(filter_categories(cats, [selector])) == set(cats)
    for cat in cats:
        assert set(filter_categories(cats, [cat])) == {cat}
    assert not set(filter_categories(cats, ['Fe-reduction']))
    assert set(filter_categories(cats, ['metallophores'])) == {'Siderophore-Synthesis'}


def test_default_and_supplied_json_context_retention_matches_prior_grouping():
    rng = random.Random(20261119)
    other = [('synth_a','Siderophore-Synthesis'),('synth_b','Siderophore-Synthesis'),
             ('transport_a','Siderophore-Transport'),('transport_b','Siderophore-Transport'),
             ('copper','MetalResistance-Copper'),('FoxY','Iron-Oxidation')]
    pool = [(s, REG[s]['category']) for s in CHANGED] + other
    json_rules = json.loads((ROOT / 'operon_rules.json').read_text())['rules']
    engines = [lambda rows: filter_cluster_fegenie(rows, _REPORT_ALL_PATTERNS),
               lambda rows: filter_cluster_json(rows, _DEFAULT_OPERON_RULES, _REPORT_ALL_PATTERNS),
               lambda rows: filter_cluster_json(rows, json_rules, _REPORT_ALL_PATTERNS),
               lambda rows: filter_cluster_json(rows, json_rules, _REPORT_ALL_PATTERNS, contig_len=1200)]
    cases = [[(s, REG[s]['category'])] for s in CHANGED]
    cases += [rng.sample(pool, rng.randint(1, len(pool))) for _ in range(200)]
    for case in cases:
        new = [hit(s, c, i) for i, (s, c) in enumerate(case)]
        old = [dict(r, cat=BEFORE[r['hmm_stem']]['category']) if r['hmm_stem'] in BEFORE else dict(r) for r in new]
        for engine in engines:
            assert [r['orf'] for r in engine(new)] == [r['orf'] for r in engine(old)]
    # Neither a new label nor an unresolved domain can bypass its former minimum.
    for s in CHANGED:
        assert filter_cluster_fegenie([hit(s, REG[s]['category'])], _REPORT_ALL_PATTERNS) == []


def test_custom_json_rule_can_address_the_reported_function_and_evidence_agrees():
    stem = next(s for s in CHANGED if s.startswith('Sid_BesA_'))
    rows = [hit(stem, REG[stem]['category'])]
    custom = [dict(name='utilization', categories=['Siderophore-Utilization'],
                   rule='require_n_cat', min_genes=2, on_fail='drop')]
    assert filter_cluster_json(rows, custom, _REPORT_ALL_PATTERNS) == []
    records = context_evidence(rows, _DEFAULT_OPERON_RULES, _REPORT_ALL_PATTERNS,
                               mode='normal', coords={'0':dict(contig='c',start=10,end=200)})
    assert len(records) == 1
    assert records[0]['rule_id'] == 'SIDERO_SYNTH'
    assert records[0]['category_basis'] == 'context'
    assert records[0]['rule_status'] == 'not_satisfied'
    assert context_category(hit(stem, 'Custom-Function')) == 'Custom-Function'


def test_curated_labels_and_categories_survive_reimport():
    spec = importlib.util.spec_from_file_location('sidero_importer', ROOT / 'scripts/curate_hmm_library.py')
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    for stem in CHANGED:
        model = dict(stem=stem, category=BEFORE[stem]['category'], gene_name=BEFORE[stem]['gene_name'])
        result = module._apply_category_curations([model])[0]
        assert result['category'] == REG[stem]['category']
        assert result['gene_name'] == REG[stem]['gene_name']


def test_reference_domain_scope_with_real_hmmer(tmp_path):
    exe = shutil.which('hmmsearch')
    if not exe:
        pytest.skip('HMMER is needed for the real protein scope regression')
    for prefix, expected in [('PF06276-', {'P39405','A0A1A9TA96','Q93AT8'}),
                             ('Sid_YqjH_', {'Q46871'})]:
        row = next(r for s,r in REG.items() if s.startswith(prefix))
        output = tmp_path / (prefix + '.tbl')
        subprocess.run([exe,'--cpu','1','--noali','--tblout',str(output),str(LIB / row['hmm_file']),str(DATA / 'references.faa')],
                       check=True,stdout=subprocess.DEVNULL)
        matched={line.split()[0].split('|')[1] for line in output.read_text().splitlines()
                 if not line.startswith('#') and float(line.split()[5]) >= float(row['cutoff'])}
        assert expected <= matched
    assert next(r for s,r in REG.items() if s.startswith('PF06276-'))['category'] == 'Siderophore-FunctionUnresolved'


def test_no_identical_active_profile_matrices_in_this_sector():
    # Ignore names, dates, consensus residues and alignment coordinates.
    # Overlapping biological targets are a separate issue from identical matrices.
    groups = {}
    for stem, row in REG.items():
        if row['status'] not in {'active', 'experimental'}:
            continue
        text = (LIB / row['hmm_file']).read_text().split('\nHMM ', 1)[1]
        matrix = []
        for line in text.splitlines()[2:]:
            cells = line.split()
            if not cells or cells[0] == '//':
                continue
            if cells[0] == 'COMPO' or cells[0].isdigit():
                cells = cells[1:21]
            else:
                cells = cells[:20]
            matrix.append(tuple(None if x == '*' else float(x) for x in cells))
        signature = tuple(matrix)
        previous = groups.setdefault(signature, stem)
        if stem in BEFORE or previous in BEFORE:
            assert previous == stem, (previous, stem)


def test_unresolved_function_models_are_experimental_and_still_enabled():
    unresolved = [r for r in REG.values() if r['stem'] == 'K14759' or r['stem'].startswith('PF06276-')]
    assert len(unresolved) == 2
    assert all(r['status'] == 'experimental' for r in unresolved)
    assert REG['K02014']['status'] == 'active'


def test_output_reclassification_table_matches_current_scope_and_previous_ids():
    with (LIB / 'siderophore_reclassification.tsv').open() as handle:
        mapping = list(csv.DictReader(handle, delimiter='\t'))
    assert len(mapping) == len({r['stem'] for r in mapping}) == len(BEFORE)
    scope = {r['stem']: r for r in csv.DictReader((LIB / 'siderophore_metadata.tsv').open(), delimiter='\t')}
    for row in mapping:
        stem = row['stem']
        assert row['category_previous'] == BEFORE[stem]['category']
        assert row['category_current'] == REG[stem]['category']
        assert row['gene_current'] == REG[stem]['gene_name']
        assert row['status_current'] == REG[stem]['status']
        assert row['interpretation'] == scope[stem]['interpretation']
    yqjh = next(r for r in mapping if r['stem'].startswith('Sid_YqjH_'))
    assert 'including siderophore utilization' in yqjh['interpretation']
