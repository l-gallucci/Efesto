"""Specific Mbn subfamilies must not turn isolated MNIO hits into a pathway."""
import csv
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path
import shutil
import subprocess

import pytest

from efesto.io import load_registry, read_cutoffs
from efesto.operon import (filter_cluster_fegenie, filter_cluster_json,
                           _DEFAULT_OPERON_RULES, _REPORT_ALL_PATTERNS)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data/validation/methanobactin'
GENES = ['mbnA_TIGR04071', 'mbnB_TIGR04159', 'mbnC_TIGR04160']


def hit(stem, cat='Metallophore-Synthesis'):
    return dict(hmm_stem=stem, cat=cat, genome='test', orf=stem)


def engines(rows):
    yield filter_cluster_fegenie(rows, _REPORT_ALL_PATTERNS)
    yield filter_cluster_json(rows, _DEFAULT_OPERON_RULES, _REPORT_ALL_PATTERNS)
    yield filter_cluster_json(rows, _DEFAULT_OPERON_RULES, _REPORT_ALL_PATTERNS,
                              contig_len=1200)
    rules = json.loads((ROOT / 'operon_rules.json').read_text())['rules']
    yield filter_cluster_json(rows, rules, _REPORT_ALL_PATTERNS)


@pytest.mark.parametrize('core', [list(c) for n in range(3)
                                 for c in itertools.combinations(GENES, n)])
def test_partial_mbn_is_removed_without_hiding_unrelated_evidence(core):
    other = hit('copper', 'MetalResistance-Copper')
    for result in engines([hit(s) for s in core] + [other]):
        assert result == [other]


def test_abc_families_pass_and_duplicates_do_not_rescue_missing_precursor():
    rows = [hit(s) for s in GENES]
    for result in engines(rows):
        assert result == rows
    for result in engines([hit(GENES[1]), hit(GENES[1]), hit(GENES[2])]):
        assert result == []


def test_cnt_and_mbn_require_their_own_core():
    cnt = [hit('K23120'), hit('K23122')]
    for result in engines(cnt + [hit(GENES[1])]):
        assert result == cnt
    mbn = [hit(s) for s in GENES]
    for result in engines(mbn + [hit('K23120'), hit('EetA', 'Iron-Oxidation')]):
        assert result == mbn


def test_individual_family_evidence_remains_available_without_genomic_context():
    rows = [hit(GENES[1])]
    assert filter_cluster_fegenie(rows, _REPORT_ALL_PATTERNS, catalog_mode=True) == rows
    assert filter_cluster_fegenie(rows, _REPORT_ALL_PATTERNS, all_results=True) == rows
    assert filter_cluster_json(rows, _DEFAULT_OPERON_RULES, _REPORT_ALL_PATTERNS,
                               all_results=True) == rows


def test_official_mbn_bytes_cutoffs_and_import_categories():
    registry = {r['stem']: r for r in load_registry(ROOT / 'hmm_library')}
    cutoffs = read_cutoffs(ROOT / 'hmm_library/HMM-bitcutoffs.txt')
    spec = importlib.util.spec_from_file_location('mbn_importer', ROOT / 'scripts/curate_hmm_library.py')
    importer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(importer)
    for profile in json.loads((DATA / 'profile_sources.json').read_text())['profiles']:
        row = registry[profile['stem']]
        raw = (ROOT / 'hmm_library' / row['hmm_file']).read_bytes()
        assert hashlib.sha256(raw).hexdigest() == profile['sha256']
        assert float(row['cutoff']) == profile['sequence_cutoff'] == cutoffs[row['stem']]
        assert float(next(l.split()[1] for l in raw.decode().splitlines() if l.startswith('GA '))) == float(row['cutoff'])
        assert importer.infer_category(row['stem'], row['name'], '')[0] == 'Metallophore-Synthesis'
        assert row['status'] == 'active'


def test_real_mbn_references_reject_related_tgl_and_chr_families(tmp_path):
    exe = shutil.which('hmmsearch')
    if not exe:
        pytest.skip('HMMER is required for accession-verified protein regression')
    registry = {r['stem']: r for r in load_registry(ROOT / 'hmm_library')}
    references = (DATA / 'references.faa').read_text()
    for stem, expected in json.loads((DATA / 'reference_expectations.json').read_text()).items():
        row = registry[stem]
        tbl = tmp_path / (stem + '.tbl')
        subprocess.run([exe, '--cpu', '1', '--noali', '--tblout', str(tbl),
                        str(ROOT / 'hmm_library' / row['hmm_file']), str(DATA / 'references.faa')],
                       check=True, stdout=subprocess.DEVNULL)
        passed = {l.split()[0] for l in tbl.read_text().splitlines()
                  if not l.startswith('#') and float(l.split()[5]) >= float(row['cutoff'])}
        assert set(expected['positive']) <= passed
        assert not set(expected['negative']) & passed
        for acc in expected['positive'] + expected['negative']:
            assert '>' + acc + '\n' in references
