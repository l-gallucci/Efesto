# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Efesto contributors.
"""UniOP cache provenance and Python scripts discovered through PATH."""
import json
import os
import pytest
from efesto.uniop import run_uniop, make_uniop_faa
from efesto.eggnog import run_eggnog
from tests._tool_stubs import write_uniop_stub, invocation_count


@pytest.fixture
def uniop_case(tmp_path):
    faa=tmp_path/'sample.faa';faa.write_text('>a # 1 # 90 # 1 # ID=1\nMAA\n>b # 100 # 190 # 1 # ID=2\nMKK\n')
    tool=write_uniop_stub(tmp_path/'UniOP')
    def run(**options):
        return run_uniop([faa],None,tmp_path,str(tool),prodigal_faa_dir=tmp_path,
                         return_status=True,**options)
    return faa,tool,run


def test_uniop_reuses_only_the_exact_ordered_input(uniop_case):
    faa,tool,run=uniop_case
    assert run()[2][faa.name]=='completed'
    assert run()[2][faa.name]=='cached' and invocation_count(tool)==1
    faa.write_text('>new_B # 100 # 190 # 1 # ID=2\nMKK\n>new_A # 1 # 90 # 1 # ID=1\nMAA\n')
    ops,probs,status=run()
    assert status[faa.name]=='completed' and invocation_count(tool)==2
    assert set(ops[faa.name])=={'new_A','new_B'}
    assert probs[faa.name]=={('new_A','new_B'):0.8}
    # An unrelated FAA in the work directory cannot replace the preserved index.
    work=faa.parent/'_uniop/sample';(work/'aaa.faa').write_text('>wrong\nM\n')
    assert run()[0]==ops and invocation_count(tool)==2
    (work/'uniop.pred').write_text('partial')
    assert run()[2][faa.name]=='completed' and invocation_count(tool)==3
    tool.write_text(tool.read_text()+'\n# new tool version\n')
    assert run()[2][faa.name]=='completed' and invocation_count(tool)==4


def test_uniop_legacy_or_failed_output_is_never_reused(uniop_case):
    faa,tool,run=uniop_case
    work=faa.parent/'_uniop/sample';work.mkdir(parents=True)
    (work/'uniop.operon').write_text('idx_genes,idx_op\n"[1,2]",99\n')
    (work/'uniop.pred').write_text('Gene A\tGene B\tPrediction\n1\t2\t0.99\n')
    tool.with_suffix('.control.json').write_text(json.dumps({'exit':2}))
    for attempt in (1,2):
        ops,probs,status=run()
        assert ops==probs=={} and status[faa.name]=='failed'
        assert invocation_count(tool)==attempt
        assert not (work/'efesto-cache.json').exists()
    tool.with_suffix('.control.json').write_text('{}')
    assert run()[1][faa.name][('a','b')]==0.8 and invocation_count(tool)==3


def test_bakta_coordinates_invalidate_prepared_faa(tmp_path):
    faa=tmp_path/'sample.faa';faa.write_text('>a\nMAA\n>b\nMKK\n')
    gff=tmp_path/'sample.gff'
    def coordinates(start):
        gff.write_text(f'ctg\tBakta\tCDS\t{start}\t90\t.\t+\t0\tID=a\nctg\tBakta\tCDS\t100\t190\t.\t+\t0\tID=b\n')
    coordinates(1);tool=write_uniop_stub(tmp_path/'UniOP')
    def run():return run_uniop([faa],None,tmp_path,str(tool),bakta_gff_dir=tmp_path,return_status=True)
    assert run()[2][faa.name]=='completed'
    coordinates(12)
    assert run()[2][faa.name]=='completed' and invocation_count(tool)==2
    assert '>a # 12 # 90 # 1 #' in (tmp_path/'_uniop/sample/received.txt').read_text()
    faa.write_text(faa.read_text()+'>missing\nMAA\n')
    assert run()[2][faa.name]=='input_unavailable' and invocation_count(tool)==2
    assert not make_uniop_faa(faa,gff,tmp_path/'partial.faa')


def test_fna_fallback_preserves_the_generated_index(tmp_path):
    faa=tmp_path/'sample.faa';faa.write_text('>internal_A\nMAA\n>internal_B\nMKK\n')
    fna=tmp_path/'sample.fna';fna.write_text('>ctg\nATGATG\n')
    tool=write_uniop_stub(tmp_path/'UniOP')
    def run():return run_uniop([faa],tmp_path,tmp_path,str(tool),return_status=True)
    first=run();assert first[2][faa.name]=='completed'
    assert set(first[0][faa.name])=={'internal_A','internal_B'}
    assert run()[2][faa.name]=='cached' and invocation_count(tool)==1
    fna.write_text('>ctg\nATGATGATG\n')
    assert run()[2][faa.name]=='completed' and invocation_count(tool)==2


def test_python_tools_are_resolved_through_path(tmp_path,monkeypatch):
    bindir=tmp_path/'bin';bindir.mkdir()
    monkeypatch.setenv('PATH',str(bindir)+os.pathsep+os.environ['PATH'])
    tool=write_uniop_stub(bindir/'UniOP')
    faa=tmp_path/'sample.faa';faa.write_text('>a\nMAA\n>b\nMKK\n')
    assert run_uniop([faa],None,tmp_path,'UniOP',prodigal_faa_dir=tmp_path)[1]
    assert invocation_count(tool)==1
    egg=bindir/'emapper.py'
    egg.write_text('import sys\nfrom pathlib import Path\nout=Path(sys.argv[sys.argv.index("--output_dir")+1]);prefix=sys.argv[sys.argv.index("-o")+1]\n(out/(prefix+".emapper.annotations")).write_text("#query\\tPreferred_name\\tKEGG_ko\\na\\tsufA\\tK07168\\n")\n')
    egg.chmod(0o755)
    assert run_eggnog(faa,tmp_path/'eggnog',eggnog_path='emapper.py')
    assert run_eggnog(faa,tmp_path/'missing',eggnog_path='missing-emapper.py') is None
