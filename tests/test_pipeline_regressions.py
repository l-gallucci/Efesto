# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Efesto contributors.

"""Small complete CLI runs for the pre-release identity/filter/error fixes."""
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import pytest
from efesto.io import load_registry, read_fasta

ROOT=Path(__file__).resolve().parents[1]


def setup_library(tmp_path, stem):
    row=next(r for r in load_registry(ROOT/'hmm_library') if r['stem']==stem)
    lib=tmp_path/'library';(lib/row['category']).mkdir(parents=True)
    shutil.copy2(ROOT/'hmm_library'/row['hmm_file'],lib/row['hmm_file'])
    with (lib/'hmm_registry.tsv').open('w') as h:
        w=csv.DictWriter(h,fieldnames=list(row),delimiter='\t');w.writeheader();w.writerow(row)
    (lib/'HMM-bitcutoffs.txt').write_text(stem+'\t'+row['cutoff']+'\n')
    (lib/'FeGenie-map.txt').write_text(stem+'\t'+stem+'\n')
    (tmp_path/'faa').mkdir();(tmp_path/'gff').mkdir()
    return lib


def run_cli(tmp_path, lib, *extra, path_prefix=None):
    env=os.environ.copy();env['PYTHONPATH']=str(ROOT/'src')
    if path_prefix:env['PATH']=str(path_prefix)+os.pathsep+env['PATH']
    return subprocess.run([sys.executable,'-m','efesto.cli','--faa_dir',str(tmp_path/'faa'),
                           '--gff_dir',str(tmp_path/'gff'),'--hmm_dir',str(lib),
                           '--out',str(tmp_path/'out'),'--threads','1',*extra],
                          capture_output=True,text=True,env=env)


@pytest.mark.skipif(not shutil.which('hmmsearch'),reason='HMMER required')
def test_bakta_short_contig_filtered_by_gff_contig(tmp_path):
    lib=setup_library(tmp_path,'lanM_NF053806')
    panel=ROOT/'data/validation/metal_expansion'
    key=json.loads((panel/'reference_expectations.json').read_text())['lanM_NF053806']['positive'][0]
    seq=read_fasta(panel/'references.faa')[key]
    (tmp_path/'faa/sample.faa').write_text('>BAKTA_00001\n'+seq+'\n')
    (tmp_path/'gff/sample.gff').write_text('##gff-version 3\n##sequence-region short_contig 1 1000\nshort_contig\tBakta\tCDS\t100\t'+str(99+len(seq)*3)+'\t.\t+\t0\tID=BAKTA_00001\n')
    control=run_cli(tmp_path,lib,'--all_results','--keep_tblout')
    assert control.returncode==0,control.stderr
    path=tmp_path/'out/Efesto-results-long.tsv'
    assert len(list(csv.DictReader(path.open(),delimiter='\t')))==1
    filtered=run_cli(tmp_path,lib,'--all_results','--keep_tblout','--min_contig_len','2000')
    assert filtered.returncode==0,filtered.stderr
    assert list(csv.DictReader(path.open(),delimiter='\t'))==[]
    assert '1 ORFs removed (short contig)' in filtered.stdout


@pytest.mark.skipif(not shutil.which('hmmsearch'),reason='HMMER required')
def test_two_genomes_same_orf_receive_independent_topology(tmp_path):
    lib=setup_library(tmp_path,'MtrA')
    seed=next(iter(read_fasta(ROOT/'data/seeds/curated/MtrA_seeds.faa').values()))
    for name,seq in [('a',seed),('b',seed+'AAA')]:
        (tmp_path/f'faa/{name}.faa').write_text('>contig_1\n'+seq+'\n')
        (tmp_path/f'gff/{name}.gff').write_text('ctg\tcaller\tCDS\t5000\t'+str(4999+3*len(seq))+'\t.\t+\t0\tID=contig_1\n')
    stub=tmp_path/'phobius'
    stub.write_text('#!'+sys.executable+'\nimport sys\nfrom pathlib import Path\nrecords=Path(sys.argv[-1]).read_text().split(">")\n'+
                    'print("SEQENCE_ID TM SP PREDICTION")\nfor record in records:\n if not record:continue\n lines=record.splitlines();key=lines[0].split()[0];seq="".join(lines[1:])\n print(key,"0 Y test" if len(seq)=='+str(len(seed))+' else "4 0 test")\n')
    stub.chmod(0o755)
    result=run_cli(tmp_path,lib,'--all_results','--topology_check','--phobius_path',str(stub),'--anvio')
    assert result.returncode==0,result.stderr
    rr=list(csv.DictReader((tmp_path/'out/Efesto-results-long.tsv').open(),delimiter='\t'))
    assert {r['genome']:r['topology_verdict'] for r in rr}=={'a.faa':'confirmed','b.faa':'contradicted'}
    assert (tmp_path/'out/efesto_topology_candidates.faa').read_text().count('>')==2
    aa=list(csv.DictReader((tmp_path/'out/Efesto-anvio-functions.tsv').open(),delimiter='\t'))
    assert len(aa)==2 and len({r['gene_callers_id'] for r in aa})==2


def test_failed_hmmer_stops_cli_and_is_retried(tmp_path):
    lib=setup_library(tmp_path,'MtrA')
    (tmp_path/'faa/a.faa').write_text('>orf_1\nMAAAA\n')
    bindir=tmp_path/'bin';bindir.mkdir()
    stub=bindir/'hmmsearch'
    stub.write_text('#!'+sys.executable+'\nimport sys\nfrom pathlib import Path\np=Path(sys.argv[sys.argv.index("--tblout")+1]);p.write_text("# partial table\\n")\nprint("simulated HMMER failure",file=sys.stderr)\nsys.exit(2)\n')
    stub.chmod(0o755)
    for _ in range(2):
        result=run_cli(tmp_path,lib,'--keep_tblout',path_prefix=bindir)
        assert result.returncode!=0
        assert 'simulated HMMER failure' in result.stderr
        assert not (tmp_path/'out/Efesto-results-long.tsv').exists()
        assert not list((tmp_path/'out/_tblout_cache').glob('*.tblout'))


@pytest.mark.skipif(not shutil.which('hmmsearch'),reason='HMMER required')
@pytest.mark.parametrize('probability',[0.0,0.2])
def test_raw_uniop_pairs_and_hit_evidence_survive_normalisation(tmp_path,probability):
    lib=setup_library(tmp_path,'MtrA')
    seed=next(iter(read_fasta(ROOT/'data/seeds/curated/MtrA_seeds.faa').values()))
    (tmp_path/'faa/a.faa').write_text('>ctg_1\n'+seed+'\n>ctg_2\n'+seed+'\n')
    (tmp_path/'gff/a.gff').write_text('##sequence-region ctg 1 10000\nctg\tcaller\tCDS\t100\t300\t.\t+\t0\tID=ctg_1\nctg\tcaller\tCDS\t400\t600\t.\t+\t0\tID=ctg_2\n')
    from tests._tool_stubs import write_uniop_stub
    tool=write_uniop_stub(tmp_path/'UniOP',probability=probability)
    depth=tmp_path/'depth.tsv';depth.write_text('ctg\t10\n')
    options=('--all_results','--keep_tblout','--anvio','--operon_prediction','--uniop_path',str(tool),'--bakta_gff_dir',str(tmp_path/'gff'),'--depth',str(depth))
    result=run_cli(tmp_path,lib,*options)
    assert result.returncode==0,result.stderr
    out=tmp_path/'out'
    rows=list(csv.DictReader((out/'Efesto-results-long.tsv').open(),delimiter='\t'))
    assert len(rows)==2
    assert all(r['threshold_type']=='model_specific' for r in rows)
    assert all(r['uniop_status']=='completed' for r in rows)
    pairs=list(csv.DictReader((out/'Efesto-UniOP-pair-evidence.tsv').open(),delimiter='\t'))
    assert len(pairs)==1 and float(pairs[0]['operon_probability'])==probability
    for name in ('Efesto-results-long.tsv','Efesto-summary.csv','Efesto-geneSummary-clusters.csv',
                 'Efesto-hits.gff3','Efesto-hits.faa','Efesto-summary-stats.tsv',
                 'Efesto-anvio-gene-evidence.tsv'):
        content=(out/name).read_text()
        assert 'cluster_confidence' not in content
        assert 'co_occ_score' not in content
    normalised=run_cli(tmp_path,lib,*options,'--norm_coverage')
    assert normalised.returncode==0,normalised.stderr
    repeated=list(csv.DictReader((out/'Efesto-results-long.tsv').open(),delimiter='\t'))
    assert all(r['uniop_status']=='cached' for r in repeated)
    assert [{**r,'uniop_status':'cached'} for r in rows]==repeated
    assert 'coverage_per_million_sum' in (out/'Efesto-coverage-heatmap.csv').read_text()


@pytest.mark.skipif(not shutil.which('hmmsearch'),reason='HMMER required')
@pytest.mark.parametrize('supported',[False,True])
def test_confirmation_candidates_report_dropped_and_rescued_without_changing_hit_scope(tmp_path,supported):
    lib=setup_library(tmp_path,'MtrA')
    registry=lib/'hmm_registry.tsv'
    with registry.open() as handle:
        reader=csv.DictReader(handle,delimiter='\t');fields=reader.fieldnames;row=next(reader)
    row['needs_confirmation']='yes';row['kegg_ko']='K00001'
    with registry.open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t');writer.writeheader();writer.writerow(row)
    (lib/'operon_rules.json').write_text(json.dumps({'rules':[{'name':'minimum2','genes':['MtrA','other'],
        'rule':'require_n_of','min_genes':2,'on_fail':'drop'}]}))
    seed=next(iter(read_fasta(ROOT/'data/seeds/curated/MtrA_seeds.faa').values()))
    (tmp_path/'faa/a.faa').write_text('>ctg_1\n'+seed+'\n')
    (tmp_path/'gff/a.gff').write_text('ctg\tcaller\tCDS\t5000\t7000\t.\t+\t0\tID=ctg_1\n')
    annotations=tmp_path/'annotations.tsv'
    annotations.write_text('#query\tPreferred_name\tKEGG_ko\nctg_1\t'+('MtrA\tK00001' if supported else 'other\tK99999')+'\n')
    result=run_cli(tmp_path,lib,'--eggnog_annotations',str(annotations))
    assert result.returncode==0,result.stderr
    out=tmp_path/'out'
    retained=list(csv.DictReader((out/'Efesto-results-long.tsv').open(),delimiter='\t'))
    candidates=list(csv.DictReader((out/'Efesto-confirmation-candidates.tsv').open(),delimiter='\t'))
    assert len(candidates)==1 and len(retained)==int(supported)
    assert candidates[0]['eggnog_verdict']==('confirmed' if supported else 'contradicted')
    assert candidates[0]['retention_mode']==('eggnog_rescue' if supported else 'not_retained')
    if supported:
        assert retained[0]['cluster_id']==candidates[0]['cluster_id']
        assert retained[0]['original_orf']==candidates[0]['original_orf']=='ctg_1'
    else:
        stats=(out/'Efesto-summary-stats.tsv').read_text()
        assert 'RUN\ttotal_orf_hits\t0' in stats


@pytest.mark.skipif(not shutil.which('hmmsearch'),reason='HMMER required')
def test_fna_run_ignores_orphan_genomes_and_stops_on_prediction_failure(tmp_path):
    lib=setup_library(tmp_path,'MtrA')
    seed=next(iter(read_fasta(ROOT/'data/seeds/curated/MtrA_seeds.faa').values()))
    fna=tmp_path/'fna';fna.mkdir();(fna/'sample.fna').write_text('>ctg\n'+'ATG'*400+'\n')
    stub_dir=tmp_path/'python-stubs';stub_dir.mkdir()
    (stub_dir/'pyrodigal.py').write_text('''__version__='test'
class Gene:
    begin=1;end=1200;strand=1;start_type='ATG';rbs_motif=None;partial_end=False
    def translate(self):return SEED
class GeneFinder:
    def __init__(self,**kwargs):pass
    def find_genes(self,seq):
        if 'CCC' in seq:raise ValueError('controlled prediction failure')
        return [Gene()]
'''.replace('SEED',repr(seed)))
    orphan=tmp_path/'out/_prodigal/faa/orphan.faa';orphan.parent.mkdir(parents=True)
    orphan.write_text('>orphan_1\n'+seed+'\n')
    env=os.environ.copy();env['PYTHONPATH']=str(stub_dir)+os.pathsep+str(ROOT/'src')
    cmd=[sys.executable,'-m','efesto.cli','--fna_dir',str(fna),'--meta','--gene_caller','pyrodigal',
         '--hmm_dir',str(lib),'--out',str(tmp_path/'out'),'--threads','1','--all_results']
    first=subprocess.run(cmd,env=env,capture_output=True,text=True)
    assert first.returncode==0,first.stderr
    output=tmp_path/'out/Efesto-results-long.tsv'
    rows=list(csv.DictReader(output.open(),delimiter='\t'))
    assert len(rows)==1 and rows[0]['genome']=='sample.faa'
    assert 'Genomes input : 1' in (tmp_path/'out/Efesto-run.log').read_text()
    # A failed prediction must stop before any new annotation table is emitted.
    output.unlink();(fna/'sample.fna').write_text('>ctg\n'+'CCC'*400+'\n')
    failed=subprocess.run(cmd,env=env,capture_output=True,text=True)
    assert failed.returncode!=0 and 'controlled prediction failure' in failed.stderr
    assert not output.exists()
    # The same failed input is retried rather than accepted as cached success.
    retry=subprocess.run(cmd,env=env,capture_output=True,text=True)
    assert retry.returncode!=0 and 'controlled prediction failure' in retry.stderr
