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
    assert {r['genome']:float(r['topology_weight']) for r in rr}=={'a.faa':1.2,'b.faa':.5}
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
def test_final_confidence_consistent_in_csv_tsv_and_anvio(tmp_path,probability):
    lib=setup_library(tmp_path,'MtrA')
    seed=next(iter(read_fasta(ROOT/'data/seeds/curated/MtrA_seeds.faa').values()))
    (tmp_path/'faa/a.faa').write_text('>ctg_1\n'+seed+'\n>ctg_2\n'+seed+'\n')
    (tmp_path/'gff/a.gff').write_text('##sequence-region ctg 1 10000\nctg\tcaller\tCDS\t100\t300\t.\t+\t0\tID=ctg_1\nctg\tcaller\tCDS\t400\t600\t.\t+\t0\tID=ctg_2\n')
    # Use an independently low UniOP score, even though an operon exists.
    work=tmp_path/'out/_uniop/a';work.mkdir(parents=True)
    (work/'cached.faa').write_text((tmp_path/'faa/a.faa').read_text())
    (work/'uniop.operon').write_text('idx_genes,idx_op\n"[1, 2]",0\n')
    (work/'uniop.pred').write_text(f'Gene A\tGene B\tPrediction\n1\t2\t{probability}\n')
    depth=tmp_path/'depth.tsv';depth.write_text('ctg\t10\n')
    options=('--all_results','--keep_tblout','--anvio','--operon_prediction','--depth',str(depth))
    control=run_cli(tmp_path,lib,*options)
    assert control.returncode==0,control.stderr
    out=tmp_path/'out'
    tsv=list(csv.DictReader((out/'Efesto-results-long.tsv').open(),delimiter='\t'))
    baseline={r['orf']:r['cluster_confidence'] for r in tsv}
    assert len(baseline)==2
    if probability:
        assert all(0 < float(v) < probability for v in baseline.values())
    else:
        assert all(float(v)==0.0 for v in baseline.values())
    for name in ('Efesto-summary.csv','Efesto-geneSummary-clusters.csv'):
        rr=[r for r in csv.DictReader((out/name).open()) if r.get('orf') in baseline]
        assert {r['orf']:r['cluster_confidence'] for r in rr}==baseline
    aa=list(csv.DictReader((out/'Efesto-anvio-gene-scores.tsv').open(),delimiter='\t'))
    assert all(float(r['uniop_weight'])==probability for r in aa)
    assert all(abs(float(r['cluster_confidence'])-float(baseline[r['gene_callers_id']]))<.001 for r in aa)
    normalised=run_cli(tmp_path,lib,*options,'--norm_coverage')
    assert normalised.returncode==0,normalised.stderr
    rr=list(csv.DictReader((out/'Efesto-results-long.tsv').open(),delimiter='\t'))
    assert {r['orf']:r['cluster_confidence'] for r in rr}==baseline
    assert 'coverage_per_million_sum' in (out/'Efesto-coverage-heatmap.csv').read_text()
