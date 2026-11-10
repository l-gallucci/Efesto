#!/usr/bin/env python3
"""Offline development workflow: grouped nested HMM calibration, never deployment.

Inputs are curated FASTA and TSV (id,label,group,partition,evidence,reference).
Use background for labelled non-targets, unknown for unlabelled candidates;
partition is development, external, or unlabelled. Near-homolog/alias groups
must be assigned before running, independently of labels/model scores.
External and unknown records never enter alignment/build/cutoff selection.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import random
import statistics
import subprocess

BACKGROUND = {'background', 'unknown'}


def fasta(path):
    result = {}; key = None
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            key = line[1:].split()[0]
            if key in result: raise ValueError('Duplicate FASTA ID: ' + key)
            result[key] = ''
        elif line.strip():
            if key is None: raise ValueError('FASTA sequence before header')
            result[key] += line.strip().upper()
    if any(not s or set(s) - set('ACDEFGHIKLMNPQRSTVWYBXZJUO') for s in result.values()):
        raise ValueError('Empty or non-protein sequence')
    return result


def read_inputs(faa, manifest):
    seqs = fasta(faa)
    with Path(manifest).open() as f: rows = list(csv.DictReader(f, delimiter='\t'))
    ids = [r['id'] for r in rows]
    if len(set(ids)) != len(ids): raise ValueError('Duplicate manifest IDs')
    if set(ids) != set(seqs): raise ValueError('FASTA/manifest ID mismatch')
    group_parts = {}; hash_groups = {}; hash_labels = {}
    for r in rows:
        if not r['group']: raise ValueError('Missing preassigned similarity group')
        if r['partition'] not in {'development', 'external', 'unlabelled'}: raise ValueError('Invalid partition')
        if r['partition'] != 'unlabelled' and r['label'] == 'unknown': raise ValueError('Unknown labels cannot be calibration/test truth')
        if r['label'] not in BACKGROUND and not r['reference']: raise ValueError('Target labels need an evidence reference')
        if r['partition']=='development' and r['label'] not in BACKGROUND and set(seqs[r['id']])-set('ACDEFGHIKLMNPQRSTVWYBXZ'):
            raise ValueError('Training seed needs an explicit alignment policy for U/O/J residues')
        h = hashlib.sha256(seqs[r['id']].encode()).hexdigest()
        if r.get('sequence_sha256') and h != r['sequence_sha256']: raise ValueError('Sequence hash mismatch')
        hash_groups.setdefault(h, set()).add(r['group']); hash_labels.setdefault(h, set()).add(r['label'])
        if r['partition'] != 'unlabelled': group_parts.setdefault(r['group'], set()).add(r['partition'])
    if any(len(v) > 1 for v in hash_groups.values()): raise ValueError('Identical sequences assigned to different groups')
    if any(len(v - {'unknown'}) > 1 for v in hash_labels.values()): raise ValueError('Conflicting labels for identical sequences')
    if any(len(v) > 1 for v in group_parts.values()): raise ValueError('Similarity group leaks across development/external')
    # Collapse exact sequence aliases in labelled inputs; preserve IDs elsewhere in manifest.
    seen = set(); selected = []
    for r in sorted(rows, key=lambda r: r['id']):
        h = hashlib.sha256(seqs[r['id']].encode()).hexdigest()
        key = (h, r['partition'], r['label'])
        if key not in seen: selected.append(r); seen.add(key)
    return seqs, selected


def make_folds(rows, n, seed=17):
    groups = {}
    for r in rows: groups.setdefault(r['group'], []).append(r)
    targets = sorted({r['label'] for r in rows} - BACKGROUND)
    support = {c: len({r['group'] for r in rows if r['label'] == c}) for c in targets}
    if n < 2 or any(v < n for v in support.values()):
        raise ValueError('Insufficient independent positive groups for folds: ' + repr(support))
    rng = random.Random(seed); keys = sorted(groups); rng.shuffle(keys)
    # Rare target groups first; each target is balanced before abundant background.
    keys.sort(key=lambda g: min([support[r['label']] for r in groups[g] if r['label'] in targets] or [math.inf]))
    folds = [[] for _ in range(n)]; counts = [dict.fromkeys(targets, 0) for _ in range(n)]
    for g in keys:
        labels = {r['label'] for r in groups[g]} & set(targets)
        k = min(range(n), key=lambda j: (sum(counts[j][c] / support[c] for c in labels), len(folds[j]), j))
        folds[k].extend(groups[g])
        for r in groups[g]:
            if r['label'] in targets: counts[k][r['label']] += 1
    for fold in folds:
        if not all(any(r['label'] == c for r in fold) for c in targets):
            raise ValueError('Cannot stratify mixed-label groups safely; revise fold count/group design')
    return folds


def metrics(truth, predicted):
    tp = sum(t and p for t,p in zip(truth,predicted)); tn = sum(not t and not p for t,p in zip(truth,predicted))
    fp = sum(not t and p for t,p in zip(truth,predicted)); fn = sum(t and not p for t,p in zip(truth,predicted))
    div = lambda a,b: a/b if b else None
    den = (tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)
    return dict(tp=tp,fp=fp,fn=fn,tn=tn,precision=div(tp,tp+fp),recall=div(tp,tp+fn),fpr=div(fp,fp+tn),
                f1=div(2*tp,2*tp+fp+fn),mcc=(tp*tn-fp*fn)/math.sqrt(den) if den else 0.0,mcc_defined=bool(den))


def accepts(score, cutoff, max_evalue):
    return score is not None and score['score'] >= cutoff and score['evalue'] <= max_evalue


def choose_cutoff(rows, scores, target, max_evalue=1e-5):
    vals = sorted({s['score'] for s in scores.values() if s['evalue'] <= max_evalue})
    if not vals: return 0.0, metrics([r['label']==target for r in rows], [False]*len(rows))
    vals.append(math.nextafter(vals[-1], math.inf))
    best = None
    for cutoff in vals:
        m = metrics([r['label']==target for r in rows], [accepts(scores.get(r['id']), cutoff, max_evalue) for r in rows])
        # Deterministic tie: specificity, then lower cutoff. This is predefined,
        # never selected after inspecting the external data.
        key = (m['mcc'], -m['fp'], -cutoff)
        if best is None or key > best[0]: best = (key, cutoff, m)
    return best[1], best[2]


def run_command(cmd, log, stdout=None):
    with Path(log).open('w') as err:
        subprocess.run([str(x) for x in cmd], stdout=stdout or subprocess.DEVNULL, stderr=err, check=True)


def write_fasta(path, rows, seqs):
    Path(path).write_text(''.join('>'+r['id']+'\n'+seqs[r['id']]+'\n' for r in rows))


def build_and_score(train, evaluate, seqs, targets, out, tools):
    out.mkdir(parents=True, exist_ok=False)
    query = out/'evaluate.faa'; write_fasta(query, evaluate, seqs)
    result = {}
    for target in targets:
        positives = [r for r in train if r['label']==target]
        if not positives: raise ValueError('Fold has no target training sequence: ' + target)
        prefix = out/target; seed=prefix.with_suffix('.faa'); aln=prefix.with_suffix('.aln'); hmm=prefix.with_suffix('.hmm')
        write_fasta(seed, positives, seqs)
        with aln.open('w') as f: run_command([tools['mafft'],'--thread','1','--auto','--quiet',seed],str(prefix)+'.mafft.log',f)
        run_command([tools['hmmbuild'],'--cpu','1','--seed','17','-n',target,hmm,aln],str(prefix)+'.hmmbuild.log')
        tbl=prefix.with_suffix('.tbl')
        run_command([tools['hmmsearch'],'--cpu','1','--max','-T','-1000','--domT','-1000','--noali','--tblout',tbl,hmm,query],str(prefix)+'.hmmsearch.log')
        result[target] = {}
        for line in tbl.read_text().splitlines():
            if not line or line.startswith('#'):continue
            c=line.split();result[target][c[0]]=dict(score=float(c[5]),evalue=float(c[4]))
    (out/'provenance.json').write_text(json.dumps(dict(train=[r['id'] for r in train],evaluate=[r['id'] for r in evaluate],train_groups=sorted({r['group'] for r in train}),evaluate_groups=sorted({r['group'] for r in evaluate})),indent=2)+'\n')
    return result


def inner_calibrate(rows, seqs, targets, out, n, tools, max_evalue):
    folds=make_folds(rows,n); scores={t:{} for t in targets}
    for i, held in enumerate(folds):
        held_groups={r['group'] for r in held};train=[r for r in rows if r['group'] not in held_groups]
        fold_scores=build_and_score(train,held,seqs,targets,out/f'inner_{i}',tools)
        for t in targets:scores[t].update(fold_scores[t])
    cutoffs={t:choose_cutoff(rows,scores[t],t,max_evalue)[0] for t in targets}
    return cutoffs,scores


def predictions(rows, scores, cutoffs, max_evalue, method='margin'):
    output=[]
    for r in rows:
        eligible=[(scores[t].get(r['id']),t) for t in cutoffs if accepts(scores[t].get(r['id']),cutoffs[t],max_evalue)]
        if eligible:
            s,t=max(eligible,key=lambda st:((st[0]['score']-cutoffs[st[1]]) if method=='margin' else st[0]['score'],st[1]));pred=t
        else:pred='background'
        output.append(dict(id=r['id'],group=r['group'],truth=r['label'],predicted=pred,evidence=r['evidence']))
    return output


def summarize(pred, targets):
    return {t:metrics([r['truth']==t for r in pred],[r['predicted']==t for r in pred]) for t in targets}


def bootstrap_groups(pred, target, seed=17, iterations=500):
    groups={}
    for r in pred:groups.setdefault(r['group'],[]).append(r)
    keys=sorted(groups);rng=random.Random(seed);values={k:[] for k in ['mcc','precision','recall','fpr']}
    for _ in range(iterations):
        sample=[r for g in rng.choices(keys,k=len(keys)) for r in groups[g]]
        m=summarize(sample,[target])[target]
        for k in values:
            if m[k] is not None and (k!='mcc' or m['mcc_defined']):values[k].append(m[k])
    result={}
    for k,v in values.items():
        v.sort();result[k]=dict(interval95=[v[int(.025*(len(v)-1))],v[int(.975*(len(v)-1))]] if v else None,valid_resamples=len(v))
    return result


def workflow(faa, manifest, output, outer=3, inner=2, mafft='mafft', hmmbuild='hmmbuild',hmmsearch='hmmsearch', max_evalue=1e-5):
    seqs,rows=read_inputs(faa,manifest);dev=[r for r in rows if r['partition']=='development'];external=[r for r in rows if r['partition']=='external'];unknown=[r for r in rows if r['partition']=='unlabelled']
    targets=sorted({r['label'] for r in dev}-BACKGROUND)
    if not targets or any(not all(c.isalnum() or c=='_' for c in t) for t in targets):raise ValueError('Target names must be nonempty alphanumeric/underscore IDs')
    out=Path(output);out.mkdir(parents=True,exist_ok=False);tools=dict(mafft=mafft,hmmbuild=hmmbuild,hmmsearch=hmmsearch)
    folds=make_folds(dev,outer);outer_predictions=[];fold_results=[]
    for i,held in enumerate(folds):
        held_groups={r['group'] for r in held};train=[r for r in dev if r['group'] not in held_groups]
        folder=out/f'outer_{i}';folder.mkdir()
        cuts,_=inner_calibrate(train,seqs,targets,folder,inner,tools,max_evalue)
        scores=build_and_score(train,held,seqs,targets,folder/'held_out',tools)
        pred=predictions(held,scores,cuts,max_evalue);outer_predictions.extend(pred)
        fold_results.append(dict(fold=i,cutoffs=cuts,metrics=summarize(pred,targets),train_groups=sorted({r['group'] for r in train}),test_groups=sorted(held_groups),raw_score_competition=summarize(predictions(held,scores,cuts,max_evalue,'raw'),targets)))
    final_folder=out/'final';final_folder.mkdir();cuts,_=inner_calibrate(dev,seqs,targets,final_folder,inner,tools,max_evalue)
    # Only now score external and unknown records with the frozen final build.
    evaluate=external+unknown
    if not evaluate:evaluate=dev  # build artifacts only; development rescores are never reported as a test
    scores=build_and_score(dev,evaluate,seqs,targets,final_folder/'model',tools)
    external_pred=predictions(external,scores,cuts,max_evalue)
    unknown_pred=predictions(unknown,scores,cuts,max_evalue)
    support={t:dict(development_records=sum(r['label']==t for r in dev),development_groups=len({r['group'] for r in dev if r['label']==t}),external_records=sum(r['label']==t for r in external)) for t in targets}
    reasons=[]
    if any(v['development_groups']<10 for v in support.values()):reasons.append('sparse_independent_positive_groups')
    if not external or any(v['external_records']==0 for v in support.values()):reasons.append('no_independent_external_positive_test')
    if not any(r['evidence']=='experimental_hard_negative' for r in dev):reasons.append('no_experimental_homologous_hard_negatives')
    result=dict(protocol='Efesto hybrid curation/grouped nested HMM calibration v1',objective='MCC',score_type='full_sequence',max_evalue=max_evalue,outer_folds=outer,inner_folds=inner,support=support,folds=fold_results,outer_pooled=summarize(outer_predictions,targets),group_bootstrap={t:bootstrap_groups(outer_predictions,t) for t in targets},outer_predictions=outer_predictions,final_cutoffs=cuts,external_test=summarize(external_pred,targets) if external else None,external_predictions=external_pred,unlabelled_predictions=unknown_pred,promotion=dict(eligible=False,reasons=reasons or ['manual_evidence_and_context_review_required']),inputs={str(faa):hashlib.sha256(Path(faa).read_bytes()).hexdigest(),str(manifest):hashlib.sha256(Path(manifest).read_bytes()).hexdigest()})
    (out/'result.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--fasta',required=True);p.add_argument('--manifest',required=True);p.add_argument('--out',required=True);p.add_argument('--outer',type=int,default=3);p.add_argument('--inner',type=int,default=2)
    p.add_argument('--mafft',default='mafft');p.add_argument('--hmmbuild',default='hmmbuild');p.add_argument('--hmmsearch',default='hmmsearch')
    a=p.parse_args();r=workflow(a.fasta,a.manifest,a.out,a.outer,a.inner,a.mafft,a.hmmbuild,a.hmmsearch)
    print(json.dumps({k:r[k] for k in ['support','outer_pooled','final_cutoffs','external_test','promotion']},indent=2))
if __name__=='__main__':main()
