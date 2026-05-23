"""HMMER search execution and tblout parsing."""

import os
import sys
import subprocess
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path


def _hmmsearch_job(args_tuple):
    hmm_file, faa_file, tblout_path, bitscore_cutoff, threads = args_tuple
    if Path(tblout_path).exists():
        return tblout_path, True, ""
    cmd = ["hmmsearch", "--cpu", str(threads),
           "-T", str(max(bitscore_cutoff, 0)),
           "--tblout", tblout_path, "--noali",
           "-o", "/dev/null", hmm_file, faa_file]
    r = subprocess.run(cmd, capture_output=True, text=True)
    return tblout_path, r.returncode == 0, r.stderr if r.returncode else ""


def parse_tblout(path):
    hits = []
    if not os.path.isfile(path):
        return hits
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            try:
                e = float(parts[4])
                b = float(parts[5])
            except ValueError:
                continue
            if e < 0.1:
                hits.append((parts[0], e, b))
    return hits


def run_all_hmmsearches(faa_files, cat_hmms, cutoffs, out_tmp,
                        threads_total, hmm_threads=1):
    jobs = [
        (str(hmm_path), str(faa),
         str(out_tmp / f"{faa.name}__{stem}.tblout"),
         cutoffs.get(stem, 0), hmm_threads)
        for faa in faa_files
        for cat, hmm_list in cat_hmms.items()
        for stem, hmm_path in hmm_list
    ]
    n_workers = max(1, threads_total // hmm_threads)
    total = len(jobs)
    done = errors = 0
    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        futures = {pool.submit(_hmmsearch_job, j): j for j in jobs}
        for fut in as_completed(futures):
            done += 1
            _, ok, err = fut.result()
            if not ok:
                errors += 1
                print(f"\n  [WARN] hmmsearch: {err[:80]}", file=sys.stderr)
            sys.stdout.write(
                f"\r[INFO] hmmsearch: {done}/{total}  ({errors} errors)  ")
            sys.stdout.flush()
    print()


def collect_best_hits(faa_files, cat_hmms, out_tmp):
    bh = defaultdict(dict)
    for faa in faa_files:
        genome = faa.name
        for cat, hmm_list in cat_hmms.items():
            for stem, _ in hmm_list:
                for orf, ev, bs in parse_tblout(
                        str(out_tmp / f"{genome}__{stem}.tblout")):
                    prev = bh[genome].get(orf)
                    if prev is None or bs > prev["bitscore"]:
                        bh[genome][orf] = {
                            "hmm_stem": stem, "cat": cat,
                            "evalue": ev, "bitscore": bs,
                        }
    return bh
