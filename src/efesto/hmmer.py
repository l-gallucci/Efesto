"""HMMER search execution, tblout parsing, and library normalization."""

import os
import sys
import subprocess
import tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

_CURRENT_HMM_TAG = "HMMER3/f"


def normalize_hmm_library(hmm_dir):
    """Convert any pre-HMMER3/f profiles to current format using hmmconvert.

    Safe to call repeatedly — skips files already in HMMER3/f format.
    Converts in-place. Prints a one-line summary.
    """
    hmm_dir = Path(hmm_dir)
    # Only category dirs — skip work/staging dirs (._ prefixed) so their scratch
    # HMMs are never mutated or counted.
    all_hmms = sorted(
        hf for d in hmm_dir.iterdir()
        if d.is_dir() and not d.name.startswith((".", "_"))
        for hf in d.glob("*.hmm")
    )
    needs = []
    for hmm in all_hmms:
        try:
            with open(hmm) as fh:
                first = next((l.strip() for l in fh if l.strip()), "")
        except OSError:
            continue
        if not first.startswith(_CURRENT_HMM_TAG):
            needs.append(hmm)

    if not needs:
        print(f"[INFO] HMM library already normalized ({len(all_hmms)} files)")
        return

    print(f"[INFO] Normalizing {len(needs)}/{len(all_hmms)} HMM profiles "
          f"to HMMER3/f via hmmconvert…")
    converted, failed = 0, []
    for hmm in needs:
        r = subprocess.run(["hmmconvert", str(hmm)],
                           capture_output=True, text=True)
        if r.returncode == 0 and r.stdout.strip():
            hmm.write_text(r.stdout)
            converted += 1
        else:
            failed.append(hmm.name)

    if failed:
        print(f"[WARN] hmmconvert failed for: {', '.join(failed)}", file=sys.stderr)
    print(f"[INFO] Normalized {converted} profiles"
          + (f"  ({len(failed)} failed)" if failed else ""))


def _hmmsearch_job(args_tuple):
    hmm_file, faa_file, tblout_path, bitscore_cutoff, threads, fallback_bitscore = args_tuple
    if Path(tblout_path).exists():
        return tblout_path, True, ""
    effective_cutoff = bitscore_cutoff if bitscore_cutoff > 0 else fallback_bitscore
    cmd = ["hmmsearch", "--cpu", str(threads),
           "-T", str(max(effective_cutoff, 0)),
           "--tblout", tblout_path, "--noali",
           "-o", "/dev/null", hmm_file, faa_file]
    r = subprocess.run(cmd, capture_output=True, text=True)
    return tblout_path, r.returncode == 0, r.stderr if r.returncode else ""


def parse_tblout(path, max_evalue=1e-5):
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
            if e <= max_evalue:
                hits.append((parts[0], e, b))
    return hits


def run_all_hmmsearches(faa_files, cat_hmms, cutoffs, out_tmp,
                        threads_total, hmm_threads=1, zero_cutoff_min_bitscore=30):
    jobs = [
        (str(hmm_path), str(faa),
         str(out_tmp / f"{faa.name}__{stem}.tblout"),
         cutoffs.get(stem, 0), hmm_threads, zero_cutoff_min_bitscore)
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


def collect_best_hits(faa_files, cat_hmms, out_tmp, cutoffs=None, max_evalue=1e-5):
    if cutoffs is None:
        cutoffs = {}
    bh = defaultdict(dict)
    for faa in faa_files:
        genome = faa.name
        for cat, hmm_list in cat_hmms.items():
            for stem, _ in hmm_list:
                calibrated = cutoffs.get(stem, 0) > 0
                for orf, ev, bs in parse_tblout(
                        str(out_tmp / f"{genome}__{stem}.tblout"),
                        max_evalue=max_evalue):
                    prev = bh[genome].get(orf)
                    if prev is None or bs > prev["bitscore"]:
                        bh[genome][orf] = {
                            "hmm_stem": stem, "cat": cat,
                            "evalue": ev, "bitscore": bs,
                            "confidence": "calibrated" if calibrated else "low_confidence",
                        }
    return bh
