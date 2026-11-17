# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Efesto contributors.

"""HMMER search execution, tblout parsing, and library normalization."""

import os
import sys
import subprocess
import tempfile
import hashlib
import json
import shutil
from functools import lru_cache
from collections import defaultdict
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from itertools import islice
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


class HMMSearchError(RuntimeError):
    """One or more required searches failed; results must not be published."""


@lru_cache(maxsize=4096)
def _digest_at_stat(path, size, mtime_ns, ctime_ns):
    # The stat tuple invalidates this memo within a run; large FAA inputs are
    # hashed once rather than once for each profile.
    with open(path, "rb") as handle:
        digest = hashlib.sha256()
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
        return digest.hexdigest()


def _fingerprint(path):
    path = Path(path).resolve()
    stat = path.stat()
    return {"path": str(path), "sha256": _digest_at_stat(
        str(path), stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)}


def _cache_signature(hmm, faa, cutoff, threads):
    binary = shutil.which("hmmsearch")
    if not binary:
        raise FileNotFoundError("hmmsearch executable not found in PATH")
    return {"schema": 1, "hmm": _fingerprint(hmm), "faa": _fingerprint(faa),
            "binary": _fingerprint(binary), "cutoff": cutoff, "threads": threads}


def _complete_table(path):
    if not path.is_file():
        return False
    with path.open("rb") as handle:
        handle.seek(max(0, path.stat().st_size - 1024))
        return handle.read().rstrip().endswith(b"# [ok]")


def _hmmsearch_job(args_tuple):
    hmm_file, faa_file, tblout_path, bitscore_cutoff, threads, fallback_bitscore = args_tuple
    target = Path(tblout_path)
    metadata = target.with_name(target.name + ".cache.json")
    temporary = None
    effective_cutoff = max(bitscore_cutoff if bitscore_cutoff > 0 else fallback_bitscore, 0)
    try:
        signature = _cache_signature(hmm_file, faa_file, effective_cutoff, threads)
        try:
            cached = json.loads(metadata.read_text())
            if (isinstance(cached, dict) and _complete_table(target) and cached.get("inputs") == signature
                    and cached.get("tblout_sha256") == _fingerprint(target)["sha256"]):
                return tblout_path, True, ""
        except (OSError, ValueError):
            pass
        # Invalidate stale output before launching: failures cannot leave a
        # final file that a later collector might mistake for a valid result.
        target.unlink(missing_ok=True)
        metadata.unlink(missing_ok=True)
        with tempfile.NamedTemporaryFile(dir=target.parent, prefix=target.name + ".",
                                         suffix=".tmp", delete=False) as handle:
            temporary = Path(handle.name)
        cmd = ["hmmsearch", "--cpu", str(threads),
               "-T", str(effective_cutoff), "--tblout", str(temporary), "--noali",
               "-o", os.devnull, hmm_file, faa_file]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            return tblout_path, False, f"{hmm_file} / {faa_file}: exit {result.returncode}: {result.stderr.strip()}"
        if not _complete_table(temporary):
            return tblout_path, False, f"{hmm_file} / {faa_file}: incomplete HMMER tblout (missing # [ok])"
        digest = _fingerprint(temporary)["sha256"]
        temporary.replace(target)
        with tempfile.NamedTemporaryFile(mode="w", dir=target.parent,
                                         prefix=metadata.name + ".", suffix=".tmp",
                                         delete=False) as handle:
            temporary = Path(handle.name)
            json.dump({"inputs": signature, "tblout_sha256": digest}, handle)
        temporary.replace(metadata)
        return tblout_path, True, ""
    except OSError as error:
        return tblout_path, False, f"{hmm_file} / {faa_file}: {error}"
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


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
    # Python only launches/waits for external HMMER processes. Threads avoid
    # an extra Python interpreter per concurrent search (or inherited parent
    # memory on fork-based platforms).
    if threads_total < 1 or hmm_threads < 1:
        raise ValueError("Thread counts must be positive")
    faa_files = tuple(faa_files)
    jobs = (
        (str(hmm_path), str(faa),
         str(out_tmp / f"{faa.name}__{stem}.tblout"),
         cutoffs.get(stem, 0), hmm_threads, zero_cutoff_min_bitscore)
        for faa in faa_files
        for cat, hmm_list in cat_hmms.items()
        for stem, hmm_path in hmm_list
    )
    total = len(faa_files) * sum(len(hmms) for hmms in cat_hmms.values())
    # HMMER's per-search working matrices dominate peak RSS. Limit their
    # overlap by default, while respecting a smaller requested CPU budget.
    n_workers = min(2, max(1, threads_total // hmm_threads), max(1, total))
    print(f"[INFO] HMMsearch concurrency: {n_workers} searches, "
          f"{hmm_threads} HMMER workers each")
    done = errors = 0
    with ThreadPoolExecutor(max_workers=n_workers) as pool:
        # Keep at most two jobs per worker in flight, rather than retaining
        # a Future and a job tuple for every genome/profile combination.
        pending = {pool.submit(_hmmsearch_job, job)
                   for job in islice(jobs, 2 * n_workers)}
        while pending:
            completed, pending = wait(pending, return_when=FIRST_COMPLETED)
            for fut in completed:
                done += 1
                _, ok, err = fut.result()
                if not ok:
                    errors += 1
                    print(f"\n  [WARN] hmmsearch: {err[:4000]}", file=sys.stderr)
                sys.stdout.write(
                    f"\r[INFO] hmmsearch: {done}/{total}  ({errors} errors)  ")
                sys.stdout.flush()
                job = next(jobs, None)
                if job is not None:
                    pending.add(pool.submit(_hmmsearch_job, job))
    print()
    if errors:
        raise HMMSearchError(f"{errors}/{total} HMMER searches failed; no final annotation was generated. "
                             "Fix the reported errors and rerun; valid searches can be reused with --keep_tblout.")


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
                            "threshold_type": "model_specific" if calibrated else "generic",
                        }
    return bh
