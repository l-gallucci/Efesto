#!/usr/bin/env python3
"""
layer_b_dedup_all.py
====================
Layer B deduplication: sequence-level clustering of all active HMM profiles
using hmmemit consensus sequences + MMseqs2.

Identifies redundant models (≥70% identity, ≥80% bidirectional coverage),
selects a representative per cluster, and marks others as
'deprecated_dedup_sequence_cluster' in hmm_registry.tsv.

Usage:
    python scripts/layer_b_dedup_all.py [options]

Options:
    --hmm_dir     HMM library root (default: hmm_library/)
    --work_dir    Scratch directory for intermediates (default: hmm_library/_dedup_work/)
    --min_id      Minimum sequence identity for clustering (default: 0.70)
    --min_cov     Minimum bidirectional coverage (default: 0.80)
    --dry_run     Report clusters without modifying registry
    --skip_emit   Skip hmmemit step; reuse existing all_consensus.faa in work_dir
"""

import argparse
import csv
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path


# ── Representative selection ──────────────────────────────────────────────────
# Priority: calibrated cutoff > nseq > source tier
_SOURCE_TIER = {"interpro": 3, "fegenie": 2, "tabuteau": 2, "methmmdb": 1}


def _rep_score(row):
    calibrated = 1 if float(row.get("cutoff", 0) or 0) > 0 else 0
    try:
        nseq = int(row.get("nseq", 0) or 0)
    except ValueError:
        nseq = 0
    tier = _SOURCE_TIER.get(row.get("source", ""), 0)
    return (calibrated, nseq, tier)


# ── Registry I/O ──────────────────────────────────────────────────────────────

def load_registry(hmm_dir):
    path = Path(hmm_dir) / "hmm_registry.tsv"
    if not path.exists():
        sys.exit(f"[ERROR] Registry not found: {path}")
    rows = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = reader.fieldnames
        for row in reader:
            rows.append(row)
    return rows, fieldnames


def save_registry(hmm_dir, rows, fieldnames):
    path = Path(hmm_dir) / "hmm_registry.tsv"
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t",
                           extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


# ── hmmemit step ─────────────────────────────────────────────────────────────

def build_consensus_faa(active_rows, faa_path, hmm_dir):
    """Run hmmemit -c for each active model; write pooled FAA.

    Sequence header format: {category}__{stem}
    (double underscore chosen because single underscore appears in stems)
    """
    hmm_dir = Path(hmm_dir)
    written = 0
    failed = []
    with open(faa_path, "w") as out:
        for row in active_rows:
            stem = row["stem"]
            cat  = row.get("category", "unknown")
            # Prefer the local repo path (portable); fall back to registry path
            local = hmm_dir / cat / f"{stem}.hmm"
            hmm_file = str(local) if local.exists() else row.get("hmm_file", "")
            if not hmm_file or not Path(hmm_file).exists():
                failed.append(stem)
                continue
            r = subprocess.run(
                ["hmmemit", "-c", hmm_file],
                capture_output=True, text=True)
            if r.returncode != 0 or not r.stdout.strip():
                failed.append(stem)
                continue
            # Replace the default ">modelname-consensus" header
            seq_lines = r.stdout.splitlines()
            header = f">{cat}__{stem}"
            out.write(header + "\n")
            out.write("\n".join(
                l for l in seq_lines if not l.startswith(">")) + "\n")
            written += 1
    if failed:
        print(f"[WARN] hmmemit failed for {len(failed)} models: "
              f"{', '.join(failed[:10])}" +
              (f" …and {len(failed)-10} more" if len(failed) > 10 else ""),
              file=sys.stderr)
    print(f"[INFO] hmmemit: {written} consensus sequences written → {faa_path}")
    return written


# ── MMseqs2 step ──────────────────────────────────────────────────────────────

def run_mmseqs(faa_path, work_dir, min_id, min_cov):
    prefix = work_dir / "layer_b_clusters"
    tmp    = work_dir / "mmseqs_tmp"
    tmp.mkdir(exist_ok=True)
    cmd = [
        "mmseqs", "easy-cluster",
        str(faa_path), str(prefix), str(tmp),
        "--min-seq-id", str(min_id),
        "-c", str(min_cov),
        "--cov-mode", "0",
        "-s", "7.5",
        "--threads", "4",
        "-v", "1",
    ]
    print(f"[INFO] Running: {' '.join(cmd)}")
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stderr[-2000:], file=sys.stderr)
        sys.exit("[ERROR] mmseqs easy-cluster failed")
    tsv = Path(str(prefix) + "_cluster.tsv")
    if not tsv.exists():
        sys.exit(f"[ERROR] Expected cluster file not found: {tsv}")
    print(f"[INFO] MMseqs2 done → {tsv}")
    return tsv


# ── Parse clusters ────────────────────────────────────────────────────────────

def parse_clusters(tsv_path):
    """Return dict: rep_header → [member_headers] (rep included in list)."""
    clusters = defaultdict(list)
    with open(tsv_path) as fh:
        for line in fh:
            parts = line.rstrip().split("\t")
            if len(parts) < 2:
                continue
            rep, member = parts[0], parts[1]
            clusters[rep].append(member)
    return dict(clusters)


def header_to_stem(header):
    """Extract stem from '{cat}__{stem}' header."""
    parts = header.split("__", 1)
    return parts[1] if len(parts) == 2 else header


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--hmm_dir",  default="hmm_library")
    p.add_argument("--work_dir", default=None,
                   help="Scratch dir (default: hmm_library/_dedup_work/)")
    p.add_argument("--min_id",   type=float, default=0.70)
    p.add_argument("--min_cov",  type=float, default=0.80)
    p.add_argument("--dry_run",  action="store_true")
    p.add_argument("--skip_emit", action="store_true",
                   help="Reuse existing all_consensus.faa (skip hmmemit)")
    args = p.parse_args()

    hmm_dir  = Path(args.hmm_dir)
    work_dir = Path(args.work_dir) if args.work_dir else hmm_dir / "_dedup_work"
    work_dir.mkdir(parents=True, exist_ok=True)
    faa_path = work_dir / "all_consensus.faa"

    rows, fieldnames = load_registry(hmm_dir)
    active = [r for r in rows if r.get("status", "active") == "active"]
    print(f"[INFO] Registry: {len(rows)} total, {len(active)} active")

    # ── Step 1: hmmemit ───────────────────────────────────────────────────────
    if args.skip_emit and faa_path.exists():
        print(f"[INFO] Skipping hmmemit — reusing {faa_path}")
    else:
        build_consensus_faa(active, faa_path, hmm_dir)

    # ── Step 2: MMseqs2 ───────────────────────────────────────────────────────
    cluster_tsv = run_mmseqs(faa_path, work_dir, args.min_id, args.min_cov)

    # ── Step 3: parse clusters ────────────────────────────────────────────────
    clusters = parse_clusters(cluster_tsv)
    multi = {rep: members for rep, members in clusters.items()
             if len(members) > 1}
    print(f"[INFO] {len(clusters)} clusters total; "
          f"{len(multi)} have >1 member (redundant)")

    if not multi:
        print("[INFO] No redundant pairs found — library already deduplicated at "
              f"{args.min_id:.0%} identity / {args.min_cov:.0%} coverage.")
        return

    # Build stem → row lookup
    stem_to_row = {r["stem"]: r for r in rows}

    # ── Step 4: select representatives, collect to-deprecate ─────────────────
    to_deprecate = {}   # stem → reason
    cluster_report = []

    for rep_header, members in multi.items():
        stems = [header_to_stem(h) for h in members]
        member_rows = [stem_to_row[s] for s in stems if s in stem_to_row]
        if not member_rows:
            continue

        # Sort by score descending; ties broken by nseq then stem
        member_rows.sort(
            key=lambda r: (_rep_score(r), int(r.get("nseq", 0) or 0), r["stem"]),
            reverse=True)
        keeper = member_rows[0]
        deprecated = member_rows[1:]

        cluster_report.append({
            "keep": keeper["stem"],
            "keep_cat": keeper.get("category", ""),
            "keep_src": keeper.get("source", ""),
            "keep_cutoff": keeper.get("cutoff", 0),
            "keep_nseq": keeper.get("nseq", 0),
            "deprecated": [(r["stem"], r.get("category", ""),
                            r.get("source", ""), r.get("cutoff", 0),
                            r.get("nseq", 0))
                           for r in deprecated],
        })
        for r in deprecated:
            to_deprecate[r["stem"]] = "deprecated_dedup_sequence_cluster"

    # ── Step 5: report ────────────────────────────────────────────────────────
    print(f"\n{'─'*70}")
    print(f"  Layer B results  ({args.min_id:.0%} id / {args.min_cov:.0%} cov, bidirectional)")
    print(f"  Clusters with redundancy : {len(cluster_report)}")
    print(f"  Models to deprecate      : {len(to_deprecate)}")
    print(f"{'─'*70}")

    for c in cluster_report:
        print(f"\n  KEEP  {c['keep']}  [{c['keep_cat']}]  "
              f"src={c['keep_src']}  cutoff={c['keep_cutoff']}  nseq={c['keep_nseq']}")
        for stem, cat, src, cutoff, nseq in c["deprecated"]:
            print(f"  drop  {stem}  [{cat}]  "
                  f"src={src}  cutoff={cutoff}  nseq={nseq}")

    print(f"\n{'─'*70}")

    if args.dry_run:
        print("[DRY RUN] Registry not modified.")
        return

    # ── Step 6: update registry ───────────────────────────────────────────────
    n_updated = 0
    for row in rows:
        if row["stem"] in to_deprecate:
            row["status"] = to_deprecate[row["stem"]]
            n_updated += 1

    save_registry(hmm_dir, rows, fieldnames)
    active_after = sum(1 for r in rows if r.get("status", "active") == "active")
    print(f"[INFO] Registry updated: {n_updated} models deprecated")
    print(f"[INFO] Active models: {len(active)} → {active_after}")
    print(f"[INFO] Commit message suggestion:")
    print(f"       Layer B dedup all categories: "
          f"{len(to_deprecate)} models deprecated_dedup_sequence_cluster")


if __name__ == "__main__":
    main()
