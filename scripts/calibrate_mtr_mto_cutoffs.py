#!/usr/bin/env python3
"""
calibrate_mtr_mto_cutoffs.py
=============================
Calibrate TC (Trusted Cutoff), NC (Noise Cutoff), and GA (Gathering) for the
curated MtrA and MtoA HMM profiles.

Methodology:
  1. Download all "DmsE family decaheme c-type cytochrome" sequences from NCBI
     (~3800 sequences = the full TIGR03508 superfamily, our calibration universe).
  2. Search both HMMs against this full family.
  3. Classify each hit by taxonomy into one of:
       MtrA_confirmed  — title contains "MtrA" (NCBI-annotated iron-reduction)
       MtoA_candidate  — organism is Gallionellaceae or Sideroxydans
       Shewanella_other — Shewanella sp. annotated as DmsE, not MtrA
       other           — all other organisms (DmsE, uncharacterized)
  4. Report score distributions per class; suggest TC and NC values.
  5. IMPORTANT: TC = lowest confirmed TRUE POSITIVE. NC = highest confirmed
     non-member. If NC >= TC, the model cannot discriminate — rebuild needed.

Usage:
    export NCBI_EMAIL=you@example.com
    python scripts/calibrate_mtr_mto_cutoffs.py [--skip_dl] [--out_dir PATH]

Outputs:
    <out_dir>/calibration_universe.faa     — all downloaded sequences
    <out_dir>/mtrA_vs_universe.tbl         — hmmsearch tblout
    <out_dir>/mtoA_vs_universe.tbl         — hmmsearch tblout
    <out_dir>/calibration_report.tsv       — per-sequence scores + classification
    <out_dir>/calibration_summary.txt      — TC/NC/GA recommendations
"""

import argparse
import csv
import os
import subprocess
import sys
import time
from pathlib import Path

try:
    from Bio import Entrez, SeqIO
except ImportError:
    sys.exit("[ERROR] Biopython required: conda install -c conda-forge biopython")

BATCH_SIZE = 200
SEARCH_TERM = 'decaheme[Title] AND "DmsE"[Title]'

# Taxonomy-based classification rules (checked in order, first match wins)
GALLIONELLACEAE = {"Gallionella", "Sideroxydans", "Ferriphaselus",
                   "Litorivivens", "Gallionellaceae"}
MTRA_TITLE_TOKENS = {"MtrA", "mtrA", "metal-reducing", "metal reducing"}
MTOA_TITLE_TOKENS = {"MtoA", "mtoA", "metal-oxidizing", "metal oxidizing"}


def classify(seq_id, description, organism=""):
    """Return class label for a sequence based on description + organism."""
    desc_upper = description.upper()
    org_words = set(organism.split())

    if any(t.upper() in desc_upper for t in MTRA_TITLE_TOKENS):
        return "MtrA_confirmed"
    if any(t.upper() in desc_upper for t in MTOA_TITLE_TOKENS):
        return "MtoA_confirmed"
    if any(g in organism for g in GALLIONELLACEAE):
        return "MtoA_candidate"
    if "Shewanella" in organism:
        return "Shewanella_DmsE"
    return "other"


# ── Download ──────────────────────────────────────────────────────────────────

def fetch_all(out_faa, email):
    Entrez.email = email
    h = Entrez.esearch(db="protein", term=SEARCH_TERM, retmax=1, usehistory="y")
    rec = Entrez.read(h)
    total = int(rec["Count"])
    webenv = rec["WebEnv"]
    query_key = rec["QueryKey"]
    print(f"[INFO] {total} sequences to download in batches of {BATCH_SIZE}")

    written = 0
    with open(out_faa, "w") as fh:
        for start in range(0, total, BATCH_SIZE):
            attempt = 0
            while attempt < 3:
                try:
                    h2 = Entrez.efetch(
                        db="protein", rettype="fasta", retmode="text",
                        retstart=start, retmax=BATCH_SIZE,
                        webenv=webenv, query_key=query_key)
                    batch = h2.read()
                    fh.write(batch)
                    written += batch.count(">")
                    print(f"\r[INFO] Downloaded {written}/{total}", end="", flush=True)
                    time.sleep(0.35)  # NCBI rate limit: 3 req/s
                    break
                except Exception as e:
                    attempt += 1
                    print(f"\n[WARN] Batch {start}: {e} — retry {attempt}/3",
                          file=sys.stderr)
                    time.sleep(2)
    print()
    print(f"[INFO] Wrote {written} sequences → {out_faa}")
    return written


# ── hmmsearch ────────────────────────────────────────────────────────────────

def run_hmmsearch(hmm_path, faa_path, tbl_path, threads=4):
    cmd = ["hmmsearch", "--tblout", str(tbl_path), "--noali",
           "--cpu", str(threads), "-E", "0.01",
           str(hmm_path), str(faa_path)]
    print(f"[INFO] hmmsearch: {hmm_path.name} vs universe...")
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stderr[-1000:], file=sys.stderr)
        sys.exit(f"[ERROR] hmmsearch failed for {hmm_path}")
    hits = sum(1 for l in open(tbl_path) if not l.startswith("#") and l.strip())
    print(f"[INFO] {hits} hits → {tbl_path}")


# ── Parse tblout ──────────────────────────────────────────────────────────────

def parse_tblout(tbl_path):
    hits = {}
    with open(tbl_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split(None, 18)
            if len(parts) < 6:
                continue
            seq_id = parts[0]
            try:
                score = float(parts[5])
                evalue = float(parts[4])
            except ValueError:
                continue
            desc = parts[18].strip() if len(parts) > 18 else ""
            hits[seq_id] = {"score": score, "evalue": evalue, "desc": desc}
    return hits


# ── Classify all hits ─────────────────────────────────────────────────────────

def build_report(mtrA_hits, mtoA_hits, all_seq_ids_desc):
    rows = []
    all_ids = set(mtrA_hits) | set(mtoA_hits)
    for seq_id in all_ids:
        desc = all_seq_ids_desc.get(seq_id, "")
        # Extract organism from description: "acc desc [organism]"
        org = ""
        if "[" in desc and desc.endswith("]"):
            org = desc[desc.rfind("[") + 1:-1]
        cls = classify(seq_id, desc, org)
        mtrA_score = mtrA_hits.get(seq_id, {}).get("score", 0.0)
        mtoA_score = mtoA_hits.get(seq_id, {}).get("score", 0.0)
        rows.append({
            "seq_id":     seq_id,
            "class":      cls,
            "organism":   org,
            "mtrA_score": mtrA_score,
            "mtoA_score": mtoA_score,
            "desc":       desc[:120],
        })
    rows.sort(key=lambda r: -r["mtrA_score"])
    return rows


def print_summary(rows, out_txt):
    by_class = {}
    for r in rows:
        by_class.setdefault(r["class"], []).append(r)

    lines = []
    lines.append("=" * 70)
    lines.append("  MtrA / MtoA calibration summary")
    lines.append("=" * 70)

    for model in ("mtrA", "mtoA"):
        score_key = f"{model}_score"
        lines.append(f"\n── {model.upper()} HMM ──────────────────────────────────────────────")
        for cls in ("MtrA_confirmed", "MtoA_confirmed", "MtoA_candidate",
                    "Shewanella_DmsE", "other"):
            seqs = by_class.get(cls, [])
            scores = sorted([r[score_key] for r in seqs if r[score_key] > 0], reverse=True)
            if not scores:
                lines.append(f"  {cls:20s}  n=0")
                continue
            mn, mx = min(scores), max(scores)
            mean = sum(scores) / len(scores)
            lines.append(f"  {cls:20s}  n={len(scores):4d}  "
                         f"min={mn:7.1f}  max={mx:7.1f}  mean={mean:7.1f}")
            # Show top 5 and bottom 5
            if len(scores) > 5:
                for sc in scores[:3]:
                    lines.append(f"    top   {sc:7.1f}")
                lines.append("    ...")
                for sc in scores[-3:]:
                    lines.append(f"    bot   {sc:7.1f}")
            else:
                for sc in scores:
                    lines.append(f"          {sc:7.1f}")

    # TC / NC recommendations
    lines.append("\n── RECOMMENDED CUTOFFS ─────────────────────────────────────────────")
    for model in ("mtrA", "mtoA"):
        score_key = f"{model}_score"
        positive_cls = "MtrA_confirmed" if model == "mtrA" else "MtoA_candidate"
        positive_seqs = by_class.get(positive_cls, [])
        if model == "mtoA":
            positive_seqs += by_class.get("MtoA_confirmed", [])

        positive_scores = [r[score_key] for r in positive_seqs if r[score_key] > 0]
        # Noise = everything NOT in the positive class
        noise_scores = [r[score_key] for r in rows
                        if r["class"] != positive_cls
                        and (model == "mtrA" or r["class"] != "MtoA_confirmed")
                        and r[score_key] > 0]

        if positive_scores:
            tc = min(positive_scores)
            lines.append(f"\n  {model.upper()} TC (lowest confirmed true positive): {tc:.1f}")
        if noise_scores:
            nc = max(noise_scores)
            lines.append(f"  {model.upper()} NC (highest non-member):             {nc:.1f}")
        if positive_scores and noise_scores:
            ga = min(positive_scores) - 10
            gap = min(positive_scores) - max(noise_scores)
            lines.append(f"  {model.upper()} GA (TC - 10, recommended):           {ga:.1f}")
            lines.append(f"  {model.upper()} discrimination gap (TC - NC):        {gap:.1f} bits")
            if gap < 0:
                lines.append(f"  *** WARNING: NC > TC — model cannot discriminate! "
                              f"Needs more/better seeds. ***")
            elif gap < 50:
                lines.append(f"  *** WARNING: gap < 50 bits — twilight zone is wide. "
                              f"Manual inspection of borderline sequences required. ***")

    lines.append("\n" + "=" * 70)
    text = "\n".join(lines)
    print(text)
    with open(out_txt, "w") as fh:
        fh.write(text + "\n")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out_dir",  default="hmm_library/_calibration/mtr_mto")
    p.add_argument("--mtrA_hmm", default="hmm_library/iron_reduction/MtrA.hmm")
    p.add_argument("--mtoA_hmm", default="hmm_library/iron_oxidation/MtoA.hmm")
    p.add_argument("--skip_dl",  action="store_true",
                   help="Reuse existing calibration_universe.faa")
    p.add_argument("--threads",  type=int, default=4)
    args = p.parse_args()

    email = os.environ.get("NCBI_EMAIL", "")
    if not email and not args.skip_dl:
        sys.exit("[ERROR] Set NCBI_EMAIL environment variable")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    mtrA_hmm = Path(args.mtrA_hmm)
    mtoA_hmm = Path(args.mtoA_hmm)
    faa_path  = out_dir / "calibration_universe.faa"
    mtrA_tbl  = out_dir / "mtrA_vs_universe.tbl"
    mtoA_tbl  = out_dir / "mtoA_vs_universe.tbl"
    report    = out_dir / "calibration_report.tsv"
    summary   = out_dir / "calibration_summary.txt"

    # Step 1 — download
    if args.skip_dl and faa_path.exists():
        n = sum(1 for l in open(faa_path) if l.startswith(">"))
        print(f"[INFO] Reusing {faa_path} ({n} sequences)")
    else:
        fetch_all(faa_path, email)

    # Build seq_id → description map
    print("[INFO] Indexing sequences...")
    id_to_desc = {}
    for rec in SeqIO.parse(faa_path, "fasta"):
        id_to_desc[rec.id] = rec.description

    # Step 2 — hmmsearch
    run_hmmsearch(mtrA_hmm, faa_path, mtrA_tbl, args.threads)
    run_hmmsearch(mtoA_hmm, faa_path, mtoA_tbl, args.threads)

    # Step 3 — parse and classify
    mtrA_hits = parse_tblout(mtrA_tbl)
    mtoA_hits = parse_tblout(mtoA_tbl)
    print(f"[INFO] MtrA hits: {len(mtrA_hits)} | MtoA hits: {len(mtoA_hits)}")

    rows = build_report(mtrA_hits, mtoA_hits, id_to_desc)

    # Write TSV report
    fields = ["seq_id", "class", "organism", "mtrA_score", "mtoA_score", "desc"]
    with open(report, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"[INFO] Full report → {report}")

    # Step 4 — print summary with TC/NC/GA
    print_summary(rows, summary)
    print(f"[INFO] Summary saved → {summary}")


if __name__ == "__main__":
    main()
