#!/usr/bin/env python3
"""
build_mtr_mto_subfamily_hmms.py
================================
Rebuild MtrA (iron reduction) and MtoA (iron oxidation) HMM profiles from
curated seed sequences.

Problem: the current FeGenie-derived MtrA and MtoA models cross-hit each other
at ~2× the calibrated cutoff (265-289 bits vs cutoff 140). Both are built from
only 9-10 sequences and cannot discriminate the two subfamilies.

Cross-hit evidence:
  MtrA HMM vs MtoA consensus: 265.8 bits (cutoff 140)
  MtoA HMM vs MtrA consensus: 289.8 bits (cutoff 140)

This script:
  1. Reads curated seed accessions from data/seeds/mtr_mto_seeds.tsv
  2. Downloads sequences via NCBI Entrez (requires Biopython + NCBI_EMAIL env var)
  3. Aligns each subfamily with MAFFT
  4. Builds HMMs with hmmbuild
  5. Runs cross-validation to confirm discrimination gap
  6. Writes new HMMs to a staging directory for review before deploying

Usage:
    export NCBI_EMAIL=you@example.com
    python scripts/build_mtr_mto_subfamily_hmms.py [--seed_tsv PATH] [--out_dir PATH]

Options:
    --seed_tsv   Path to seed TSV (default: data/seeds/mtr_mto_seeds.tsv)
    --out_dir    Directory for output HMMs (default: hmm_library/_staging/mtr_mto/)
    --skip_dl    Reuse sequences already in out_dir (skip Entrez download)
    --dry_run    Validate inputs and show what would run; do not execute
"""

import argparse
import csv
import os
import subprocess
import sys
from pathlib import Path

try:
    from Bio import Entrez, SeqIO
except ImportError:
    sys.exit("[ERROR] Biopython required: conda install -c conda-forge biopython")


# ── Seed TSV schema ───────────────────────────────────────────────────────────
# accession  subfamily  organism  role  reference_doi
# Fields:
#   accession  : NCBI protein accession (e.g. WP_011071012.1)
#   subfamily  : 'MtrA' or 'MtoA'
#   organism   : organism name for provenance
#   role       : 'iron_reduction' or 'iron_oxidation'
#   reference  : DOI of paper confirming function

SEED_TSV_HEADER = ["accession", "subfamily", "organism", "role", "reference"]


def load_seeds(tsv_path):
    seeds = {"MtrA": [], "MtoA": []}
    with open(tsv_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sub = row["subfamily"]
            if sub not in seeds:
                print(f"[WARN] Unknown subfamily '{sub}' in seed TSV — skipping",
                      file=sys.stderr)
                continue
            seeds[sub].append(row)
    return seeds


def fetch_sequences(accessions, out_faa, email, db="protein"):
    Entrez.email = email
    print(f"[INFO] Fetching {len(accessions)} sequences from NCBI...")
    with Entrez.efetch(db=db, id=",".join(accessions),
                       rettype="fasta", retmode="text") as handle:
        seqs = list(SeqIO.parse(handle, "fasta"))
    if not seqs:
        sys.exit("[ERROR] No sequences returned from Entrez")
    if len(seqs) != len(accessions):
        print(f"[WARN] Requested {len(accessions)}, received {len(seqs)}",
              file=sys.stderr)
    SeqIO.write(seqs, out_faa, "fasta")
    print(f"[INFO] Wrote {len(seqs)} sequences → {out_faa}")
    return seqs


def run_mafft(in_faa, out_aln):
    cmd = ["mafft", "--auto", "--quiet", str(in_faa)]
    with open(out_aln, "w") as fh:
        r = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE)
    if r.returncode != 0:
        sys.exit(f"[ERROR] mafft failed:\n{r.stderr.decode()}")
    print(f"[INFO] MAFFT alignment → {out_aln}")


def run_hmmbuild(aln, hmm_out, name):
    cmd = ["hmmbuild", "-n", name, str(hmm_out), str(aln)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"[ERROR] hmmbuild failed:\n{r.stderr}")
    print(f"[INFO] hmmbuild → {hmm_out}")


def run_hmmcalibrate(hmm_path):
    cmd = ["hmmcalibrate", str(hmm_path)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"[WARN] hmmcalibrate failed (E-value stats may be absent): "
              f"{r.stderr.strip()}", file=sys.stderr)


def _parse_tblout(tbl_path):
    """Parse hmmsearch --tblout file; return list of (seq_id, bitscore)."""
    hits = []
    with open(tbl_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            try:
                hits.append((parts[0], float(parts[5])))
            except ValueError:
                pass
    return hits


def cross_validate(mtrA_hmm, mtoA_hmm, mtrA_seqs, mtoA_seqs):
    """Run each HMM against both sequence sets; print bitscore summary."""
    import tempfile
    print(f"\n{'─'*65}")
    print("  Cross-validation: discrimination gap")
    print(f"{'─'*65}")
    combos = [
        ("MtrA HMM", mtrA_hmm, "MtrA seqs (correct)", mtrA_seqs),
        ("MtrA HMM", mtrA_hmm, "MtoA seqs (cross-hit)", mtoA_seqs),
        ("MtoA HMM", mtoA_hmm, "MtrA seqs (cross-hit)", mtrA_seqs),
        ("MtoA HMM", mtoA_hmm, "MtoA seqs (correct)", mtoA_seqs),
    ]
    with tempfile.NamedTemporaryFile(suffix=".tbl", mode="w", delete=False) as tf:
        tbl_path = tf.name
    for hmm_label, hmm_path, seq_label, seq_path in combos:
        subprocess.run(
            ["hmmsearch", "--tblout", tbl_path, "--noali", "-E", "1e-3",
             str(hmm_path), str(seq_path)],
            capture_output=True, text=True)
        hits = _parse_tblout(tbl_path)
        if hits:
            scores = [s for _, s in hits]
            mn, mx, avg = min(scores), max(scores), sum(scores) / len(scores)
            print(f"  {hmm_label} vs {seq_label}:")
            print(f"    n={len(scores)}  min={mn:.1f}  max={mx:.1f}  mean={avg:.1f}")
            for seq_id, sc in sorted(hits, key=lambda x: -x[1]):
                print(f"      {sc:7.1f}  {seq_id}")
        else:
            print(f"  {hmm_label} vs {seq_label}: no hits")
    print(f"{'─'*65}")
    print("\n  Discrimination is adequate if correct-hit scores are")
    print("  consistently higher than cross-hit scores with a clear gap.")
    print("  If the gap is < 50 bits, set cutoff at correct-hit min − 10.")
    print("  If cross-hits exceed correct-hit min, models need more seeds.")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--seed_tsv", default="data/seeds/mtr_mto_seeds.tsv")
    p.add_argument("--out_dir",  default="hmm_library/_staging/mtr_mto")
    p.add_argument("--skip_dl",  action="store_true",
                   help="Reuse existing FASTA files in out_dir")
    p.add_argument("--dry_run",  action="store_true")
    args = p.parse_args()

    seed_path = Path(args.seed_tsv)
    out_dir   = Path(args.out_dir)

    if not seed_path.exists():
        sys.exit(f"[ERROR] Seed TSV not found: {seed_path}\n"
                 f"Create it with columns: {', '.join(SEED_TSV_HEADER)}\n"
                 f"See data/seeds/mtr_mto_seeds.tsv.example for a template.")

    email = os.environ.get("NCBI_EMAIL", "")
    if not email and not args.skip_dl:
        sys.exit("[ERROR] Set NCBI_EMAIL environment variable before running.")

    seeds = load_seeds(seed_path)
    for sub, rows in seeds.items():
        print(f"[INFO] {sub}: {len(rows)} seeds")
        if len(rows) < 5:
            print(f"[WARN] {sub} has only {len(rows)} seeds — "
                  f"recommend ≥ 10 for reliable HMM boundaries", file=sys.stderr)

    if args.dry_run:
        print("[DRY RUN] Would download, align, and build HMMs for:")
        for sub, rows in seeds.items():
            accs = [r["accession"] for r in rows]
            print(f"  {sub}: {', '.join(accs)}")
        return

    out_dir.mkdir(parents=True, exist_ok=True)

    # ── Download or reuse ─────────────────────────────────────────────────────
    faa = {}
    for sub, rows in seeds.items():
        faa_path = out_dir / f"{sub}_seeds.faa"
        faa[sub]  = faa_path
        if args.skip_dl and faa_path.exists():
            print(f"[INFO] Reusing {faa_path}")
        else:
            accs = [r["accession"] for r in rows]
            fetch_sequences(accs, faa_path, email)

    # ── Align ─────────────────────────────────────────────────────────────────
    aln = {}
    for sub in seeds:
        aln_path = out_dir / f"{sub}_seeds.aln"
        aln[sub]  = aln_path
        run_mafft(faa[sub], aln_path)

    # ── Build HMMs ────────────────────────────────────────────────────────────
    hmm = {}
    for sub in seeds:
        hmm_path = out_dir / f"{sub}.hmm"
        hmm[sub]  = hmm_path
        role = "iron_reduction" if sub == "MtrA" else "iron_oxidation"
        run_hmmbuild(aln[sub], hmm_path, name=sub)
        print(f"[INFO] HMM category: {role} → {hmm_path}")

    # ── Cross-validate ────────────────────────────────────────────────────────
    cross_validate(hmm["MtrA"], hmm["MtoA"], faa["MtrA"], faa["MtoA"])

    print(f"\n[INFO] Staging HMMs in {out_dir}/")
    print("[INFO] Review cross-validation scores above before deploying.")
    print("[INFO] To deploy:")
    print(f"       cp {out_dir}/MtrA.hmm hmm_library/iron_reduction/MtrA.hmm")
    print(f"       cp {out_dir}/MtoA.hmm hmm_library/iron_oxidation/MtoA.hmm")
    print("[INFO] Then update hmm_registry.tsv: nseq, cutoff, source='curated'")
    print("[INFO] Suggested commit message:")
    print("       Rebuild MtrA/MtoA subfamily HMMs: curated seeds, cross-validated")


if __name__ == "__main__":
    main()
