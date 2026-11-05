#!/usr/bin/env python3
"""
benchmark_fegenie.py — head-to-head benchmark of Efesto vs FeGenie.

Measures three things the paper needs:
  1. Speed / peak RAM   (wall-clock + maxRSS per tool)
  2. Accuracy           (precision / recall / F1 per category vs benchmark/truth.tsv)
  3. MtrA/MtoA cross-hits (the sharp claim: FeGenie mis-calls Mtr on Mto genomes
                           and vice-versa; Efesto resolves by calibrated cutoff)

Two comparison modes (see benchmark/README.md):

  accuracy  — feed BOTH tools the SAME ORF calls (one prodigal run). Isolates
              HMM library + operon logic from gene-calling differences.
  speed     — each tool runs its native path end-to-end from nucleotide FASTA.

The script shells out to `efesto` and `FeGenie.py`; it never imports them, so the
two tools stay in their own conda envs. Point at them with --efesto / --fegenie.

Stdlib only. Python 3.8+.
"""
import argparse, csv, glob, os, re, shutil, subprocess, sys, time
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from efesto.categories import canonical_category

# ─────────────────────────────────────────────────────────────────────────────
# timed subprocess with peak-RAM capture (GNU `time -v` or BSD `time -l`)
# ─────────────────────────────────────────────────────────────────────────────
def _detect_time_bin(explicit=None):
    """Return (path, flavor) where flavor in {'gnu','bsd',None}."""
    cand = [explicit] if explicit else ["/usr/bin/time", "gtime", "time"]
    for c in cand:
        if not c:
            continue
        path = shutil.which(c) or (c if os.path.exists(c) else None)
        if not path:
            continue
        # GNU time understands -v; BSD time errors on it.
        try:
            r = subprocess.run([path, "-v", "true"], capture_output=True, text=True)
            if "Maximum resident set size" in r.stderr:
                return path, "gnu"
        except Exception:
            pass
        try:
            r = subprocess.run([path, "-l", "true"], capture_output=True, text=True)
            if "maximum resident set size" in r.stderr:
                return path, "bsd"
        except Exception:
            pass
    return None, None


def _parse_maxrss(stderr, flavor):
    """Return peak RSS in MB, or None."""
    if flavor == "gnu":
        m = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", stderr)
        if m:
            return int(m.group(1)) / 1024.0          # KB -> MB
    elif flavor == "bsd":
        m = re.search(r"(\d+)\s+maximum resident set size", stderr)
        if m:
            return int(m.group(1)) / (1024.0 * 1024) # bytes -> MB
    return None


def timed_run(cmd, log_path, time_bin, flavor):
    """Run cmd, return dict(wall_s, maxrss_mb, rc). Stdout+stderr -> log_path."""
    wrapped = ([time_bin, "-v"] if flavor == "gnu"
               else [time_bin, "-l"] if flavor == "bsd" else []) + cmd
    t0 = time.perf_counter()
    r = subprocess.run(wrapped, capture_output=True, text=True)
    wall = time.perf_counter() - t0
    with open(log_path, "w") as fh:
        fh.write("### CMD\n" + " ".join(cmd) + "\n\n### STDOUT\n" + r.stdout +
                 "\n### STDERR\n" + r.stderr)
    return {"wall_s": wall,
            "maxrss_mb": _parse_maxrss(r.stderr, flavor),
            "rc": r.returncode}


# ─────────────────────────────────────────────────────────────────────────────
# gene calling (shared ORFs for accuracy mode)
# ─────────────────────────────────────────────────────────────────────────────
def call_orfs(fna_dir, fna_ext, out_dir, meta, threads):
    """Prodigal on every *.fna -> <out>/faa/<stem>.faa + <out>/gff/<stem>.gff."""
    faa_dir, gff_dir = out_dir / "faa", out_dir / "gff"
    faa_dir.mkdir(parents=True, exist_ok=True)
    gff_dir.mkdir(parents=True, exist_ok=True)
    prodigal = shutil.which("prodigal") or shutil.which("pyrodigal")
    if not prodigal:
        sys.exit("[ERROR] prodigal/pyrodigal not found in PATH")
    fnas = sorted(glob.glob(str(Path(fna_dir) / f"*.{fna_ext}")))
    if not fnas:
        sys.exit(f"[ERROR] no *.{fna_ext} in {fna_dir}")
    for fna in fnas:
        stem = Path(fna).stem
        faa, gff = faa_dir / f"{stem}.faa", gff_dir / f"{stem}.gff"
        if faa.exists() and gff.exists():
            continue
        cmd = [prodigal, "-i", fna, "-a", str(faa), "-f", "gff", "-o", str(gff), "-q"]
        if meta:
            cmd += ["-p", "meta"]
        print(f"[prodigal] {stem}")
        subprocess.run(cmd, check=True, capture_output=True)
    return faa_dir, gff_dir


# ─────────────────────────────────────────────────────────────────────────────
# output parsers -> list of hits: dict(genome, category, gene, orf, bitscore)
# ─────────────────────────────────────────────────────────────────────────────
def _norm_genome(name):
    """faa filename / assembly label -> truth label (strip ext + common suffixes)."""
    n = re.sub(r"\.(faa|fna|fasta|fa|gff)$", "", name.strip())
    return n


def parse_efesto(out_dir):
    """Efesto-results-long.tsv (tab): category genome contig orf gene bitscore ..."""
    f = Path(out_dir) / "Efesto-results-long.tsv"
    if not f.exists():
        sys.exit(f"[ERROR] Efesto output not found: {f}")
    hits = []
    with open(f) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            if not row.get("gene"):
                continue
            hits.append({"genome": _norm_genome(row["genome"]),
                         "category": canonical_category(row["category"].strip()),
                         "gene": row["gene"].strip(),
                         "orf": row.get("orf", "").strip(),
                         "bitscore": float(row["bitscore"] or 0)})
    return hits


def parse_fegenie(out_dir):
    """FeGenie-geneSummary*.csv — column names vary by version; autodetect."""
    cands = glob.glob(str(Path(out_dir) / "*geneSummary*.csv"))
    if not cands:
        cands = glob.glob(str(Path(out_dir) / "*.csv"))
    if not cands:
        sys.exit(f"[ERROR] no FeGenie geneSummary csv in {out_dir}")
    path = cands[0]
    with open(path) as fh:
        rows = list(csv.reader(fh))
    # find header row (contains a cell 'category' or 'process')
    hdr_i, hdr = None, None
    for i, row in enumerate(rows):
        low = [c.strip().lower() for c in row]
        if "category" in low or "process" in low:
            hdr_i, hdr = i, low
            break
    if hdr_i is None:
        sys.exit(f"[ERROR] could not find header in {path}")

    def col(*names):
        for nm in names:
            if nm in hdr:
                return hdr.index(nm)
        return None
    ci = col("category", "process")
    gi = col("genome/assembly", "genome", "assembly", "bin")
    oi = col("orf")
    ni = col("gene", "hmm")
    bi = col("bitscore", "bit_score")
    hits = []
    for row in rows[hdr_i + 1:]:
        if not row or all(c.strip() in ("", "#") for c in row):
            continue
        if row[ci].strip() == "#":
            continue
        try:
            gene = row[ni].strip() if ni is not None else ""
            if not gene or gene == "#":
                continue
            bs = 0.0
            if bi is not None and row[bi].strip() not in ("", "#"):
                try: bs = float(row[bi])
                except ValueError: bs = 0.0
            hits.append({"genome": _norm_genome(row[gi]) if gi is not None else "",
                         "category": canonical_category(row[ci].strip()),
                         "gene": gene,
                         "orf": row[oi].strip() if oi is not None else "",
                         "bitscore": bs})
        except IndexError:
            continue
    return hits


# ─────────────────────────────────────────────────────────────────────────────
# truth + scoring
# ─────────────────────────────────────────────────────────────────────────────
def load_truth(path):
    """-> list of dict(genome, category, gene, label, note)."""
    truth = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if parts[0].strip() == "genome":       # header
                continue
            if len(parts) < 4:
                continue
            truth.append({"genome": parts[0].strip(), "category": canonical_category(parts[1].strip()),
                          "gene": parts[2].strip(), "label": parts[3].strip(),
                          "note": parts[4].strip() if len(parts) > 4 else ""})
    return truth


def _matches(hits, genome, category, gene):
    """Does the tool report >=1 hit for (genome, category, gene)?
       gene '*' = any gene in category. Otherwise case-insensitive substring."""
    g = gene.lower()
    for h in hits:
        if h["genome"] != genome or h["category"] != category:
            continue
        if gene == "*" or g in h["gene"].lower():
            return True
    return False


def score(hits, truth):
    """Per-tool confusion vs truth. Only rows labelled present/absent/diagnostic score.
       present/diagnostic-present -> recall (want match).
       absent/diagnostic-absent   -> precision (want NO match)."""
    genomes = sorted({t["genome"] for t in truth})
    per_cat = defaultdict(lambda: {"TP": 0, "FN": 0, "FP": 0, "TN": 0})
    detail = []            # per truth row: outcome
    crosshit = []          # MtrA/MtoA specific
    for t in truth:
        got = _matches(hits, t["genome"], t["category"], t["gene"])
        want_present = t["label"] in ("present", "diagnostic")
        # 'diagnostic' can be present-or-absent; decide by note keyword 'absent'/'must be absent'
        if t["label"] == "diagnostic" and re.search(r"\babsent\b|must be absent|NO ", t["note"]):
            want_present = False
        cat = t["category"]
        if want_present:
            key = "TP" if got else "FN"
        else:
            key = "FP" if got else "TN"
        per_cat[cat][key] += 1
        detail.append({**t, "reported": got, "want_present": want_present, "outcome": key})
        # cross-hit tracking: Mtr on an oxidizer, Mto on a reducer
        gl = t["gene"].lower()
        if gl in ("mtra", "mtoa") and not want_present and got:
            crosshit.append({"genome": t["genome"], "gene": t["gene"],
                             "category": cat, "note": t["note"]})
    return per_cat, detail, crosshit


def prf(c):
    tp, fp, fn = c["TP"], c["FP"], c["FN"]
    p = tp / (tp + fp) if (tp + fp) else float("nan")
    r = tp / (tp + fn) if (tp + fn) else float("nan")
    f = 2 * p * r / (p + r) if (p and r and p == p and r == r) else float("nan")
    return p, r, f


# ─────────────────────────────────────────────────────────────────────────────
# writers
# ─────────────────────────────────────────────────────────────────────────────
def _fmt(x):
    return "NA" if x != x else f"{x:.3f}"   # x!=x  -> NaN


def write_accuracy(out, tool, per_cat, detail, crosshit):
    with open(out / f"accuracy-{tool}.tsv", "w") as fh:
        fh.write("category\tTP\tFP\tFN\tTN\tprecision\trecall\tf1\n")
        agg = {"TP": 0, "FP": 0, "FN": 0, "TN": 0}
        for cat in sorted(per_cat):
            c = per_cat[cat]
            for k in agg: agg[k] += c[k]
            p, r, f = prf(c)
            fh.write(f"{cat}\t{c['TP']}\t{c['FP']}\t{c['FN']}\t{c['TN']}\t"
                     f"{_fmt(p)}\t{_fmt(r)}\t{_fmt(f)}\n")
        p, r, f = prf(agg)
        fh.write(f"OVERALL\t{agg['TP']}\t{agg['FP']}\t{agg['FN']}\t{agg['TN']}\t"
                 f"{_fmt(p)}\t{_fmt(r)}\t{_fmt(f)}\n")
    with open(out / f"detail-{tool}.tsv", "w") as fh:
        fh.write("genome\tcategory\tgene\tlabel\twant_present\treported\toutcome\tnote\n")
        for d in detail:
            fh.write(f"{d['genome']}\t{d['category']}\t{d['gene']}\t{d['label']}\t"
                     f"{d['want_present']}\t{d['reported']}\t{d['outcome']}\t{d['note']}\n")
    with open(out / f"crosshits-{tool}.tsv", "w") as fh:
        fh.write("genome\tgene\tcategory\tnote\n")
        for c in crosshit:
            fh.write(f"{c['genome']}\t{c['gene']}\t{c['category']}\t{c['note']}\n")
    return agg


def write_speed(out, rows):
    with open(out / "speed.tsv", "w") as fh:
        fh.write("tool\tn_genomes\twall_s\tmaxrss_mb\trc\n")
        for r in rows:
            ram = "" if r["maxrss_mb"] is None else f"{r['maxrss_mb']:.0f}"
            fh.write(f"{r['tool']}\t{r['n']}\t{r['wall_s']:.1f}\t{ram}\t{r['rc']}\n")


# ─────────────────────────────────────────────────────────────────────────────
# main
# ─────────────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(description="Benchmark Efesto vs FeGenie")
    ap.add_argument("--mode", choices=["accuracy", "speed", "both"], default="accuracy",
                    help="accuracy=shared ORFs; speed=native end-to-end from FNA")
    ap.add_argument("--fna_dir", help="nucleotide assemblies (*.fna) — required for speed / to call ORFs")
    ap.add_argument("--faa_dir", help="pre-called ORFs (*.faa); skips prodigal in accuracy mode")
    ap.add_argument("--gff_dir", help="matching GFFs (Efesto bp-clustering); optional")
    ap.add_argument("--fna_ext", default="fna")
    ap.add_argument("--hmm_dir", required=True, help="Efesto hmm_library/")
    ap.add_argument("--truth", default="benchmark/truth.tsv")
    ap.add_argument("--out", default="benchmark/results")
    ap.add_argument("--efesto", default="efesto", help="efesto executable (its conda env)")
    ap.add_argument("--fegenie", default="FeGenie.py", help="FeGenie.py (its conda env)")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--meta", action="store_true", help="prodigal -p meta (metagenome contigs)")
    ap.add_argument("--time_bin", help="path to GNU/BSD time (autodetected)")
    ap.add_argument("--skip_fegenie", action="store_true", help="run/score Efesto only")
    ap.add_argument("--compare_only", action="store_true",
                    help="skip tool runs; re-parse existing <out>/efesto_acc + "
                         "fegenie_acc and re-score vs --truth (fast truth/parser iteration)")
    args = ap.parse_args()

    if args.compare_only and args.mode == "speed":
        sys.exit("[ERROR] --compare_only only applies to accuracy scoring, not speed mode")

    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    time_bin, flavor = _detect_time_bin(args.time_bin)
    if not time_bin:
        print("[WARN] no GNU/BSD time found — RAM not measured, wall-clock only")

    speed_rows = []

    # ---- accuracy mode: shared ORFs ----
    if args.mode in ("accuracy", "both"):
        ef_out, fg_out = out / "efesto_acc", out / "fegenie_acc"

        if args.compare_only:
            # re-score from existing tool outputs — no runs, no gene calling
            if not ef_out.exists():
                sys.exit(f"[ERROR] --compare_only: {ef_out} not found — run without it first")
            print("[compare-only] re-parsing existing outputs, skipping tool runs")
            ef_hits = parse_efesto(ef_out)
            fg_hits = parse_fegenie(fg_out) if (not args.skip_fegenie and fg_out.exists()) else []
        else:
            if args.faa_dir:
                faa_dir, gff_dir = Path(args.faa_dir), (Path(args.gff_dir) if args.gff_dir else None)
            elif args.fna_dir:
                faa_dir, gff_dir = call_orfs(args.fna_dir, args.fna_ext, out / "orfs",
                                             args.meta, args.threads)
            else:
                sys.exit("[ERROR] accuracy mode needs --faa_dir or --fna_dir")
            n = len(glob.glob(str(faa_dir / "*.faa")))

            # Efesto
            cmd = [args.efesto, "--faa_dir", str(faa_dir), "--hmm_dir", args.hmm_dir,
                   "--out", str(ef_out), "--threads", str(args.threads)]
            if gff_dir:
                cmd += ["--gff_dir", str(gff_dir)]
            print("[run] Efesto (accuracy)…")
            m = timed_run(cmd, out / "log-efesto-acc.txt", time_bin, flavor)
            speed_rows.append({"tool": "efesto", "n": n, **m})
            if m["rc"] != 0:
                print(f"[WARN] efesto rc={m['rc']} — see {out}/log-efesto-acc.txt")
            ef_hits = parse_efesto(ef_out)

            # FeGenie
            fg_hits = []
            if not args.skip_fegenie:
                cmd = [args.fegenie, "-bin_dir", str(faa_dir), "-bin_ext", "faa",
                       "-out", str(fg_out), "-t", str(args.threads), "--orfs"]
                print("[run] FeGenie (accuracy)…")
                m = timed_run(cmd, out / "log-fegenie-acc.txt", time_bin, flavor)
                speed_rows.append({"tool": "fegenie", "n": n, **m})
                if m["rc"] != 0:
                    print(f"[WARN] fegenie rc={m['rc']} — see {out}/log-fegenie-acc.txt")
                fg_hits = parse_fegenie(fg_out)

        truth = load_truth(args.truth)
        pc, det, cross = score(ef_hits, truth)
        ef_agg = write_accuracy(out, "efesto", pc, det, cross)
        print(f"[efesto ] TP={ef_agg['TP']} FP={ef_agg['FP']} FN={ef_agg['FN']} "
              f"cross-hits={len(cross)}")
        if fg_hits:
            pc, det, cross = score(fg_hits, truth)
            fg_agg = write_accuracy(out, "fegenie", pc, det, cross)
            print(f"[fegenie] TP={fg_agg['TP']} FP={fg_agg['FP']} FN={fg_agg['FN']} "
                  f"cross-hits={len(cross)}")

    # ---- speed mode: native end-to-end ----
    if args.mode in ("speed", "both") and not args.compare_only:
        if not args.fna_dir:
            sys.exit("[ERROR] speed mode needs --fna_dir")
        n = len(glob.glob(str(Path(args.fna_dir) / f"*.{args.fna_ext}")))
        ef_out = out / "efesto_speed"
        cmd = [args.efesto, "--fna_dir", args.fna_dir, "--fna_ext", args.fna_ext,
               "--hmm_dir", args.hmm_dir, "--out", str(ef_out),
               "--threads", str(args.threads)]
        if args.meta:
            cmd.append("--meta")
        print("[run] Efesto (speed, native)…")
        m = timed_run(cmd, out / "log-efesto-speed.txt", time_bin, flavor)
        speed_rows.append({"tool": "efesto-native", "n": n, **m})
        if not args.skip_fegenie:
            fg_out = out / "fegenie_speed"
            cmd = [args.fegenie, "-bin_dir", args.fna_dir, "-bin_ext", args.fna_ext,
                   "-out", str(fg_out), "-t", str(args.threads)]
            if args.meta:
                cmd.append("-meta")
            print("[run] FeGenie (speed, native)…")
            m = timed_run(cmd, out / "log-fegenie-speed.txt", time_bin, flavor)
            speed_rows.append({"tool": "fegenie-native", "n": n, **m})

    if speed_rows:
        write_speed(out, speed_rows)
        print(f"\n[done] tables in {out}/  (accuracy-*.tsv, crosshits-*.tsv, speed.tsv)")


if __name__ == "__main__":
    main()
