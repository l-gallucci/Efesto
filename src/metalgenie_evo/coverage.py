"""BAM/depth-file coverage loading and TPM normalisation."""

import os
import subprocess
import sys
from collections import defaultdict


def compute_coverage_from_bam(bam_path, out_dir, label):
    result = {}
    df = out_dir / f"{label}.coverage"
    r = subprocess.run(
        ["samtools", "coverage", "-H", "-o", str(df), bam_path],
        capture_output=True, text=True)
    if r.returncode == 0 and df.exists():
        with open(df) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                p = line.rstrip().split("\t")
                if len(p) >= 7:
                    try:
                        result[p[0]] = {
                            "depth":  float(p[6]),
                            "reads":  int(p[3]),
                            "length": int(p[2]) - int(p[1]) + 1,
                        }
                    except (ValueError, IndexError):
                        pass
        return result

    print(f"  [WARN] samtools coverage failed for {label}, trying depth",
          file=sys.stderr)
    r = subprocess.run(["samtools", "depth", "-a", bam_path],
                       capture_output=True, text=True)
    if r.returncode != 0:
        print(f"  [ERROR] samtools depth failed: {r.stderr[:80]}", file=sys.stderr)
        return {}
    sums   = defaultdict(float)
    counts = defaultdict(int)
    for line in r.stdout.splitlines():
        p = line.split("\t")
        if len(p) >= 3:
            try:
                sums[p[0]]   += float(p[2])
                counts[p[0]] += 1
            except ValueError:
                pass
    for c in sums:
        n = counts[c]
        result[c] = {"depth": sums[c] / n if n else 0.0, "reads": 0, "length": n}
    return result


def load_depth_file(path):
    result = {}
    if not os.path.isfile(path):
        print(f"  [WARN] Depth file not found: {path}", file=sys.stderr)
        return result
    with open(path) as fh:
        lines = [l.rstrip() for l in fh if l.strip()]
    if not lines:
        return result
    cols = lines[0].split("\t")

    def _row(d=0.0, r=0, l=0):
        return {"depth": float(d), "reads": int(r), "length": int(l)}

    if cols[0].lower() in ("contigname", "#contigname"):
        for line in lines[1:]:
            p = line.split("\t")
            if len(p) >= 3 and p[0] != "contigName":
                try:
                    result[p[0]] = _row(p[2], 0, int(p[1]))
                except (ValueError, IndexError):
                    pass
        return result

    if cols[0].lower() == "rname" or (len(cols) >= 7 and cols[6].lower() == "meandepth"):
        for line in lines[1:]:
            p = line.split("\t")
            if len(p) >= 7:
                try:
                    result[p[0]] = _row(p[6], p[3], int(p[2]) - int(p[1]) + 1)
                except (ValueError, IndexError):
                    pass
        return result

    if cols[0].startswith("#"):
        for line in lines[1:]:
            if line.startswith("#"):
                continue
            p = line.split("\t")
            if len(p) >= 3:
                try:
                    result[p[0]] = _row(p[1], 0, int(p[2]))
                except (ValueError, IndexError):
                    pass
        return result

    for line in lines:
        p = line.split("\t")
        if len(p) >= 2:
            try:
                result[p[0]] = _row(p[1])
            except ValueError:
                pass
    return result


def load_bams_tsv(path):
    d = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line or line.startswith("#"):
                continue
            p = line.split("\t")
            if len(p) >= 2:
                d[p[0].strip()] = p[1].strip()
    return d


def build_contig_coverage(faa_files, bam_map, depth_map, out_dir):
    cov_dir = out_dir / "_coverage"
    cov_dir.mkdir(exist_ok=True)
    gc = {}
    for faa in faa_files:
        label = faa.name
        stem  = faa.stem
        dp = depth_map.get(label) or depth_map.get(stem)
        if dp:
            print(f"  [INFO] Loading depth for {stem}…")
            gc[label] = load_depth_file(dp)
            continue
        bp = bam_map.get(label) or bam_map.get(stem)
        if bp:
            print(f"  [INFO] Computing coverage for {stem}…")
            gc[label] = compute_coverage_from_bam(bp, cov_dir, stem)
    return gc


def compute_tpm(genome_cov, contig_lens):
    rates = {}
    for c, info in genome_cov.items():
        length  = contig_lens.get(c, info.get("length", 1)) or 1
        rates[c] = info["depth"] / length
    total = sum(rates.values()) or 1.0
    return {c: (r / total) * 1e6 for c, r in rates.items()}
