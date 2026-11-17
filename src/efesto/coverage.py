# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Efesto contributors.

"""Streaming BAM/depth-file coverage loading and TPM normalisation."""

import math
import subprocess
import sys
import tempfile
from collections import defaultdict


def _number(value):
    value = float(value)
    if not math.isfinite(value) or value < 0:
        raise ValueError("coverage values must be finite and non-negative")
    return value


def _parse_depth_lines(lines, source):
    """Read one sample: native depth positions, mean-depth tables, or coverage."""
    result, sums, counts, last_pos = {}, defaultdict(float), defaultdict(int), {}
    mode = None
    columns = None
    for number, line in enumerate(lines, 1):
        if not line.strip():
            continue
        parts = line.strip().split()
        header = [p.lstrip("#").lower() for p in parts]
        first = header[0]
        if first == "contigname" and "totalavgdepth" in header:
            mode, columns = "mean", (header.index("totalavgdepth"), header.index("contiglen"))
            continue
        if "meandepth" in header and first == "rname":
            mode = "coverage"
            continue
        if first in ("chrom", "chr", "contig", "#chrom") and len(header) >= 3 and header[1] in ("pos", "position"):
            if len(header) != 3:
                raise ValueError(f"{source}:{number}: expected depth from one BAM (three columns)")
            mode = "positions"
            continue
        if first in ("id", "contig", "contigname", "rname") and len(header) >= 2 and header[1] in ("avg_fold", "depth", "mean_depth", "meandepth"):
            mode = "mean"
            columns = (1, header.index("length") if "length" in header else None)
            continue
        if line.lstrip().startswith("#"):
            continue
        try:
            if mode is None:
                if len(parts) == 2:
                    mode, columns = "mean", (1, None)
                elif len(parts) == 3:
                    mode = "positions"
                elif len(parts) >= 7:
                    mode = "coverage"
                else:
                    raise ValueError("unsupported depth format; expected one sample")
            contig = parts[0]
            if mode == "positions":
                if len(parts) != 3:
                    raise ValueError("expected contig, position, depth for one BAM")
                pos = int(parts[1])
                if pos < 1 or pos <= last_pos.get(contig, 0):
                    raise ValueError("positions must be positive and increasing within each contig")
                depth = _number(parts[2])
                last_pos[contig] = pos
                sums[contig] += depth
                counts[contig] += 1
            else:
                if contig in result:
                    raise ValueError("duplicate contig in mean-depth table")
                if mode == "coverage":
                    if len(parts) < 7:
                        raise ValueError("incomplete samtools coverage row")
                    start, end = int(parts[1]), int(parts[2])
                    if start < 1 or end < start:
                        raise ValueError("invalid coverage interval")
                    depth, reads, length = _number(parts[6]), int(parts[3]), end - start + 1
                else:
                    dcol, lcol = columns
                    depth, reads = _number(parts[dcol]), 0
                    length = int(parts[lcol]) if lcol is not None else 0
                if length < 0 or reads < 0:
                    raise ValueError("length and read counts cannot be negative")
                result[contig] = {"depth": depth, "reads": reads, "length": length}
        except (ValueError, IndexError) as error:
            raise ValueError(f"{source}:{number}: {error}") from error
    for contig, total in sums.items():
        result[contig] = {"depth": total / counts[contig], "reads": 0,
                          "length": last_pos[contig]}
    return result


def compute_coverage_from_bam(bam_path, out_dir, label):
    table = out_dir / f"{label}.coverage"
    result = subprocess.run(
        ["samtools", "coverage", "-H", "-o", str(table), str(bam_path)],
        capture_output=True, text=True)
    if result.returncode == 0 and table.is_file():
        return load_depth_file(table)
    print(f"  [WARN] samtools coverage failed for {label}, trying depth", file=sys.stderr)
    # Stream stdout and spool stderr to disk: neither the full per-base table
    # nor a full pipe of diagnostic output can accumulate in Python memory.
    with tempfile.TemporaryFile(mode="w+") as errors:
        process = subprocess.Popen(["samtools", "depth", "-aa", str(bam_path)],
                                   stdout=subprocess.PIPE, stderr=errors, text=True)
        try:
            with process.stdout:
                coverage = _parse_depth_lines(process.stdout, str(bam_path))
        except BaseException:
            process.kill()
            process.wait()
            raise
        code = process.wait()
        if code:
            errors.seek(0)
            raise RuntimeError(f"samtools depth failed for {label} (exit {code}): {errors.read(2000).strip()}")
    return coverage


def load_depth_file(path):
    with open(path) as handle:
        return _parse_depth_lines(handle, str(path))


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


def compute_coverage_per_million(genome_cov):
    """Scale contig mean depths to 1e6 within one genome/sample input.

    Mean depth is already length-normalised. This is a DNA coverage metric,
    not transcript TPM, and must not divide by contig length a second time.
    Category outputs sum the assigned value once per retained gene hit.
    """
    depths = {contig: _number(info["depth"]) for contig, info in genome_cov.items()}
    total = math.fsum(depths.values())
    if total == 0:
        return {contig: 0.0 for contig in depths}
    return {contig: depth / total * 1e6 for contig, depth in depths.items()}


def compute_tpm(genome_cov, contig_lens=None):
    """Compatibility alias; returns coverage per million, independent of length."""
    return compute_coverage_per_million(genome_cov)
