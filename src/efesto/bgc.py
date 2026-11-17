# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Efesto contributors.

"""antiSMASH BGC region parsing for siderophore/metallophore overlap evidence.

Parses antiSMASH GFF3 output (produced with --output-format gff3) to identify
siderophore-type BGC regions. Overlapping retained ORFs are reported as observations, without numerical weights.

Only siderophore/metallophore-type regions are retained. Other antiSMASH region types
(NRPS, PKS, terpene, ...) are ignored.

Usage:
    bgc_map = parse_antismash_gff(bgc_dir, faa_files)
"""

from pathlib import Path


_SIDEROPHORE_KEYWORDS = ("siderophore", "metallophore")


def _is_siderophore_type(bgc_type):
    t = bgc_type.lower()
    return any(k in t for k in _SIDEROPHORE_KEYWORDS)


def parse_antismash_gff(bgc_dir, faa_files, return_status=False):
    """
    Find siderophore BGC regions from antiSMASH GFF3 files.

    Searches for <stem>.gff3 or <stem>.gff in bgc_dir, then in bgc_dir/<stem>/.
    Parses only features of type 'region' or 'cluster'. Keeps only those whose
    'product' or 'type' attribute matches siderophore keywords.

    Args:
        bgc_dir:   path to directory containing antiSMASH output GFF3 files
        faa_files: list of Path objects for genome FAA files (used for stem matching)

    Returns:
        {genome_faa_name: [(contig, start, end, bgc_type)]}
        Empty dict if bgc_dir is None or no matching files found.
    """
    if bgc_dir is None:
        return ({}, {}) if return_status else {}
    bgc_dir = Path(bgc_dir)
    result  = {}
    statuses = {}

    for faa in faa_files:
        stem     = faa.stem
        gff_path = None
        for ext in (".gff3", ".gff"):
            for candidate in (bgc_dir / (stem + ext),
                              bgc_dir / stem / (stem + ext)):
                if candidate.exists():
                    gff_path = candidate
                    break
            if gff_path:
                break

        if gff_path is None:
            statuses[faa.name] = "input_unavailable"
            continue

        statuses[faa.name] = "parsed"
        regions = []
        try:
            with open(gff_path) as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    parts = line.rstrip().split("\t")
                    if len(parts) < 9:
                        continue
                    if parts[2] not in ("region", "cluster"):
                        continue
                    try:
                        start = int(parts[3])
                        end   = int(parts[4])
                    except ValueError:
                        continue
                    contig   = parts[0]
                    bgc_type = _extract_bgc_type(parts[8])
                    if _is_siderophore_type(bgc_type):
                        regions.append((contig, start, end, bgc_type))
        except OSError as e:
            statuses[faa.name] = "failed"
            import sys
            print(f"  [WARN] BGC: could not read {gff_path}: {e}", file=sys.stderr)
            continue

        if regions:
            result[faa.name] = regions

    n_genomes = len(result)
    n_regions = sum(len(v) for v in result.values())
    if n_genomes:
        import sys
        print(f"[INFO] BGC: {n_regions} siderophore region(s) in {n_genomes} genome(s)")
    return (result, statuses) if return_status else result


def _extract_bgc_type(attr_str):
    """Extract BGC type from GFF3 attribute string. Checks 'product=' then 'type='."""
    for attr in attr_str.split(";"):
        attr = attr.strip()
        if attr.startswith("product="):
            return attr[8:].strip()
    for attr in attr_str.split(";"):
        attr = attr.strip()
        if attr.startswith("type="):
            return attr[5:].strip()
    return ""
