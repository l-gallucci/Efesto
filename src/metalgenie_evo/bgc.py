"""antiSMASH BGC region parsing for siderophore cluster boost.

Parses antiSMASH GFF3 output (produced with --output-format gff3) to identify
siderophore-type BGC regions. Clusters whose ORFs overlap a siderophore region
receive a bgc_boost=1.2 multiplier in cluster_confidence.

Only siderophore-type regions trigger the boost. Other antiSMASH region types
(NRPS, PKS, terpene, ...) are ignored to avoid inflating non-siderophore clusters.

Usage:
    bgc_map = parse_antismash_gff(bgc_dir, faa_files)
    boost   = bgc_boost_for_cluster(cluster_orfs, orf_coords, bgc_map.get(genome))
"""

from pathlib import Path


_SIDEROPHORE_KEYWORDS = ("siderophore", "metallophore")


def _is_siderophore_type(bgc_type):
    t = bgc_type.lower()
    return any(k in t for k in _SIDEROPHORE_KEYWORDS)


def parse_antismash_gff(bgc_dir, faa_files):
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
        return {}
    bgc_dir = Path(bgc_dir)
    result  = {}

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
            continue

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
    return result


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


def bgc_boost_for_cluster(cluster_orfs, orf_coords, bgc_regions, boost=1.2):
    """
    Return bgc_boost value for a cluster.

    Returns boost (default 1.2) if any ORF in the cluster overlaps a siderophore
    BGC region. Returns 1.0 otherwise.

    Args:
        cluster_orfs: iterable of ORF names
        orf_coords:   {orf: {"start": int, "end": int, "contig": str}}
        bgc_regions:  [(contig, start, end, bgc_type)] for this genome, or None/[]
        boost:        multiplier applied when overlap found
    """
    if not bgc_regions or not orf_coords:
        return 1.0

    by_contig = {}
    for contig, rs, re, _ in bgc_regions:
        by_contig.setdefault(contig, []).append((rs, re))

    for orf in cluster_orfs:
        c = orf_coords.get(orf)
        if c is None:
            continue
        for rs, re in by_contig.get(c["contig"], []):
            if c["start"] <= re and c["end"] >= rs:
                return boost
    return 1.0
