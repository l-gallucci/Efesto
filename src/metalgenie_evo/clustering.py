"""ORF clustering by genomic index or base-pair coordinates."""

from collections import defaultdict


def _index_from_name(orf):
    parts = orf.rsplit("_", 1)
    if len(parts) == 2:
        try:
            return parts[0], int(parts[1])
        except ValueError:
            pass
    return orf, 0


def _orf_to_contig(orf):
    c, _ = _index_from_name(orf)
    return c


def cluster_by_index(orf_set, max_gap=5):
    by_c = defaultdict(list)
    for orf in orf_set:
        c, i = _index_from_name(orf)
        by_c[c].append((i, orf))
    clusters = []
    for c, entries in by_c.items():
        entries.sort(key=lambda x: x[0])
        group = [entries[0][1]]
        for i in range(1, len(entries)):
            if entries[i][0] - entries[i - 1][0] <= max_gap:
                group.append(entries[i][1])
            else:
                clusters.append(group)
                group = [entries[i][1]]
        clusters.append(group)
    return clusters


def cluster_by_coordinates(orf_set, orf_coords, max_bp_gap=5000, strand_aware=False,
                           orf_stem_map=None, stem_gap_map=None):
    """
    Cluster ORFs by base-pair proximity.

    When stem_gap_map is provided, the effective gap threshold for joining two
    consecutive ORFs is min(gap_a, gap_b) where gap_x is stem_gap_map[stem_x]
    (falling back to max_bp_gap when the stem has no per-rule value). This makes
    tight gene systems (MtrMto: 2000 bp, FoxABC: 1000 bp) stricter than the
    global window without breaking detection of loosely encoded systems
    (SIDERO_SYNTH: 5000 bp).
    """
    by_c = defaultdict(list)
    for orf in orf_set:
        c = orf_coords.get(orf)
        if c is None:
            contig, idx = _index_from_name(orf)
            by_c[contig].append((idx * 300, idx * 300, "+", orf))
        else:
            by_c[c["contig"]].append((c["start"], c["end"], c["strand"], orf))
    clusters = []
    for contig, entries in by_c.items():
        entries.sort(key=lambda x: x[0])
        if strand_aware:
            sg_dict = defaultdict(list)
            for e in entries:
                sg_dict[e[2]].append(e)
            strand_groups = list(sg_dict.values())
        else:
            strand_groups = [entries]
        for sg in strand_groups:
            if not sg:
                continue
            group = [sg[0][3]]
            for i in range(1, len(sg)):
                if stem_gap_map is not None and orf_stem_map is not None:
                    gap_a = stem_gap_map.get(
                        orf_stem_map.get(sg[i - 1][3], ""), max_bp_gap)
                    gap_b = stem_gap_map.get(
                        orf_stem_map.get(sg[i][3], ""), max_bp_gap)
                    effective_gap = min(gap_a, gap_b)
                else:
                    effective_gap = max_bp_gap
                if sg[i][0] - sg[i - 1][1] <= effective_gap:
                    group.append(sg[i][3])
                else:
                    clusters.append(group)
                    group = [sg[i][3]]
            clusters.append(group)
    return clusters


def build_clusters(genome, orf_hits, orf_coords, max_gap, max_bp_gap, strand_aware,
                   stem_gap_map=None):
    """
    Entry point for clustering ORF hits within a genome.

    stem_gap_map: {hmm_stem → max_bp_gap} from operon rules. When provided,
    overrides the global max_bp_gap per ORF pair using the more restrictive of
    the two ORFs' per-rule gaps. ORFs whose stem has no rule entry keep the
    global max_bp_gap.
    """
    if orf_coords:
        orf_stem_map = (
            {orf: hit["hmm_stem"] for orf, hit in orf_hits.items()}
            if stem_gap_map else None)
        return cluster_by_coordinates(
            orf_hits.keys(), orf_coords,
            max_bp_gap=max_bp_gap, strand_aware=strand_aware,
            orf_stem_map=orf_stem_map, stem_gap_map=stem_gap_map)
    return cluster_by_index(orf_hits.keys(), max_gap=max_gap)
