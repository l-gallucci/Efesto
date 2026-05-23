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


def cluster_by_coordinates(orf_set, orf_coords, max_bp_gap=5000, strand_aware=False):
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
                if sg[i][0] - sg[i - 1][1] <= max_bp_gap:
                    group.append(sg[i][3])
                else:
                    clusters.append(group)
                    group = [sg[i][3]]
            clusters.append(group)
    return clusters


def build_clusters(genome, orf_hits, orf_coords, max_gap, max_bp_gap, strand_aware):
    if orf_coords:
        return cluster_by_coordinates(
            orf_hits.keys(), orf_coords,
            max_bp_gap=max_bp_gap, strand_aware=strand_aware)
    return cluster_by_index(orf_hits.keys(), max_gap=max_gap)
