"""Cluster confidence scoring: co-occurrence, UniOP pairwise probability, BGC boost.

Independent evidence sources combined multiplicatively:

    cluster_confidence = hmm_weight × co_occ_weight × uniop_weight
                          × bgc_boost × annot_weight × struct_weight × topo_weight

- hmm_weight:   fraction of hits with calibrated (equivalog) HMM cutoffs
- co_occ_weight: gene completeness × distance-decay penalty × contig-edge penalty
- uniop_weight:  minimum pairwise UniOP operon probability in the cluster (weakest link)
- bgc_boost:    optional ×1.2 if cluster overlaps an antiSMASH siderophore BGC
- annot_weight: eggNOG-mapper confirmation for needs_confirmation-flagged hits
                (tier 2 — see docs/hmm_library_curation.md, "SUF operon deployment")
- struct_weight: Baktfold structural confirmation for the same flagged hits
                (tier 3 — placeholder parameter, always 1.0 until implemented)
- topo_weight:  Phobius topology + heme-motif-count consistency, computed for
                EVERY hit (not just needs_confirmation-flagged ones) — see
                docs/cluster_confidence_scoring.md, Component 7

annot_weight, struct_weight, and topo_weight all default to 1.0 (neutral) —
same no-data-no-penalty policy as uniop_weight: they only affect genes with
an applicable expectation (needs_confirmation flag for the first two,
expected_heme_c/expected_topology registry columns for topo_weight), and
only when the corresponding tool actually ran and returned an informative
result for that specific ORF.

Default call threshold: cluster_confidence ≥ 0.5
"""

import math


# Half-life for distance decay (bp).
# 500 bp reflects the ~90th-percentile of within-operon intergenic distances
# in bacteria (Zhang et al. 2006, doi:10.1016/j.compbiolchem.2006.03.002).
# Gaps beyond this are increasingly penalized; fragments in MAGs are
# penalized but not eliminated.
_HALF_LIFE_BP = 500


def distance_decay(gap_bp, half_life=_HALF_LIFE_BP):
    """Exponential decay: 1.0 at gap=0, 0.5 at gap=half_life bp."""
    return math.exp(-math.log(2) * max(gap_bp, 0) / half_life)


def co_occurrence_score(cluster_orfs, orf_coords, contig_lens,
                        canonical_size=None, edge_margin=3000):
    """
    Score co-occurrence quality of a cluster:

        score = completeness × distance_factor × edge_factor

    Args:
        cluster_orfs:   iterable of ORF names in cluster
        orf_coords:     {orf: {"start": int, "end": int, "contig": str}} or {}
        contig_lens:    {contig: int} or {}
        canonical_size: expected number of genes in the complete system;
                        None → completeness = 1.0 (not penalized)
        edge_margin:    bp from contig edge to flag possible truncation

    Returns float in [0, 1].
    Single-gene clusters return 1.0 (co-occurrence not applicable).
    """
    orfs = list(cluster_orfs)
    if len(orfs) <= 1:
        return 1.0

    completeness = (min(len(orfs) / canonical_size, 1.0)
                    if canonical_size and canonical_size > 0 else 1.0)

    if orf_coords:
        def _start(o):
            c = orf_coords.get(o)
            return c["start"] if c else 0
        sorted_orfs = sorted(orfs, key=_start)
    else:
        sorted_orfs = orfs

    gaps = []
    edge_hit = False
    for i, orf in enumerate(sorted_orfs):
        c = orf_coords.get(orf) if orf_coords else None
        if c is None:
            continue
        if i > 0:
            prev = orf_coords.get(sorted_orfs[i - 1])
            if prev:
                gap = max(0, c["start"] - prev["end"])
                gaps.append(gap)
        length = contig_lens.get(c["contig"], 0) if contig_lens else 0
        if length and (c["start"] < edge_margin or (length - c["end"]) < edge_margin):
            edge_hit = True

    distance_factor = (
        sum(distance_decay(g) for g in gaps) / len(gaps) if gaps else 1.0)
    edge_factor = 0.7 if edge_hit else 1.0

    return completeness * distance_factor * edge_factor


def uniop_pair_score(cluster_orfs, pair_probs):
    """
    Minimum pairwise UniOP operon probability among all pairs in the cluster.

    Uses the minimum (weakest-link principle): a single low-probability pair
    indicates the cluster spans an operon boundary. Returns 1.0 (neutral) when
    UniOP was not run or no pair scores are available for this cluster.

    Args:
        cluster_orfs: iterable of ORF names
        pair_probs:   {(orf_a, orf_b): float} from _parse_uniop_pred;
                      empty when uniop.operon (not uniop.pred) was parsed
    """
    if not pair_probs or len(list(cluster_orfs)) <= 1:
        return 1.0
    orfs = list(cluster_orfs)
    scores = []
    for i, a in enumerate(orfs):
        for b in orfs[i + 1:]:
            p = pair_probs.get((a, b)) or pair_probs.get((b, a))
            if p is not None:
                scores.append(p)
    return min(scores) if scores else 1.0


def hmm_weight(cluster_rows):
    """
    Mean confidence weight across all rows in a cluster.

    calibrated hit → 1.0,  low_confidence hit → 0.5.
    Returns 1.0 for empty clusters.
    """
    if not cluster_rows:
        return 1.0
    weights = [1.0 if r.get("confidence") == "calibrated" else 0.5
               for r in cluster_rows]
    return sum(weights) / len(weights)


def annot_weight(cluster_rows, eggnog_hits, confirmation_map):
    """
    eggNOG-mapper tier-2 confirmation weight for needs_confirmation-flagged hits.

    Weakest-link, same principle as uniop_weight: a single contradicted hit in
    the cluster is enough to lower confidence, since the whole point of this
    check is catching a specific false-positive mode (e.g. sufE HMM firing on
    csdE, a related but non-operon paralog — see hmm_library_curation.md).

    Args:
        cluster_rows: rows with "hmm_stem" and "orf" keys
        eggnog_hits:  {orf_id: eggnog_row_dict} from _parse_eggnog_annotations,
                      or {} if eggNOG wasn't run / no file was provided
        confirmation_map: {stem: {"gene_name": str, "kegg_ko": [str, ...]}}
                      from io.build_confirmation_map — only stems in this map
                      are ever checked; everything else is untouched by design

    Returns 1.0 (neutral) when:
        - no row in the cluster is on a needs_confirmation-flagged stem
        - eggNOG wasn't run at all (eggnog_hits is empty)
        - the flagged ORF has no eggNOG hit, or an uninformative one

    Otherwise: 1.2 if every flagged hit in the cluster is eggNOG-confirmed,
    0.5 if any flagged hit is eggNOG-contradicted (weakest link).
    """
    from efesto.eggnog import confirm_hit

    verdicts = []
    for r in cluster_rows:
        stem = r.get("hmm_stem", "")
        info = confirmation_map.get(stem)
        if info is None:
            continue
        eggnog_row = eggnog_hits.get(r.get("orf", ""))
        verdicts.append(confirm_hit(info["gene_name"], info["kegg_ko"], eggnog_row))

    if not verdicts or all(v == "neutral" for v in verdicts):
        return 1.0
    if any(v == "contradicted" for v in verdicts):
        return 0.5
    return 1.2


def struct_weight(cluster_rows, foldseek_hits, confirmation_map):
    """
    Baktfold (Foldseek+ProstT5) tier-3 confirmation weight for
    needs_confirmation-flagged hits — same scope and weakest-link principle
    as annot_weight, just a structural signal instead of an orthology one.

    Args:
        cluster_rows: rows with "hmm_stem" and "orf" keys
        foldseek_hits: {orf_id: foldseek_row} from
                       baktfold._parse_foldseek_results, or {} if Foldseek
                       wasn't run / no file was provided
        confirmation_map: {stem: {"gene_name": str, "kegg_ko": [str, ...]}}
                      from io.build_confirmation_map — only stems in this map
                      are ever checked; everything else is untouched by design
                      (same map tier-2/annot_weight already uses)

    Returns 1.0 (neutral) when:
        - no row in the cluster is on a needs_confirmation-flagged stem
        - Foldseek wasn't run at all (foldseek_hits is empty)
        - the flagged ORF has no structural hit

    Otherwise: 1.2 if every flagged hit in the cluster is Foldseek-confirmed,
    0.5 if any flagged hit is Foldseek-contradicted (weakest link).
    """
    from efesto.baktfold import confirm_struct_hit

    verdicts = []
    for r in cluster_rows:
        stem = r.get("hmm_stem", "")
        info = confirmation_map.get(stem)
        if info is None:
            continue
        foldseek_row = foldseek_hits.get(r.get("orf", ""))
        verdicts.append(confirm_struct_hit(info["gene_name"], foldseek_row))

    if not verdicts or all(v == "neutral" for v in verdicts):
        return 1.0
    if any(v == "contradicted" for v in verdicts):
        return 0.5
    return 1.2


def topo_weight(cluster_rows, phobius_hits, topology_map):
    """
    Topology/motif consistency weight — weakest-link across the cluster, same
    principle as annot_weight. Unlike annot_weight/struct_weight, this checks
    EVERY hit whose stem has an expectation set in the registry, not just
    needs_confirmation-flagged ones (see io.build_topology_expectation_map).

    Two independent checks per applicable row: Phobius topology
    (topology.topology_verdict) and heme-motif count
    (topology.heme_verdict, against the heme_motifs value already computed
    by operon.count_heme for every hit).

    Args:
        cluster_rows: rows with "hmm_stem", "orf", "heme_motifs" keys
        phobius_hits: {orf_id: phobius_row} from
                      topology._parse_phobius_output, or {} if Phobius
                      wasn't run (heme_verdict still applies regardless,
                      since heme_motifs is already computed for every hit)
        topology_map: {stem: {"expected_heme_c": int|None,
                       "expected_topology": str|None}} from
                       io.build_topology_expectation_map — only stems in
                       this map are ever checked; everything else is
                       untouched by design

    Returns 1.0 (neutral) when:
        - no row in the cluster is on a stem with any expectation set
        - all applicable checks come back neutral (e.g. Phobius wasn't run
          and heme count is inconclusive)

    Otherwise: 1.2 if every applicable verdict confirms, 0.5 if any
    contradicts (weakest link).
    """
    from efesto.topology import heme_verdict, topology_verdict

    verdicts = []
    for r in cluster_rows:
        stem = r.get("hmm_stem", "")
        info = topology_map.get(stem)
        if info is None:
            continue
        phobius_row = phobius_hits.get(r.get("orf", ""))
        verdicts.append(topology_verdict(info["expected_topology"], phobius_row))
        verdicts.append(heme_verdict(info["expected_heme_c"], r.get("heme_motifs")))

    if not verdicts or all(v == "neutral" for v in verdicts):
        return 1.0
    if any(v == "contradicted" for v in verdicts):
        return 0.5
    return 1.2


def cluster_confidence(hmm_w, co_occ_w, uniop_w, bgc_boost=1.0,
                        annot_w=1.0, struct_w=1.0, topo_w=1.0):
    """
    Combined cluster confidence score.

    Returns float, capped at 1.2 (bgc_boost/annot_w/topo_w can push above
    1.0). Recommended call threshold: ≥ 0.5.

    Args:
        hmm_w:      hmm_weight(cluster_rows)
        co_occ_w:   co_occurrence_score(...)
        uniop_w:    uniop_pair_score(...), 1.0 if UniOP not run
        bgc_boost:  1.2 if in siderophore BGC, 1.0 otherwise
        annot_w:    annot_weight(...), 1.0 if eggNOG not run or not applicable
        struct_w:   Baktfold tier-3 confirmation, 1.0 until implemented
        topo_w:     topo_weight(...), 1.0 if not applicable — computed for
                    every hit, not just needs_confirmation-flagged ones
    """
    return min(1.2, hmm_w * co_occ_w * uniop_w * bgc_boost
               * annot_w * struct_w * topo_w)
