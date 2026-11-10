# Confidence scoring

```text
cluster_confidence = min(1.2,
  hmm_weight × co_occ_score × uniop_weight × bgc_boost
  × annot_weight × struct_weight × topo_weight)
```

This is a **heuristic evidence score**, not a calibrated probability. The
components can be correlated; they should not be described as seven statistically
independent tests. The normal pipeline does not apply a global confidence cutoff.

## Components

| Component | Implemented rule |
|---|---|
| HMM | Mean hit weight: 1.0 for `calibrated`, 0.5 for `low_confidence` |
| Co-occurrence | Completeness × mean distance decay × contig-edge factor |
| UniOP | Minimum available pairwise probability; neutral 1.0 if unavailable |
| BGC | 1.2 for overlap with a qualifying antiSMASH siderophore/metallophore region; otherwise 1.0 |
| Annotation | Applicable eggNOG confirmation: 1.2 supporting, 0.5 contradicted, otherwise 1.0 |
| Structure | Applicable Foldseek confirmation: 1.2 supporting, 0.5 contradicted, otherwise 1.0 |
| Topology | Applicable topology/motif verdict: 1.2 supporting, 0.5 contradicted, otherwise 1.0 |

`calibrated` is a runtime label for a positive model-specific cutoff. It does not
certify that the threshold was independently validated by Efesto. A zero-cutoff
model uses the generic search gates and is labelled `low_confidence`.

## Co-occurrence

For multi-gene clusters, completeness is `min(n_observed / canonical_size, 1)`
when an expected size is available; otherwise it is neutral. Distance decay is
`exp(-ln(2) × max(gap, 0) / 500)` and is averaged over consecutive available
gene pairs. The edge factor is 0.7 when an applicable gene lies within 3,000 bp
of a known contig boundary; otherwise it is 1.0.

Singletons return a co-occurrence score of 1.0 because genomic co-occurrence is
not defined for one gene. Missing coordinates or contig lengths make some
penalties unavailable; neutrality is not supporting evidence.

## Scope of optional checks

Annotation and structural confirmation use the `needs_confirmation` subset
that failed first-pass context filtering, rather than all matches. Informative
eggNOG support can rescue these hits. Foldseek corroborates retained/rescued
results through a weight; it does not independently rescue a dropped hit.
Topology checking uses models with explicit
`expected_topology` / `expected_heme_c` metadata and can apply outside that subset.

For each applicable evidence type, a contradiction can reduce the cluster weight.
No result, an uninformative result or a tool that did not run is neutral. A value
above 1.0 can arise from more than one supporting component; antiSMASH is not
the sole possible boost.

## Reporting

The statistics output groups scores into high (≥0.8), medium (0.5–<0.8) and low
(<0.5) bins. These are descriptive tiers, not experimentally validated probability
bands. Report optional tools, unavailable evidence and any downstream threshold
with your results.
