# Operon rules and UniOP

## Genomic clustering

`--max_bp_gap` controls the coordinate-based default window (5,000 bp).
`--strand_aware` separates strands. Named systems can use tighter rule-specific
gaps; the effective threshold between two hits is the more restrictive of their
applicable thresholds. Without coordinates, `--max_gap` controls ordinal clustering.
Use [catalog mode](workflows.md#gene-catalogs) when header adjacency has no genomic meaning.

A detected cluster is a neighborhood of retained HMM-positive genes. It is not
a demonstration that all intervening genes were considered or that the genes
are co-transcribed.

## Context engines

Without `operon_rules.json` in `--hmm_dir`, the FeGenie-derived engine is used,
including Efesto's additional rules. A JSON file in that directory selects the
configurable engine. The repository-root example must be copied into the chosen
library directory to use it; its mere presence at the repository root is not enough.

The JSON schema has `report_all_categories` and `rules`. Category bypasses preserve
specific reporting behavior and do not certify that individual hits are validated.
Rule `genes` entries identify HMM stems.

| Rule | Condition |
|---|---|
| `require_n_of` | Minimum number of distinct named models |
| `require_anchor` | Named anchor must be present |
| `require_n_cat` | Minimum number of distinct models in selected categories |
| `require_n_cat_or_lone_trusted` | Category count or an explicitly trusted lone model |
| `mtr_disambiguation` | Contextual Mtr/Mto category assignment |

`on_fail` can keep all rows, drop the affected set or pass nonmembers, depending
on the rule. `canonical_size` supplies a completeness reference for scoring;
`max_bp_gap` provides a neighborhood threshold. Inspect the supplied JSON and
[scoring reference](cluster_confidence_scoring.md) when changing them.

CntL+CntM and the MbnA/B/C set have explicit distinct-family requirements.
`--all_results` bypasses context filtering but still reports only the best accepted
HMM per ORF. `--relaxed_operons` relaxes supported minimum-count rules on short
contigs in the JSON engine; anchor requirements are not automatically relaxed.

## UniOP

```bash
efesto --fna_dir genomes/ --meta --out efesto_out \
  --operon_prediction --uniop_path /path/to/UniOP/src/UniOP
```

UniOP adds optional operon assignments and available pairwise probabilities.
Efesto passes Prodigal-compatible FAA headers with `--faa_source prodigal_faa`;
its own caller remains Pyrodigal/Pyrodigal-GV. The `prodigal_faa` source name is
an input-format declaration, not a request to switch Efesto's gene caller.

When usable pairwise scores are unavailable, the UniOP confidence contribution
is neutral. Cached or failed runs can affect evidence availability; inspect the
log and retained UniOP working directories. See [output files](output_formats.md).
