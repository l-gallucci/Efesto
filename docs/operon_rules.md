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

Use `Efesto-filter-decisions.tsv` to inspect actual rule actions on selected HMM
candidates, including exclusions and later rescue. Counts and final annotations
still use retained hits only. See [output files](output_formats.md#filter-decision-table).

Without `operon_rules.json` in `--hmm_dir`, the FeGenie-derived engine is used,
including Efesto's additional rules. The default engine evaluates every applicable
system independently and filters only that system's members when its requirement
is not met. A complete FLEET neighborhood therefore does not rescue an isolated
SufA or an unsupported siderophore-synthesis hit. A JSON file in that directory selects the
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

In the configurable JSON engine, `on_fail=keep_all` retains the current rows;
`drop` drops the current cluster when the applicable rule fails. For named-gene
count/anchor rules, `passthrough_non_members` removes the rule's members and lets
the remaining rows continue through subsequent rules. `max_bp_gap` provides a
neighborhood threshold. The engine ignores `canonical_size` metadata and does
not calculate a completeness score. Inspect
the supplied JSON and [evidence reference](evidence.md) when changing rules.

CntL+CntM and the MbnA/B/C set have explicit distinct-family requirements.
`--all_results` bypasses context filtering but still reports only the best accepted
HMM per ORF. `--relaxed_operons` relaxes supported minimum-count rules on short
contigs in the JSON engine; anchor requirements are not automatically relaxed.

## Read rule outcomes

The [context-evidence table](evidence.md#neighborhood-observations) reports rule
requirements against the retained neighborhood. It distinguishes satisfaction,
bypass and missing information, and counts model families separately from copies.
It does not replay every branch taken by the filtering engine.

For example, observing SufB/SufC/SufS can satisfy the default three-model minimum.
This is not a six-component completeness claim. The Mtr/Mto rule assigns a
context-dependent category; it is not a general completeness definition.

## UniOP

```bash
efesto --fna_dir genomes/ --meta --out efesto_out \
  --operon_prediction --uniop_path /path/to/UniOP/src/UniOP
```

UniOP adds optional operon assignments and available pairwise probabilities.
Efesto passes Prodigal-compatible FAA headers with `--faa_source prodigal_faa`;
its own caller remains Pyrodigal/Pyrodigal-GV. The `prodigal_faa` source name is
an input-format declaration, not a request to switch Efesto's gene caller.

`uniop.operon` provides membership; `uniop.pred` provides pairwise probabilities.
Efesto reads both when available, retaining membership from the operon file and
reporting the prediction file as pairwise evidence. Zero probabilities remain zero. The
prediction file is also used to assemble fallback membership when needed.

When usable pairwise probabilities are unavailable, no probabilities are invented.
UniOP caches bind results to input contents (including FAA order and Bakta GFF),
input mode and the resolved script/interpreter. Results retain the exact FAA index
used for prediction, together with content hashes of the index and output files.
Legacy or mismatched caches are recomputed; failed runs cannot publish a completed
cache. FAA reconstruction requires matching coordinates for every protein.
Inspect the log and retained UniOP working directories. See [output files](output_formats.md).
