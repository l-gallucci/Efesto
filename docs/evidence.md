# Evidence and decisions

Efesto reports HMM matches, neighborhood observations and optional checks
separately. Read each result as a family assignment with its supporting or
unavailable observations. Evidence is not collapsed into a confidence number.

## What determines whether a hit is retained?

| Stage or check | Can change retained hits? | Current behavior |
|---|---|---|
| HMMER search gates | Yes | Require the applicable full-sequence bitscore threshold and E-value limit; keep the highest-scoring accepted model per ORF |
| Context and per-protein rules | Yes | Apply the selected neighborhood rules and applicable sequence checks |
| eggNOG candidate rescue | Yes, for a restricted candidate set | An informative matching assignment can recover a flagged candidate that failed the first-pass filters |
| UniOP | No | Report operon membership and available pair probabilities |
| antiSMASH | No | Report qualifying region overlaps for each ORF |
| Foldseek | No | Report reference matches for the confirmation-candidate set |
| Optional Phobius/motif comparisons | No | Report predictions and consistency with declared expectations |

The optional motif comparison is distinct from existing sequence filters. For
example, the Cyc2 second-pass rule still uses sequence length and heme motifs.
Those filters remain in place.

## Why was a candidate excluded?

`Efesto-filter-decisions.tsv` traces the actual filtering decisions for every
selected HMM candidate, including unflagged candidates that never enter optional
confirmation. It retains the original exclusion and the final outcome after any
eggNOG rescue. For example, a DFE_0461-family candidate can exceed its HMM cutoff
but fail a DFE2 minimum of three distinct models. Nearby Mtr-family models do not
fulfil that requirement. This reports insufficient support for that source-style
system; it does not prove that the protein has no metal-related function.

Excluded candidates appear only in the decision table, not in retained-hit
counts or sequence exports. This table is distinct from
`Efesto-context-evidence.tsv`, which assesses retained neighborhoods rather than
replaying exclusions. See the [decision-table fields](output_formats.md#filter-decision-table).

## Match statistics and search thresholds

`Efesto-results-long.tsv` has one row per retained ORF. The HMM identifier,
bitscore, E-value and model cutoff describe the accepted match.

| Field | Meaning |
|---|---|
| `threshold_type=model_specific` | A positive cutoff exists in the model cutoff table |
| `threshold_type=generic` | The model has no positive cutoff; the generic search gate is used |
| `bitscore_cutoff` | The configured model cutoff, zero when absent |
| `search_bitscore_cutoff` | The effective search threshold, including the generic fallback |
| `threshold_source` | Cutoff table or generic search option used by this run |
| `model_source` | Registry source for the model, when available |
| `retention_mode` | Ordinary context filtering, all-results/catalog mode, or eggNOG rescue |

By default, the generic bitscore floor is **30** and the E-value limit is
**1e-5**. Model-specific thresholds control search acceptance. Interpret
substrate specificity using the profile's functional scope and references.
Bitscores from different models are not a common probability scale.

The [output reference](output_formats.md) describes identifiers, tables and
export formats. Read columns by header name.

## Available, missing and discordant checks

Availability and result are different fields. A tool can complete without
producing an informative result for a particular ORF.

| Status or result | How to read it |
|---|---|
| `not_applicable` | This check does not apply to this hit or was outside the selected candidate scope |
| `not_requested` | The optional analysis was not requested |
| `not_queried` / `not_run` | No query was submitted or no analysis ran for this item |
| `input_unavailable` / `unavailable` / `index_unavailable` | Required input, usable result or index mapping was unavailable |
| `failed` | The recorded tool/input operation failed |
| `provided` / `parsed` / `completed` / `cached` | Input was supplied/read, a process completed, or a result was reused; inspect its observations too |
| `no_informative_result` | No usable comparison was available for this item |
| `confirmed` / `contradicted` / `neutral` | The applicable annotation or topology/motif comparison agrees, disagrees, or is inconclusive |

A failed or unrequested check does not establish biological absence. Successful
execution does not establish function. Some parsers return no informative result
for empty or unusable output, so also inspect warnings and retained tool files.

### Scope of optional observations

| Evidence | Reported observations | Scope |
|---|---|---|
| eggNOG | `eggnog_status`, preferred name, KO assignments, `eggnog_verdict` | `needs_confirmation` candidates rejected by the first-pass filters |
| Foldseek | `structure_status`, target, bitscore, E-value, `structure_label_match` | The same candidate set, when structural search is requested |
| Phobius | `phobius_status`, helices, signal peptide, `topology_verdict` | Retained hits with declared expectations, when requested |
| Heme motifs | Observed/expected counts and `heme_verdict` | The applicable motif comparison, when topology checking is requested |
| antiSMASH | `bgc_status`, overlap and matching region coordinates/products | Retained ORFs with available coordinates and matching input |
| UniOP | `uniop_status`, operon membership/context, pair probabilities | Genomes requested for operon prediction |

Supplying a whole-proteome eggNOG file does not make this a universal eggNOG
check of every retained hit. The additional `Efesto-confirmation-candidates.tsv`
reports the candidate set, including candidates that were not recovered.
`retention_mode=not_retained` distinguishes these from `eggnog_rescue`.
Only retained hits contribute to ordinary annotation summaries.

Foldseek target labels are compared as `matches_label`, `different_label` or
`unresolved`. A matching label is not a functional confirmation, and a different
label is not a demonstrated functional contradiction. Inspect the raw match and
the reference identity. Foldseek, topology and BGC observations do not act as
additional retention filters.

## Neighborhood observations

`Efesto-context-evidence.tsv` has one row per applicable rule and retained cluster.
When no rule applies, a `not_applicable` record preserves the neighborhood's
observations. Its scope is `retained_hits_not_filter_trace`: this is an assessment
of the retained set, not a replay of every filtering decision.

The table records model identifiers, distinct models, separate copy counts,
inter-hit gaps, coordinate availability, rule identifier/type and a SHA256 of
the rule configuration. It also lists rule models not observed and models that
were not searched with the chosen library/category selection.

| Rule outcome | Meaning |
|---|---|
| `satisfied_minimum` | The configured minimum of distinct models is observed |
| `anchor_observed` | The configured anchor is observed |
| `trusted_lone_observed` | The configured lone-model exception is observed |
| `not_satisfied` | The assessed requirement is not met in the retained set |
| `category_disambiguation_only` | Mtr/Mto context is used for category assignment, without a completeness claim |
| `bypassed` | The applicable reporting mode/category bypasses the requirement |
| `not_evaluable_*` | Coordinates, searched models or post-rescue context do not support that assessment |

Observing enough distinct members can satisfy a minimum even if other alternatives
were not searched. Conversely, an unsearched model is not evidence of absence.
Without coordinates, ordinal grouping may still be used by the pipeline, but the
assessment does not certify a genomic neighborhood. Catalog mode bypasses count
rules in the default engine; custom JSON rules should be inspected separately.
JSON count relaxations are reflected when applicable.

A satisfied minimum is not a complete or active pathway. The table does not
invent a universal completeness fraction from a model list. A `not_satisfied`
record does not retrospectively remove a retained hit; review it against the
[configured filtering rules](operon_rules.md) and run conditions.

### Contig boundaries

`contig_length_status` distinguishes a known assembly/sequence-region length from
a lower bound inferred from CDS coordinates. `distance_to_contig_edge_bp` is
reported only for known lengths; otherwise it is blank. The end of the last CDS
is not treated as the right end of the contig. Available lengths are loaded
independently of the short-contig filter.

Distances are observations. A locus near a MAG boundary may be truncated; no
fixed edge penalty or inferred biological absence is assigned.

## Worked examples

The following schematic examples illustrate how to read the output fields.

| Observed neighborhood | What the result says | What it does not establish |
|---|---|---|
| Three MtrA hits | `model_copy_counts` records three copies; the Mtr rule observes one distinct model type | A complete MtrCAB system |
| MtrA + two LanM hits | Mtr and LanM matches remain separate evidence; LanM is not counted as an Mtr component | Replacement of MtrB or MtrC by LanM |
| SufB, SufC and SufS with usable coordinates | The default SUF minimum of three distinct models can be `satisfied_minimum` | Presence of every component of a particular SUF variant |
| One LanM hit | Family match retained; no applicable neighborhood requirement | A complete lanthanide uptake pathway |
| Flagged candidate recovered by matching eggNOG annotation | Candidate and retained tables show `eggnog_rescue` | Restoration of a complete genomic neighborhood |
| UniOP failed for one genome | `uniop_status=failed`, context `not_evaluated`, no invented pair probabilities | Proof that the genes are outside operons |

For MtrA + two LanM hits, a neighborhood record can contain:

```text
observed_models:       ["MtrA", "lanM_NF053806"]
model_copy_counts:     {"MtrA": 1, "lanM_NF053806": 2}
observed_rule_models:  ["MtrA"]
distinct_rule_models:  1
rule_status:           category_disambiguation_only
```

The rule-specific observations use stable HMM identifiers rather than translated
gene labels, and duplicate copies do not become distinct system components.

## UniOP pairs

`Efesto-UniOP-pair-evidence.tsv` reports available probabilities only for pairs
whose members are both retained hits. It includes genome, original/exported ORF
IDs and each member's cluster. Zero is retained as zero; missing pairs are not
assigned 1. These probabilities concern operon prediction, not metal function.
An empty table has a header and no inferred scores.

## Using results in analysis scripts

The output schema is `separate-evidence-v1`, recorded in `Efesto-run.log`.
`threshold_type` describes the search threshold used for each hit. The numeric
anvi'o export, `Efesto-anvio-gene-evidence.tsv`, contains individual observations.

Read columns by header name using the [file definitions](output_formats.md).
Apply model-specific interpretations rather than a universal bitscore cutoff,
and keep unavailable observations separate from negative results. Use a fresh
output directory when changing software revision or analysis settings.
