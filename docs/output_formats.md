# Output files

Files below `--out` use the `separate-evidence-v1` schema recorded in the run log.
Main annotations describe retained hits after search and context filtering.

## Choose a table

- **Family assignments:** use `Efesto-results-long.tsv` for one row per retained ORF.
- **Filtering decisions:** use `Efesto-filter-decisions.tsv` for retained, excluded and recovered selected HMM candidates.
- **Neighborhood requirements:** use `Efesto-context-evidence.tsv` for retained-cluster observations.
- **Candidate recovery:** use `Efesto-confirmation-candidates.tsv`; some rows are not retained annotations.
- **Operon probabilities:** use `Efesto-UniOP-pair-evidence.tsv`; a row describes a pair, not one gene's function.

The [evidence guide](evidence.md) explains which stages affect retention and how
to read unavailable or discordant results.

## Main files

| File | Contents |
|---|---|
| `Efesto-summary.csv` | Detailed per-ORF results with cluster separators |
| `Efesto-geneSummary-clusters.csv` | Compact per-ORF cluster summary |
| `Efesto-heatmap-data.csv` | Category × genome counts of distinct retained clusters |
| `Efesto-results-long.tsv` | Per-ORF match statistics, threshold type and separate evidence |
| `Efesto-filter-decisions.tsv` | Every selected HMM candidate, actual filter decisions and final retention/rescue status |
| `Efesto-context-evidence.tsv` | Observations and rule assessments for retained neighborhoods |
| `Efesto-confirmation-candidates.tsv` | First-pass rejected confirmation candidates, including recovered and unretained hits |
| `Efesto-UniOP-pair-evidence.tsv` | Available pair probabilities between retained hits; header-only when none |
| `Efesto-summary-stats.tsv` | Counts, threshold types, evidence availability and runtime |
| `Efesto-hits.faa` | Retained matched proteins |
| `Efesto-hits.gff3` | Retained CDS features with coordinates |
| `Efesto-hits.fna` | Nucleotide matches when source FNA and coordinates are available |
| `Efesto-run.log` | Run metadata and output schema |

## Optional files

| File | Condition |
|---|---|
| `Efesto-coverage-heatmap.csv` | Coverage inputs supplied |
| `Efesto-OperonStructure.tsv` | Operon-prediction results available |
| `Efesto-anvio-functions.tsv` | `--anvio` |
| `Efesto-anvio-gene-evidence.tsv` | `--anvio`; raw numeric observations |
| `_tblout_cache/` | Search tables retained with `--keep_tblout` |

## Long-format table

The initial columns are:

```text
category genome contig orf original_orf gene hmm_stem bitscore evalue bitscore_cutoff
threshold_type cluster_id heme_c_motifs contig_len model_nseq
```

They are followed by search-gate provenance, retention mode, coordinate/contig-edge
observations and the individual tool statuses and results. See the complete
[evidence reference](evidence.md). Always read columns by header name.

`threshold_type` is `model_specific` or `generic`; `bitscore_cutoff` is the model
cutoff (zero when absent), and `search_bitscore_cutoff` records the actual search
threshold, including the generic fallback. `model_nseq` counts training records,
including duplicate or closely related sequences. Unknown contig lengths are zero and unknown
edge distances are blank. `contig_length_status=lower_bound` marks a length
inferred from annotated CDS coordinates, which is not used as a known contig
edge. GFF3 supplies gene coordinates.

`retention_mode` is `context_filter`, `all_results`, `catalog` or `eggnog_rescue`.
In the separate candidate table, `not_retained` identifies candidates that were
not recovered. Candidate rows never contribute to normal annotation summaries.

## Filter-decision table

`Efesto-filter-decisions.tsv` is always written with its own `filter-decisions-v1`
schema. It has one row per highest-scoring accepted HMM candidate per ORF, including
candidates excluded by the contig-length, context or second-pass filters. It does
not enumerate subthreshold matches or the competing models discarded by best-hit
selection. Excluded rows do not contribute to annotations, counts, heatmaps or FASTA exports.

The table records genome/original/mapped ORF identity, model and threshold
provenance, `category_before_filter`, `category_after_context`, `final_category`,
`raw_cluster_id` and `final_cluster_id`. Raw clusters precede filtering; a rescued
candidate can acquire a different final cluster. Candidates removed before
clustering have an empty raw cluster ID. Genomic coordinates, index-based fallback
and catalog-singleton grouping are distinguished. `contig_length_status` identifies
known lengths versus CDS-derived lower bounds. Only known lengths define
contig-edge distances.

`contig_length_decision`, `context_decision` and `sequence_decision` describe the
initial filtering stages; the last field covers the existing second-pass rules,
including Cyc1 context and Cyc2 length/motif requirements. A stage can be excluded,
retained, bypassed, not applied or not evaluated, as applicable. `final_decision`
is `retained`, `excluded` or `rescued`; `retention_mode` identifies eggNOG recovery
where it occurred. The original exclusion remains visible after rescue.

`rule_decisions` contains a JSON array of decisions recorded by the executing
filter: stage, rule ID/type, requirement status, action on this candidate, reason
code, observed models and applicable minima/anchors or sequence measurements.
An unsatisfied JSON rule with `on_fail=keep_all` can have `status=failed` and
`action=retained`. A JSON cluster-wide drop also explains exclusions of nonmember
hits. `reason_codes` lists initial excluding rules and any subsequent eggNOG rescue.
Category assignments are recorded separately from exclusion decisions.
Interpret retained candidates using their model scope and available evidence.

The rule-version hash fingerprints the configuration, not the software source.
Keep the source revision and run log for reproducibility. `--all_results` records
context/second-pass bypasses; an explicitly requested short-contig filter still
applies. Decision records are spooled to disk to limit memory use.

## GFF3 and sequence headers

GFF3 uses CDS features and 1-based inclusive coordinates; score and phase are `.`.
Attributes include gene, category, model, cluster and `threshold_type`. Records
lacking coordinates are omitted from GFF3. Protein and nucleotide FASTA headers
also include `threshold_type`.

## Statistics and interpretation

Statistics use `section`, `metric`, `value`, `notes` and include RUN, THRESHOLD,
EVIDENCE, CATEGORY and GENOME sections. See
[interpretation](interpretation.md) for the relationship between retained hits,
context requirements and biological function.

## Protein identities across genomes

A protein is identified by **genome + original ORF ID**. The long-format table
keeps genome and the original/mapped ORF label in separate columns; the ID map
records the original ORF even when a Bakta label is displayed. Combined FASTA, GFF and Anvi'o exports preserve
original/Bakta labels when unique. If labels collide, every exported protein ID
in that run uses `genome::original_orf`, for example
`sampleA.faa::contig_1`. Reserved characters in each component are percent-encoded.

`Efesto-protein-id-map.tsv` links `query_id`, `export_id`, `genome`, `orf`,
`bakta_id` and `contig`. Dropped confirmation candidates have an empty `export_id`.
Queries sent to optional confirmation/topology tools always use scoped IDs.
External eggNOG annotations with original unscoped IDs remain usable only when
those IDs identify a single protein across **all** input genomes. Ambiguous IDs
stop the run with an explanation; use the scoped candidate FASTA to annotate them.

Merged GFF files also scope contig names when they occur in multiple genomes.
The `genome`, `original_contig` and `original_orf` attributes preserve the original
coordinate system. Use the source assembly contig names for extraction, or rename
assembly headers to the corresponding scoped names when importing this merged
GFF into a viewer. Anvi'o still requires mapping exported labels to that target
database's integer gene IDs.

## Caches and failures

With `--keep_tblout`, Efesto retains completed HMMER tables and adjacent
`.cache.json` files. Reuse requires matching protein/profile contents, score
threshold, HMMER executable and thread setting, plus a complete, unchanged table.
Legacy caches without metadata are recomputed. An HMMER failure stops annotation
with a nonzero exit status; a partial table is never installed as a valid cache.
Fix the error and rerun to reuse any completed, valid searches. Use a fresh output
directory when comparing runs so older result files cannot be mistaken for new ones.

Gene-prediction FAA files have adjacent `.cache.json` completion records that bind
assembly content, caller/version, mode and both FAA/GFF contents. Files are prepared
before publication; prediction errors stop the run. Old generated genomes outside
the current input list are not collected.

Each successful UniOP work directory preserves `efesto-index.faa` and
`efesto-cache.json`. These bind input contents, input mode, script/interpreter and
prediction/index contents. The saved index, rather than a newly ordered FAA or an
arbitrary cached FAA, decodes the numeric prediction indices. Incomplete legacy
caches are recomputed. UniOP failures are optional-tool failures and are reported
as unavailable evidence rather than valid cached predictions.

## Coverage normalisation

Without `--norm_coverage`, the category matrix sums contig mean depth once per
retained gene hit (`mean_depth_sum`). With that flag it uses
`coverage_per_million_sum`:

```text
contig_value = mean_depth / sum(all supplied contig mean depths) × 1,000,000
category_value = sum(contig_value for each retained hit in the category)
```

Mean depth is already normalised for length. Equal mean depths therefore receive
equal values even if contig lengths differ. Normalisation is performed once per
genome/sample coverage input and does not require contig lengths. All-zero inputs
produce zero values. This is a DNA coverage measure, not transcript TPM.

Several retained genes on one contig contribute separately, and genes can share
categories, so the **category matrix is not constrained to sum to one million**.
Inputs normalised separately per MAG do not estimate that MAG's abundance in the
whole community. Keep the same coverage-input scope when making comparisons.
