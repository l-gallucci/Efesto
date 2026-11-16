# Output files

Files are written below `--out`. Outputs describe retained annotations, after
search and context filtering.

## Main files

| File | Contents / condition |
|---|---|
| `Efesto-summary.csv` | Detailed per-ORF results with cluster separators |
| `Efesto-geneSummary-clusters.csv` | Compact summary compatible with the FeGenie plotting format |
| `Efesto-heatmap-data.csv` | Gene-count matrix |
| `Efesto-results-long.tsv` | Tidy per-ORF table for analysis |
| `Efesto-summary-stats.tsv` | Counts, score tiers and run statistics |
| `Efesto-hits.faa` | Retained matched protein sequences |
| `Efesto-hits.gff3` | Retained CDS features with available coordinates |
| `Efesto-hits.fna` | Matched nucleotide sequences when source FNA and coordinates are available |
| `Efesto-run.log` | Run metadata and recorded tool information |

## Optional files

| File | Condition |
|---|---|
| `Efesto-coverage-heatmap.csv` | Coverage inputs supplied |
| `Efesto-OperonStructure.tsv` | Operon-prediction results available |
| `Efesto-anvio-functions.tsv` | `--anvio` |
| `Efesto-anvio-gene-scores.tsv` | `--anvio` |
| `_tblout_cache/` | Raw search tables; retain with `--keep_tblout` |

## Long-format table

The standard columns, in order, are:

```text
category  genome  contig  orf  gene  bitscore  bitscore_cutoff  confidence
cluster_id  heme_c_motifs  contig_len  cluster_confidence  model_nseq  topology_weight
```

`uniop_context` is appended when available. Inspect the actual header when reading
a result from a different revision.

| Field | Meaning |
|---|---|
| `category`, `gene` | Assigned category and readable family/gene label |
| `genome`, `contig`, `orf` | Input and gene identifiers; `orf` can carry a mapped Bakta label |
| `bitscore`, `bitscore_cutoff` | Full-sequence match score and configured model cutoff |
| `confidence` | `calibrated` when a positive model-specific cutoff exists; otherwise `low_confidence` |
| `cluster_id` | Genomic cluster assignment |
| `heme_c_motifs` | Count of configured heme-binding motif patterns |
| `contig_len` | Available contig length; zero can mean unavailable |
| `cluster_confidence` | Multiplicative heuristic score, capped at 1.2 |
| `model_nseq` | Source/model training-record count, not number of independently validated examples |
| `topology_weight` | Contribution of the applicable topology check |

Coordinates and the HMM stem are **not standard columns of this TSV**. Use the
GFF3 and retained search tables for those details.

## GFF3

The emitted feature type is `CDS`, with 1-based inclusive coordinates. The score
and phase columns are `.`. Attributes include `ID`, `gene`, `category`,
`hmm_stem`, `cluster_id`, `cluster_confidence` and `confidence`. Records without
coordinates are omitted from GFF3, even if present in the tabular results.

## Statistics and interpretation

The statistics table uses `section`, `metric`, `value` and `notes` columns. Confidence
tiers are ≥0.8, 0.5–<0.8 and <0.5. These are reporting bins; the normal run does
not apply a global 0.5 confidence filter. Apply any analysis threshold explicitly
and report it. See [interpretation](interpretation.md).

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

## Search cache and failures

With `--keep_tblout`, Efesto retains completed HMMER tables and adjacent
`.cache.json` files. Reuse requires matching protein/profile contents, score
threshold, HMMER executable and thread setting, plus a complete, unchanged table.
Legacy caches without metadata are recomputed. An HMMER failure stops annotation
with a nonzero exit status; a partial table is never installed as a valid cache.
Fix the error and rerun to reuse any completed, valid searches. Use a fresh output
directory when comparing runs so older result files cannot be mistaken for new ones.
