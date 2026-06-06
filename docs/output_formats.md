# Output Formats

MetalGenie-Evo writes results to `--out` (default: `metalgenie_evo_out/`).

---

## Standard outputs (always written)

| File | Description |
|------|-------------|
| `MetalGenie-Evo-summary.csv` | Per-ORF detailed results with cluster separators; includes `cluster_confidence` |
| `MetalGenie-Evo-geneSummary-clusters.csv` | FeGenie R-script compatible compact summary; includes `cluster_confidence` |
| `MetalGenie-Evo-heatmap-data.csv` | Gene-count matrix (categories × genomes) |
| `MetalGenie-Evo-results-long.tsv` | Tidy long-format TSV — one row per ORF, all scoring columns |
| `MetalGenie-Evo-results.gff3` | GFF3 feature file (written when `--gff_dir` or `--fna_dir` used) |
| `MetalGenie-Evo-summary-stats.tsv` | Run-level statistics (counts, confidence tiers, per-genome summaries) |

## Optional outputs

| File | Requires |
|------|----------|
| `MetalGenie-Evo-coverage-heatmap.csv` | `--bam` / `--bams` / `--depth` / `--depths` |
| `MetalGenie-Evo-OperonStructure.tsv` | `--operon_prediction` |
| `MetalGenie-Evo-anvio-functions.tsv` | `--anvio` |
| `MetalGenie-Evo-anvio-gene-scores.tsv` | `--anvio` |

---

## Column reference — `MetalGenie-Evo-results-long.tsv`

One row per ORF hit. The tidy format for downstream analysis.

| Column | Type | Description |
|--------|------|-------------|
| `category` | string | Functional category (e.g. `iron_reduction`) |
| `genome` | string | Source genome filename |
| `contig` | string | Source contig identifier |
| `orf` | string | ORF identifier (Bakta gene ID when `--bakta_gff_dir` used) |
| `gene` | string | Readable gene name (from `FeGenie-map.txt`) |
| `hmm_stem` | string | HMM file stem |
| `bitscore` | float | hmmsearch bitscore |
| `bitscore_cutoff` | float | Per-HMM calibrated cutoff (0 = fallback applied) |
| `confidence` | string | `calibrated` or `low_confidence` |
| `cluster_id` | int | Genomic cluster index (unique within genome) |
| `contig_len` | int | Contig length in bp |
| `start` | int | ORF start coordinate (bp, 1-based) |
| `end` | int | ORF end coordinate (bp, 1-based) |
| `strand` | string | `+` or `-` |
| `heme_c_motifs` | int | Count of CXXCH heme-binding motifs in sequence |
| `cluster_confidence` | float | Composite cluster reliability score (0–1.2) |
| `model_nseq` | int | Training sequences for this HMM (from registry) |
| `uniop_context` | string | UniOP operon ID, `singleton_<orf>`, or `not_in_operon` |

---

## Column reference — `MetalGenie-Evo-results.gff3`

Standard GFF3 format. Written automatically when coordinate data is available.
Loadable directly in IGV, Artemis, and genome browsers.

**Fixed columns (GFF3 spec):**

| Col | Description |
|-----|-------------|
| seqname | Contig identifier |
| source | `MetalGenie-Evo` |
| feature | `gene` |
| start / end | 1-based coordinates |
| score | hmmsearch bitscore |
| strand | `+` or `-` |
| frame | `.` |

**Attributes (semicolon-separated):**

| Attribute | Description |
|-----------|-------------|
| `ID` | ORF identifier |
| `gene` | Readable gene name |
| `category` | Functional category |
| `hmm_stem` | HMM file stem |
| `cluster_id` | Cluster index |
| `cluster_confidence` | Composite confidence score |
| `confidence` | `calibrated` or `low_confidence` |

---

## Column reference — `MetalGenie-Evo-summary-stats.tsv`

Two-column TSV (`section<TAB>key<TAB>value`). Sections:

| Section | Metrics |
|---------|---------|
| `RUN` | `total_orf_hits`, `total_clusters`, `genomes_with_hits`, `genomes_zero_hits`, `runtime_sec` |
| `CONFIDENCE` | `mean_cluster_confidence`, `high_confidence_clusters` (≥ 0.8), `medium_confidence_clusters` (0.5–0.8), `low_confidence_clusters` (< 0.5) |
| `CATEGORY` | ORF hit count per functional category (one row per category) |
| `GENOME` | Per-genome ORF hit count and mean cluster confidence (one row per genome) |

---

## Column reference — `MetalGenie-Evo-OperonStructure.tsv`

Written when `--operon_prediction` is used.

| Column | Description |
|--------|-------------|
| `genome` | Source genome |
| `orf` | ORF identifier |
| `operon_id` | UniOP operon assignment |
| `pair_probability` | Co-operon probability with adjacent gene |
| `gene` | Readable gene name |
| `category` | Functional category |

---

## Column reference — `MetalGenie-Evo-anvio-gene-scores.tsv`

Written when `--anvio` is used. Import with `anvi-import-misc-data`:

```bash
anvi-import-misc-data \
    -c CONTIGS.db \
    --target-data-table genes \
    results/MetalGenie-Evo-anvio-gene-scores.tsv
```

| Column | Description |
|--------|-------------|
| `gene_callers_id` | Anvi'o gene caller ID (Bakta ID if `--bakta_gff_dir` used) |
| `cluster_confidence` | Composite cluster reliability score |
| `co_occ_score` | Co-occurrence component of confidence |
| `hmm_weight` | HMM calibration weight component |
| `uniop_weight` | UniOP pair probability component |
| `bgc_boost` | BGC boost factor (1.0 or 1.2) |

---

## Cluster confidence tiers

| Score | Tier | Meaning |
|-------|------|---------|
| ≥ 0.8 | High | Complete operon, calibrated HMMs, compact genomic arrangement |
| 0.5 – 0.8 | Medium | Partial operon or low-nseq models involved |
| < 0.5 | Low | Single gene, uncalibrated model, or contig-edge hit |

`cluster_confidence` is capped at 1.2 (possible only when `--bgc_dir` triggers the 1.2× BGC boost on a near-perfect cluster).

---

## Working with outputs in R

FeGenie's original R visualisation script (`scripts/plot_heatmap.R`) is
compatible with `MetalGenie-Evo-heatmap-data.csv` and `MetalGenie-Evo-geneSummary-clusters.csv`.

For custom analysis, `MetalGenie-Evo-results-long.tsv` is the recommended
starting point:

```r
library(tidyverse)
hits <- read_tsv("results/MetalGenie-Evo-results-long.tsv")

# High-confidence iron reduction clusters
hits |>
  filter(category == "iron_reduction", cluster_confidence >= 0.8) |>
  group_by(genome, cluster_id) |>
  summarise(genes = paste(gene, collapse=","), confidence = first(cluster_confidence))
```
