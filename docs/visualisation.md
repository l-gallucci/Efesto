# Visualisation

Efesto includes `scripts/plot_heatmap.R` for static and interactive heatmaps.

---

## Features

- Hierarchical clustering (Ward's method) on both axes
- White → orange → red for cluster counts; white → blue for coverage
- Count matrices with a maximum above 10 use `log1p` for display
- Static PDF/PNG via `pheatmap`
- Interactive self-contained HTML via `plotly`
- Compatible with FeGenie's R script input format

---

## Installation

```bash
conda install -c conda-forge r-pheatmap r-plotly r-htmlwidgets r-optparse r-rcolorbrewer r-scales
```

---

## Usage

```bash
# Process entire results directory (auto-detects CSV files)
Rscript scripts/plot_heatmap.R --input results/

# Static PDF only
Rscript scripts/plot_heatmap.R --input results/ --type static --format pdf --out figures/

# Interactive HTML, filter low-count categories
Rscript scripts/plot_heatmap.R --input results/ --type interactive --min_count 1

# Specify a single CSV
Rscript scripts/plot_heatmap.R --input results/Efesto-heatmap-data.csv
```

---

## Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--input` | *(required)* | CSV file or output directory |
| `--type` | `both` | `static` \| `interactive` \| `both` |
| `--format` | `both` | `pdf` \| `png` \| `both` (static only) |
| `--out` | same as input | Output directory |
| `--width` / `--height` | auto | Plot dimensions in inches |
| `--min_count` | `0` | Minimum value to include a category row |
| `--no_cluster_rows` | off | Disable row clustering |
| `--no_cluster_cols` | off | Disable column clustering |

---

## Input files

The script reads:

- `Efesto-heatmap-data.csv` — cluster-count matrix (categories × genomes)
- `Efesto-coverage-heatmap.csv` — coverage matrix (same shape; optional)

The matrices use the FeGenie-compatible CSV layout. Their units differ:

| Input | Unit |
|---|---|
| `Efesto-heatmap-data.csv` | Distinct retained clusters per category and genome |
| The same file generated with `--norm` | Clusters per 1,000 input proteins |
| `Efesto-coverage-heatmap.csv` | Sum of per-hit contig mean depths, or the configured coverage normalisation |

Two retained genes of the same category in one cluster contribute **one** to the
raw cluster matrix and **two** to a per-ORF count. Use `Efesto-results-long.tsv`
when the question is how many retained genes matched. State the unit and any
normalisation in figure captions. Coverage is not a confidence measure or
transcript TPM; see [coverage units](output_formats.md#coverage-normalisation).

---

## Custom downstream analysis

For custom plots, `Efesto-results-long.tsv` is the recommended starting
point (one row per ORF, all columns including `threshold_type`, evidence statuses and `model_nseq`):

```r
library(tidyverse)

hits <- read_tsv("results/Efesto-results-long.tsv")

# Retained family matches by category and search-gate type
hits |>
  count(category, threshold_type, name = "orf_hits")

# Inspect iron-reduction family matches alongside available evidence
hits |>
  filter(str_detect(category, "Iron-Reduction")) |>
  select(genome, contig, gene, bitscore, evalue, threshold_type,
         topology_verdict, heme_verdict, uniop_status)
```

Counts describe retained annotations; neither an optional-tool status nor a
model-specific cutoff alone establishes biological function. See
[interpretation](interpretation.md) and [separate evidence](evidence.md).
