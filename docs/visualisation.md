# Visualisation

MetalGenie-Evo includes `scripts/plot_heatmap.R` for publication-ready heatmaps.

---

## Features

- Hierarchical clustering (Ward's method) on both axes
- White → orange → red for gene counts; white → blue for coverage
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
Rscript scripts/plot_heatmap.R --input results/MetalGenie-Evo-heatmap-data.csv
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
- `MetalGenie-Evo-heatmap-data.csv` — gene-count matrix (categories × genomes)
- `MetalGenie-Evo-coverage-heatmap.csv` — coverage matrix (same shape; optional)

Both are in FeGenie's original CSV format and are compatible with FeGenie's
own R visualisation script.

---

## Custom downstream analysis

For custom plots, `MetalGenie-Evo-results-long.tsv` is the recommended starting
point (one row per ORF, all columns including `cluster_confidence` and `model_nseq`):

```r
library(tidyverse)

hits <- read_tsv("results/MetalGenie-Evo-results-long.tsv")

# Mean confidence by category
hits |>
  group_by(category) |>
  summarise(mean_conf = mean(cluster_confidence, na.rm=TRUE), n = n()) |>
  arrange(desc(mean_conf))

# Iron reduction genes with high confidence
hits |>
  filter(str_detect(category, "iron_reduction"), cluster_confidence >= 0.8) |>
  select(genome, contig, gene, bitscore, cluster_confidence)
```
