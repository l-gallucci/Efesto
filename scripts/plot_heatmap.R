#!/usr/bin/env Rscript
# =============================================================================
# plot_heatmap.R  —  Visualisation script for Efesto output
# =============================================================================
# Generates hierarchically clustered, publication-ready heatmaps from Efesto
# CSV output. Automatically detects the input type (gene counts vs coverage)
# and adapts colour scales accordingly.
#
# Usage
# -----
#   Rscript scripts/plot_heatmap.R --input results/Efesto-heatmap-data.csv
#   Rscript scripts/plot_heatmap.R --input results/ --type both --format both
#   Rscript scripts/plot_heatmap.R --input results/ --size a4 --orientation landscape
#
# Arguments
# ---------
#   --input       Path to a single CSV file OR a Efesto output directory.
#                 If a directory is given, all *heatmap*.csv files are processed.
#   --type        Plot type: "static" | "interactive" | "both"  (default: both)
#   --format      Output format for static: "pdf" | "png" | "both"  (default: both)
#   --out         Output directory for plots (default: same dir as input CSV)
#   --size        Page-size preset: "full" | "a4" | "half" | "auto"  (default: full)
#                   full = US Letter, 8.5x11 in   a4 = 8.27x11.69 in
#                   half = 8.5x5.5 in             auto = size-to-content (legacy)
#   --orientation Page orientation for full/a4/half presets: "portrait" | "landscape"
#                 (default: portrait)
#   --width       Plot width in inches  (overrides --size)
#   --height      Plot height in inches (overrides --size)
#   --font        Font family for all plot text (default: Helvetica)
#   --min_count   Minimum gene count/coverage to include a category row (default: 0)
#   --no_command  Suppress the inset caption showing the Efesto command that
#                 produced the input CSV (auto-read from Efesto-run.log)
#   --run_log     Explicit path to an Efesto-run.log to pull the command from
#                 (default: auto-detected next to the input CSV)
#   --help        Show this help message
#
# Required R packages
# -------------------
#   ggplot2, pheatmap, plotly, htmlwidgets, optparse, RColorBrewer, scales
#
#   Install with:
#     conda install -c conda-forge r-ggplot2 r-pheatmap r-plotly r-htmlwidgets
#                                  r-optparse r-rcolorbrewer r-scales
#   or:
#     Rscript -e 'install.packages(c("ggplot2","pheatmap","plotly","htmlwidgets",
#                                    "optparse","RColorBrewer","scales"))'
#
# Font note
# ---------
#   "Helvetica" is one of R's built-in PDF base-14 fonts, so --format pdf always
#   renders correctly. For --format png, R resolves "Helvetica" through the
#   system's font backend: on macOS (Quartz) it just works; on Linux it is
#   normally aliased to Nimbus Sans by fontconfig (packages gsfonts/urw-fonts).
#   If it's missing entirely, PNG rendering silently falls back to the device
#   default rather than erroring. Use --font to pick a different family
#   (e.g. "Arial") if needed.
# =============================================================================

suppressPackageStartupMessages({
  if (!require("optparse",     quietly = TRUE)) stop("Package 'optparse' required.")
  if (!require("ggplot2",      quietly = TRUE)) stop("Package 'ggplot2' required.")
  if (!require("pheatmap",     quietly = TRUE)) stop("Package 'pheatmap' required.")
  if (!require("plotly",       quietly = TRUE)) stop("Package 'plotly' required.")
  if (!require("htmlwidgets",  quietly = TRUE)) stop("Package 'htmlwidgets' required.")
  if (!require("RColorBrewer", quietly = TRUE)) stop("Package 'RColorBrewer' required.")
  if (!require("scales",       quietly = TRUE)) stop("Package 'scales' required.")
  if (!require("grid",         quietly = TRUE)) stop("Package 'grid' required.")
})

# ─────────────────────────────────────────────────────────────────────────────
# Argument parsing
# ─────────────────────────────────────────────────────────────────────────────

option_list <- list(
  make_option("--input",     type = "character", default = NULL,
              help = "CSV file or Efesto output directory"),
  make_option("--type",      type = "character", default = "both",
              help = "Plot type: static | interactive | both  [default: both]"),
  make_option("--format",    type = "character", default = "both",
              help = "Static format: pdf | png | both  [default: both]"),
  make_option("--out",       type = "character", default = NULL,
              help = "Output directory for plots [default: same as input]"),
  make_option("--size",      type = "character", default = "full",
              help = "Page-size preset: full | a4 | half | auto  [default: full]"),
  make_option("--orientation", type = "character", default = "portrait",
              help = "Page orientation for full/a4/half: portrait | landscape  [default: portrait]"),
  make_option("--width",     type = "double",    default = NULL,
              help = "Plot width in inches [default: from --size]"),
  make_option("--height",    type = "double",    default = NULL,
              help = "Plot height in inches [default: from --size]"),
  make_option("--font",      type = "character", default = "Helvetica",
              help = "Font family for plot text [default: Helvetica]"),
  make_option("--min_count", type = "double",    default = 0,
              help = "Min gene count/coverage per category to include [default: 0]"),
  make_option("--no_cluster_rows", action = "store_true", default = FALSE,
              help = "Disable hierarchical clustering of categories"),
  make_option("--no_cluster_cols", action = "store_true", default = FALSE,
              help = "Disable hierarchical clustering of genomes"),
  make_option("--no_command", action = "store_true", default = FALSE,
              help = "Suppress inset caption showing the Efesto command"),
  make_option("--run_log",   type = "character", default = NULL,
              help = "Path to Efesto-run.log to read the command from [default: auto-detect]")
)

parser <- OptionParser(
  usage       = "Rscript scripts/plot_heatmap.R --input <file_or_dir> [options]",
  option_list = option_list
)
opt <- parse_args(parser)

if (is.null(opt$input)) {
  print_help(parser)
  stop("--input is required.")
}

# ─────────────────────────────────────────────────────────────────────────────
# Page-size presets (inches)
# ─────────────────────────────────────────────────────────────────────────────

resolve_page_size <- function(size, orientation, width, height) {
  presets <- list(
    full = c(8.5,  11),
    a4   = c(8.27, 11.69),
    half = c(8.5,  5.5)
  )
  if (identical(size, "auto")) {
    # No preset — fall back to content-based sizing in make_static() unless
    # the user gave an explicit --width/--height.
    return(list(width = width, height = height))
  } else {
    if (!size %in% names(presets))
      stop("--size must be one of: full, a4, half, auto")
    dims <- presets[[size]]
    if (identical(orientation, "landscape")) dims <- rev(dims)
    else if (!identical(orientation, "portrait"))
      stop("--orientation must be 'portrait' or 'landscape'")
  }
  list(width = if (!is.null(width)) width else dims[1],
       height = if (!is.null(height)) height else dims[2])
}

page <- resolve_page_size(opt$size, opt$orientation, opt$width, opt$height)

# ─────────────────────────────────────────────────────────────────────────────
# Efesto command caption (inset)
# ─────────────────────────────────────────────────────────────────────────────

find_run_log <- function(csv_path, explicit) {
  if (!is.null(explicit)) {
    if (file.exists(explicit)) return(explicit)
    warning("--run_log not found: ", explicit)
    return(NULL)
  }
  candidate <- file.path(dirname(csv_path), "Efesto-run.log")
  if (file.exists(candidate)) return(candidate)
  NULL
}

read_command_caption <- function(csv_path, explicit, wrap_width = 100) {
  log_path <- find_run_log(csv_path, explicit)
  if (is.null(log_path)) return(NULL)
  lines <- readLines(log_path, warn = FALSE)
  cmd_line <- grep("^Command\\s*:", lines, value = TRUE)
  if (length(cmd_line) == 0) return(NULL)
  cmd <- sub("^Command\\s*:\\s*", "", cmd_line[1])
  paste(strwrap(cmd, width = wrap_width), collapse = "\n")
}

# ─────────────────────────────────────────────────────────────────────────────
# Collect input CSV files
# ─────────────────────────────────────────────────────────────────────────────

collect_csvs <- function(path) {
  if (file.info(path)$isdir) {
    files <- list.files(path, pattern = "heatmap.*\\.csv$",
                        full.names = TRUE, ignore.case = TRUE)
    # Also grab the comment-header coverage file
    cov <- list.files(path, pattern = "coverage-heatmap.*\\.csv$",
                      full.names = TRUE, ignore.case = TRUE)
    all_files <- unique(c(files, cov))
    if (length(all_files) == 0)
      stop("No *heatmap*.csv files found in: ", path)
    return(all_files)
  }
  if (!file.exists(path)) stop("File not found: ", path)
  return(path)
}

csv_files <- collect_csvs(opt$input)
cat(sprintf("[INFO] Found %d heatmap file(s)\n", length(csv_files)))

# ─────────────────────────────────────────────────────────────────────────────
# Data loading and preprocessing
# ─────────────────────────────────────────────────────────────────────────────

detect_type <- function(filepath) {
  # Peek at comment lines for coverage metric label
  lines <- readLines(filepath, n = 3)
  if (any(grepl("^#.*coverage metric", lines, ignore.case = TRUE)))
    return("coverage")
  if (any(grepl("^#", lines)))
    return("coverage")   # any comment header → coverage file
  return("count")
}

load_heatmap_csv <- function(filepath) {
  # Skip comment lines starting with #
  lines <- readLines(filepath)
  data_lines <- lines[!grepl("^#", lines)]
  con <- textConnection(data_lines)
  on.exit(close(con))
  df <- read.csv(con, check.names = FALSE, row.names = 1)
  df
}

clean_matrix <- function(df, min_count = 0) {
  # Remove rows where all values <= min_count
  keep <- apply(df, 1, function(x) max(x, na.rm = TRUE) > min_count)
  df <- df[keep, , drop = FALSE]
  if (nrow(df) == 0) stop("No rows remain after filtering (--min_count too high?)")
  # Remove columns that are all zero
  keep_col <- apply(df, 2, function(x) max(x, na.rm = TRUE) > 0)
  df[, keep_col, drop = FALSE]
}

pretty_category <- function(x) {
  x <- gsub("_", " ", x)
  x <- gsub("-", " — ", x)
  x <- tools::toTitleCase(x)
  x
}

# ─────────────────────────────────────────────────────────────────────────────
# Colour palettes
# ─────────────────────────────────────────────────────────────────────────────

palette_for_type <- function(type) {
  if (type == "coverage") {
    # Blue-purple gradient for coverage/TPM
    colorRampPalette(c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"))(100)
  } else {
    # White → warm orange → dark red for gene counts
    colorRampPalette(c("#FFFFFF", "#FFF3E0", "#FF8F00", "#BF360C"))(100)
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# Static heatmap (pheatmap)
# ─────────────────────────────────────────────────────────────────────────────

make_static <- function(mat, title, type, out_base, caption,
                        formats = "both",
                        width = NULL, height = NULL,
                        font = "Helvetica",
                        cluster_rows = TRUE, cluster_cols = TRUE) {

  pal   <- palette_for_type(type)
  nrow_ <- nrow(mat)
  ncol_ <- ncol(mat)

  # Auto-size (used only when --size auto and no explicit --width/--height)
  w <- if (!is.null(width))  width  else max(8,  2 + ncol_ * 0.5)
  h <- if (!is.null(height)) height else max(6,  2 + nrow_ * 0.35)

  # Row labels — pretty-print category names
  rownames(mat) <- pretty_category(rownames(mat))

  # Legend label
  legend_title <- if (type == "coverage") "TPM / Depth" else "Gene count"

  # Value breaks: log1p scaling for counts, linear for coverage
  if (type == "count" && max(mat) > 10) {
    mat_display <- log1p(mat)
    legend_title <- "log(count + 1)"
  } else {
    mat_display <- mat
  }

  # Build the gtable once, without drawing, so we can render it ourselves on a
  # manually-opened device — this is what lets us control font/size/DPI and
  # add the command-caption inset underneath the heatmap.
  built <- pheatmap::pheatmap(
    mat_display,
    color             = pal,
    cluster_rows      = cluster_rows && nrow_ > 1,
    cluster_cols      = cluster_cols && ncol_ > 1,
    clustering_method = "ward.D2",
    border_color      = "white",
    cellwidth         = max(12, min(30, 400 / ncol_)),
    cellheight        = max(10, min(20, 300 / nrow_)),
    fontsize_row      = 9,
    fontsize_col      = 8,
    angle_col         = 45,
    main              = title,
    legend            = TRUE,
    fontfamily        = font,
    filename          = NA,
    silent            = TRUE
  )

  render <- function(open_device) {
    open_device()
    grid::grid.newpage()
    grid::grid.draw(built$gtable)
    if (!is.null(caption)) {
      grid::grid.text(
        caption,
        x  = unit(0.5, "npc"),
        y  = unit(4, "bigpts"),
        just = "bottom",
        gp = grid::gpar(fontsize = 6.5, fontfamily = font, col = "grey35",
                         lineheight = 1.1)
      )
    }
    grDevices::dev.off()
  }

  if (formats %in% c("pdf", "both")) {
    f <- paste0(out_base, ".pdf")
    render(function() grDevices::pdf(f, width = w, height = h, family = font))
    cat(sprintf("[INFO]   Static PDF  → %s  (%.2fx%.2f in)\n", f, w, h))
  }
  if (formats %in% c("png", "both")) {
    f <- paste0(out_base, ".png")
    render(function() grDevices::png(f, width = w, height = h, units = "in",
                                     res = 300, family = font))
    cat(sprintf("[INFO]   Static PNG  → %s  (%.2fx%.2f in @300dpi)\n", f, w, h))
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# Interactive heatmap (plotly)
# ─────────────────────────────────────────────────────────────────────────────

make_interactive <- function(mat, title, type, out_base, caption,
                             font = "Helvetica",
                             cluster_rows = TRUE, cluster_cols = TRUE) {

  # Apply clustering to get dendrogram order
  if (cluster_rows && nrow(mat) > 1) {
    hc_row <- hclust(dist(mat),    method = "ward.D2")
    mat    <- mat[hc_row$order, , drop = FALSE]
  }
  if (cluster_cols && ncol(mat) > 1) {
    hc_col <- hclust(dist(t(mat)), method = "ward.D2")
    mat    <- mat[, hc_col$order,  drop = FALSE]
  }

  # Display values
  if (type == "count" && max(mat) > 10) {
    display_mat <- log1p(mat)
    cbar_title  <- "log(n + 1)"
  } else {
    display_mat <- mat
    cbar_title  <- if (type == "coverage") "TPM / Depth" else "Gene count"
  }

  colorscale <- if (type == "coverage") {
    list(c(0, "#F7FBFF"), c(0.25, "#C6DBEF"), c(0.5, "#6BAED6"),
         c(0.75, "#2171B5"), c(1, "#08306B"))
  } else {
    list(c(0, "#FFFFFF"), c(0.33, "#FFF3E0"), c(0.66, "#FF8F00"),
         c(1, "#BF360C"))
  }

  # Hover text — show original (un-transformed) values
  hover_text <- matrix(
    sprintf("<b>%s</b><br>%s<br>Value: %s",
            rep(pretty_category(rownames(mat)), ncol(mat)),
            rep(colnames(mat), each = nrow(mat)),
            format(as.vector(mat), digits = 3, big.mark = ",")),
    nrow = nrow(mat), ncol = ncol(mat)
  )

  annotations <- NULL
  if (!is.null(caption)) {
    annotations <- list(list(
      text = gsub("\n", "<br>", caption),
      x = 0, y = -0.02, xref = "paper", yref = "paper",
      xanchor = "left", yanchor = "top",
      showarrow = FALSE,
      font = list(size = 9, family = font, color = "grey35")
    ))
  }

  fig <- plotly::plot_ly(
    z           = display_mat,
    x           = colnames(mat),
    y           = pretty_category(rownames(mat)),
    type        = "heatmap",
    colorscale  = colorscale,
    text        = hover_text,
    hoverinfo   = "text",
    colorbar    = list(title = cbar_title)
  ) %>%
    plotly::layout(
      title       = list(text = title, font = list(size = 14, family = font)),
      font        = list(family = font),
      xaxis       = list(title = "", tickangle = -45,
                         tickfont = list(size = 10, family = font)),
      yaxis       = list(title = "", autorange = "reversed",
                         tickfont = list(size = 9, family = font)),
      annotations = annotations,
      margin      = list(l = 200, b = 140, t = 60, r = 60)
    ) %>%
    plotly::config(displayModeBar = TRUE,
                   toImageButtonOptions = list(
                     format = "png", filename = basename(out_base),
                     width = 1200, height = 800
                   ))

  out_file <- paste0(out_base, ".html")
  htmlwidgets::saveWidget(fig, file = out_file, selfcontained = TRUE)
  cat(sprintf("[INFO]   Interactive → %s\n", out_file))
}

# ─────────────────────────────────────────────────────────────────────────────
# Main loop
# ─────────────────────────────────────────────────────────────────────────────

for (csv_path in csv_files) {

  cat(sprintf("\n[INFO] Processing: %s\n", basename(csv_path)))

  # Detect data type
  data_type <- detect_type(csv_path)
  cat(sprintf("[INFO]   Detected type: %s\n", data_type))

  # Load and clean
  df  <- tryCatch(load_heatmap_csv(csv_path),
                  error = function(e) { warning("Could not load ", csv_path, ": ", e$message); NULL })
  if (is.null(df)) next

  mat <- tryCatch(clean_matrix(as.matrix(df), min_count = opt$min_count),
                  error = function(e) { warning(e$message); NULL })
  if (is.null(mat)) next

  cat(sprintf("[INFO]   Matrix: %d categories × %d genomes\n", nrow(mat), ncol(mat)))

  # Determine output directory and base name
  out_dir <- if (!is.null(opt$out)) opt$out else dirname(csv_path)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  stem     <- tools::file_path_sans_ext(basename(csv_path))
  out_base <- file.path(out_dir, stem)

  # Title
  title <- if (data_type == "coverage") {
    "Efesto — Coverage-based abundance"
  } else {
    "Efesto — Gene presence by category"
  }

  # Inset caption: the Efesto command that produced this CSV
  caption <- if (opt$no_command) NULL else read_command_caption(csv_path, opt$run_log)
  if (!opt$no_command && is.null(caption))
    cat("[INFO]   No Efesto-run.log found — skipping command caption\n")

  cluster_rows <- !opt$no_cluster_rows
  cluster_cols <- !opt$no_cluster_cols

  # Generate plots
  if (opt$type %in% c("static", "both")) {
    tryCatch(
      make_static(mat, title, data_type, out_base, caption,
                  formats      = opt$format,
                  width        = page$width,
                  height       = page$height,
                  font         = opt$font,
                  cluster_rows = cluster_rows,
                  cluster_cols = cluster_cols),
      error = function(e) warning("Static plot failed: ", e$message)
    )
  }

  if (opt$type %in% c("interactive", "both")) {
    tryCatch(
      make_interactive(mat, title, data_type, out_base, caption,
                       font         = opt$font,
                       cluster_rows = cluster_rows,
                       cluster_cols = cluster_cols),
      error = function(e) warning("Interactive plot failed: ", e$message)
    )
  }
}

cat("\n[DONE] All plots written.\n")
