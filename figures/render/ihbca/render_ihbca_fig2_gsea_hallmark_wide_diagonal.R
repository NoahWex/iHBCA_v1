#!/usr/bin/env Rscript
# Render GSEA Hallmark heatmap per track: cell type × pathway, color = NES,
# asterisks for padj-significant.
#
# Aesthetic anchors (Fig 2):
#   - Color scale: expression_diverging (red high, white mid, blue low)
#   - 5-7pt typography per panel_theme
#   - Asterisk significance markers in cell text
#
# Per-track invocation:
#   --gsea-csv   outputs/gsea/{track}/pathways_long.csv
#   --out-pdf    outputs/plots/gsea_heatmap_{track}.pdf
#   --top-n      number of pathways to show per cell type union (default 25)
#   --fdr-mark   FDR threshold for asterisk (default 0.10)
#
# Code pattern provenance:
#   - Heatmap idiom: ggplot + geom_tile, common in Fig 2 supp
#   - Color scale: render_fig2_beeswarm_faceted.R:274-278

suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

parser <- ArgumentParser()
parser$add_argument("--gsea-csv", required = TRUE,
                    help = "pathways_long.csv from script 14")
parser$add_argument("--out-pdf",  required = TRUE)
parser$add_argument("--publication-root", required = TRUE)
parser$add_argument("--cell-types", default = "",
                    help = "Subset (semicolons or commas). Default: all.")
parser$add_argument("--cell-types-file", default = "",
                    help = "Optional file with one cell type per line (overrides --cell-types)")
parser$add_argument("--top-n", type = "integer", default = 25,
                    help = "Top N pathways by min padj across cells (ignored if --pathways-file set)")
parser$add_argument("--pathways-file", default = "",
                    help = "Optional file with one pathway name per line (overrides --top-n; pathways must be in 'collection::name' or 'name' form)")
parser$add_argument("--fdr-mark", type = "double", default = 0.10)
parser$add_argument("--orient", default = "tall",
                    choices = c("tall", "wide"),
                    help = "tall: pathways on Y (default). wide: pathways on X, cells on Y.")
parser$add_argument("--sort", default = "mean_nes",
                    choices = c("mean_nes", "diagonal"),
                    help = "Pathway ordering. diagonal: group by argmax-cell, sort by NES magnitude within group.")
parser$add_argument("--cell-order", default = "",
                    help = "Optional comma-separated cell type order (default: alphabetical / data order)")
parser$add_argument("--pdf-width",  type = "double", default = NULL)
parser$add_argument("--pdf-height", type = "double", default = NULL)
args <- parser$parse_args()

source(file.path(args$publication_root, "publication", "config",
                 "load_aesthetics.R"))
aes_cfg <- load_aesthetics(file.path(args$publication_root, "publication",
                                     "config"))
panel_theme <- get_theme(aes_cfg)
div_scale   <- aes_cfg$scales$expression_diverging

df <- read_csv(args$gsea_csv, show_col_types = FALSE)
cat(sprintf("[1] loaded %d (cell × pathway) rows\n", nrow(df)))

wanted <- character(0)
if (nchar(args$cell_types_file) > 0L && file.exists(args$cell_types_file)) {
  wanted <- readLines(args$cell_types_file) |> trimws()
  wanted <- wanted[nzchar(wanted)]
} else if (nchar(args$cell_types) > 0L) {
  # Accept either commas or semicolons
  wanted <- strsplit(args$cell_types, "[;,]")[[1]] |> trimws()
}
if (length(wanted) > 0L) {
  df <- df |> filter(cell_type %in% wanted)
  cat(sprintf("  subset to %d cell types: %s\n",
              length(wanted), paste(wanted, collapse = ", ")))
}

# Hallmark-only filter (defensive — main use case)
df <- df |> filter(collection == "H")

# Pathway selection: explicit list (overrides) > top-N
if (nchar(args$pathways_file) > 0L && file.exists(args$pathways_file)) {
  wanted_pw <- readLines(args$pathways_file) |> trimws()
  wanted_pw <- wanted_pw[nzchar(wanted_pw)]
  # pathway_name in df has "HALLMARK_" prefix; accept either form
  top_pw <- intersect(unique(df$pathway_name), wanted_pw)
  cat(sprintf("[2] explicit pathway list: %d requested, %d matched\n",
              length(wanted_pw), length(top_pw)))
} else {
  top_pw <- df |>
    group_by(pathway_name) |>
    summarise(min_padj = min(padj, na.rm = TRUE), .groups = "drop") |>
    arrange(min_padj) |>
    slice_head(n = args$top_n) |>
    pull(pathway_name)
  cat(sprintf("[2] keeping top %d pathways by min padj\n", length(top_pw)))
}

df_sub <- df |> filter(pathway_name %in% top_pw)

# Cell-type order (used for diagonal sort grouping AND axis order)
cell_order <- if (nchar(args$cell_order) > 0L) {
  trimws(strsplit(args$cell_order, ",", fixed = TRUE)[[1]])
} else {
  sort(unique(df_sub$cell_type))
}
miss_cells <- setdiff(cell_order, unique(df_sub$cell_type))
if (length(miss_cells)) {
  stop("--cell-order entries not found in data: ",
       paste(miss_cells, collapse = ", "),
       "\nData cell_type values present: ",
       paste(sort(unique(df_sub$cell_type)), collapse = ", "))
}
df_sub$cell_type <- factor(df_sub$cell_type, levels = cell_order)
# Display cleanup: strip "str__" / "epi__" / "imm__" prefix; underscores → hyphens
clean_label <- function(x) {
  gsub("_", "-",
       sub("^(str|imm|epi)__", "", as.character(x)))
}
df_sub$cell_type_display <- factor(clean_label(df_sub$cell_type),
                                    levels = clean_label(cell_order))

# Pathway ordering
pw_order <- if (args$sort == "diagonal") {
  # Sign-pattern grouped diagonal sort:
  #   1. Pass 1: rough soft-argmax centroid sort of cells (to find a sensible row order)
  #   2. With final cell order, compute each pathway's sign pattern (+/− vector across L2s)
  #   3. Group pathways by exact sign pattern
  #   4. Order groups by median sign_centroid (positive-in-early-L2 patterns first)
  #   5. Within each group, sort pathways by mean |NES| descending (most intense first)
  BETA <- 2.0
  nes_z <- df_sub |>
    group_by(cell_type) |>
    mutate(nes_z = (NES - mean(NES, na.rm = TRUE)) / sd(NES, na.rm = TRUE)) |>
    ungroup()

  pathway_soft_pos <- function(z_df, ordr) {
    z_df |>
      group_by(pathway_name) |>
      summarise(
        soft_pos = {
          l2idx <- match(as.character(cell_type), ordr)
          ord   <- order(l2idx)
          w     <- exp(BETA * nes_z[ord])
          sum(l2idx[ord] * w) / sum(w)
        },
        .groups = "drop")
  }

  # Pass 1: pathway centroids → cell centroids → reorder rows
  pw_pos1   <- pathway_soft_pos(nes_z, cell_order)
  pw_lookup <- setNames(pw_pos1$soft_pos, pw_pos1$pathway_name)
  cell_centroids <- nes_z |>
    mutate(pw_soft_pos = pw_lookup[as.character(pathway_name)],
           w = pmax(nes_z, 0) + 1e-3) |>
    group_by(cell_type) |>
    summarise(cell_centroid = sum(pw_soft_pos * w) / sum(w),
              .groups = "drop") |>
    arrange(cell_centroid)
  new_cell_order <- as.character(cell_centroids$cell_type)
  if (!identical(new_cell_order, cell_order)) {
    cat(sprintf("[2c] re-ordering cells by centroid: %s -> %s\n",
                paste(cell_order, collapse = " → "),
                paste(new_cell_order, collapse = " → ")))
    cell_order <- new_cell_order
    df_sub$cell_type <- factor(df_sub$cell_type, levels = cell_order)
    df_sub$cell_type_display <- factor(clean_label(df_sub$cell_type),
                                        levels = clean_label(cell_order))
  } else {
    cat("[2c] cell order already optimal\n")
  }

  # Pass 2: sign-pattern grouping under final cell order
  df_for_pat <- df_sub |>
    mutate(l2_idx = match(as.character(cell_type), cell_order)) |>
    arrange(pathway_name, l2_idx)
  sign_pat <- df_for_pat |>
    group_by(pathway_name) |>
    summarise(
      pattern        = paste0(ifelse(NES > 0, "+", "-"), collapse = ""),
      mean_abs_nes   = mean(abs(NES), na.rm = TRUE),
      sign_centroid  = sum(l2_idx * (pmax(NES, 0) + 1e-3)) /
                       sum(pmax(NES, 0) + 1e-3),
      .groups = "drop")
  group_meta <- sign_pat |>
    group_by(pattern) |>
    summarise(group_centroid = median(sign_centroid),
              n_pathways = dplyr::n(),
              .groups = "drop") |>
    arrange(group_centroid)
  sign_pat <- sign_pat |>
    left_join(group_meta, by = "pattern") |>
    arrange(group_centroid, desc(mean_abs_nes))
  cat("[2b] diagonal sort (sign-pattern groups + intensity within):\n")
  print(as.data.frame(group_meta))
  cat("\nFinal pathway order:\n")
  print(as.data.frame(sign_pat[, c("pathway_name", "pattern", "mean_abs_nes")]))
  sign_pat |> pull(pathway_name)
} else {
  df_sub |>
    group_by(pathway_name) |>
    summarise(mean_nes = mean(NES, na.rm = TRUE), .groups = "drop") |>
    arrange(mean_nes) |>
    pull(pathway_name)
}
df_sub$pathway_name <- factor(df_sub$pathway_name, levels = pw_order)

# Strip "HALLMARK_" prefix for display
df_sub$pathway_display <- gsub("^HALLMARK_", "", as.character(df_sub$pathway_name))
df_sub$pathway_display <- factor(df_sub$pathway_display,
                                  levels = gsub("^HALLMARK_", "", pw_order))

# Significance marker
df_sub$sig_mark <- ifelse(is.na(df_sub$padj), "",
                           ifelse(df_sub$padj < 0.01, "**",
                                   ifelse(df_sub$padj < args$fdr_mark, "*", "")))

# Symmetric NES cap
nes_cap <- max(quantile(abs(df_sub$NES), 0.98, na.rm = TRUE), 1.0)

if (args$orient == "tall") {
  p <- ggplot(df_sub, aes(x = cell_type_display, y = pathway_display, fill = NES))
  x_label_angle <- 45; x_label_hjust <- 1; x_label_size <- 6
  y_label_size <- 5.5
} else {
  p <- ggplot(df_sub, aes(x = pathway_display, y = cell_type_display, fill = NES))
  x_label_angle <- 45; x_label_hjust <- 1; x_label_size <- 6
  y_label_size <- 7  # fewer y labels (4 L2s) → larger font OK
}

p <- p +
  geom_tile(color = "white", linewidth = 0.2) +
  geom_text(aes(label = sig_mark), size = 2.4, color = "black",
            fontface = "bold", vjust = 0.7) +
  scale_fill_gradient2(
    low = div_scale$low, mid = "grey95", high = div_scale$high,
    midpoint = 0, limits = c(-nes_cap, nes_cap),
    name = "NES", oob = scales::squish,
    guide = guide_colorbar(
      barheight = unit(if (args$orient == "wide") 0.65 else 1.2, "in"),
      barwidth  = unit(0.10, "in"),
      title.position = "top",
      title.hjust = 0.5,
      frame.colour = "grey50",
      frame.linewidth = 0.2,
      ticks.colour = "grey50",
      ticks.linewidth = 0.2)) +
  labs(x = NULL, y = NULL) +
  panel_theme +
  theme(axis.text.x = element_text(angle = x_label_angle, hjust = x_label_hjust,
                                    size = x_label_size),
        axis.text.y = element_text(size = y_label_size),
        panel.grid  = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 6),
        legend.text  = element_text(size = 5),
        legend.margin = margin(0, 0, 0, 2),
        legend.box.margin = margin(0, 0, 0, 0),
        plot.margin   = margin(2, 4, 2, 2))

# Dimensions: wide orient → wide+short; tall orient → tall+narrow
n_cells <- length(unique(df_sub$cell_type))
n_pw    <- length(top_pw)
if (args$orient == "tall") {
  w <- args$pdf_width  %||% max(2.5, 0.32 * n_cells + 1.5)
  h <- args$pdf_height %||% max(2.5, 0.18 * n_pw + 1.0)
} else {
  w <- args$pdf_width  %||% max(4.0, 0.35 * n_pw + 1.7)
  h <- args$pdf_height %||% max(1.5, 0.35 * n_cells + 0.9)
}

dir.create(dirname(args$out_pdf), showWarnings = FALSE, recursive = TRUE)
ggsave(args$out_pdf, p, width = w, height = h, units = "in",
       device = cairo_pdf)
cat(sprintf("[3] wrote %s (%.1f x %.1f in)\n", args$out_pdf, w, h))
