#!/usr/bin/env Rscript
# Generic Fig 2-aesthetic volcano renderer with red-white-blue diverging by
# sign(logFC) × −log10(padj). Replaces the Okabe-Ito CAF-class coloring in
# 06_volcano_fibros_caf_labels.R; preserves the CAF-marker labeling overlay
# (now plain text labels, no class-fill since color encodes significance).
#
# Per-track invocation:
#   --de-dir  outputs/de_results/{cohort}/{formula}      # one CSV per cell type
#   --out-dir outputs/plots/volcanos_diverging_{tag}
#   --cell-types "str__Fibro_major,str__Fibro_IGF1,..."  # subset to render
#
# Aesthetic anchors:
#   - Color: scale_color_gradient2 low=expression_diverging.low (#2166AC),
#     mid="grey90", high=expression_diverging.high (#B2182B)
#   - Dashed FDR + LFC reference lines
#   - 6pt CAF marker labels via ggrepel, only on significant marker points
#   - No on-plot title (legends go in figure caption / index.yaml)
#
# Code pattern provenance:
#   - Color scale + cap: render_fig2_beeswarm_faceted.R:151-153, 274-278
#   - CAF marker panel: 06_volcano_fibros_caf_labels.R:38-67
#   - Volcano layer + ggrepel: 06_volcano_fibros_caf_labels.R volcano-build block

suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

parser <- ArgumentParser()
parser$add_argument("--de-dir", required = TRUE)
parser$add_argument("--out-dir", required = TRUE)
parser$add_argument("--publication-root", required = TRUE)
parser$add_argument("--cell-types", default = "",
                    help = "Comma-separated subset (default: all non-skipped CSVs in de-dir)")
parser$add_argument("--fdr", type = "double", default = 0.05)
parser$add_argument("--lfc", type = "double", default = 0.5)
parser$add_argument("--fdr-floor", type = "double", default = 1e-15)
args <- parser$parse_args()

source(file.path(args$publication_root, "publication", "config",
                 "load_aesthetics.R"))
aes_cfg <- load_aesthetics(file.path(args$publication_root, "publication",
                                     "config"))
panel_theme <- get_theme(aes_cfg)
div_scale   <- aes_cfg$scales$expression_diverging

dir.create(args$out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- CAF marker panel (from 06) ----------------------------------------
ICAF   <- c("IL6","CXCL1","CXCL2","CXCL3","CXCL8","CXCL12","CXCL14",
            "CCL2","CCL7","LIF","PTGS2","C3","HAS1","HAS2","MMP1","MMP3",
            "TIMP1","SERPINE1","PDPN","DPT","DPP4","ICAM1")
MYCAF  <- c("ACTA2","TAGLN","MYH11","MCAM","RGS5","CNN1")
MATCAF <- c("POSTN","COMP","COL1A1","COL3A1","COL10A1","COL11A1","COL12A1",
            "FN1","ASPN","VCAN","LRRC15")
PANCAF <- c("FAP","THY1","S100A4","PDGFRA","PDGFRB")
OTHER  <- c("LOX","LOXL2","MMP2","MMP11","CCN2","CTGF","SFRP2","SFRP4",
            "GREM1","INHBA","MMP10","MMP12")
all_markers <- unique(c(ICAF, MYCAF, MATCAF, PANCAF, OTHER))

# ---- list DE CSVs ------------------------------------------------------
de_csvs <- list.files(args$de_dir, pattern = "\\.csv$", full.names = TRUE)
de_csvs <- de_csvs[!grepl("\\.skipped\\.csv$", de_csvs)]
if (nchar(args$cell_types) > 0L) {
  wanted <- strsplit(args$cell_types, ",", fixed = TRUE)[[1]] |> trimws()
  de_csvs <- de_csvs[basename(de_csvs) %in% paste0(wanted, ".csv")]
}
cat(sprintf("[1] rendering %d cell types from %s\n",
            length(de_csvs), args$de_dir))

# Pretty display name
pretty_l <- function(l) {
  s <- sub("^[a-z]+__", "", l)
  gsub("_", "-", s)
}

# ---- build one volcano -------------------------------------------------
build_volcano <- function(csv_path) {
  cell <- sub("\\.csv$", "", basename(csv_path))
  hdr <- readr::read_lines(csv_path, n_max = 1)
  if (grepl("^skipped", hdr)) return(NULL)

  de <- read_csv(csv_path, show_col_types = FALSE) |>
    filter(!is.na(adj.P.Val), !is.na(logFC)) |>
    mutate(
      neglog_fdr = -log10(pmax(adj.P.Val, args$fdr_floor)),
      signed_sig = sign(logFC) * neglog_fdr,
      is_sig     = adj.P.Val < args$fdr & abs(logFC) > args$lfc,
      is_marker  = symbol %in% all_markers
    )

  sig_cap <- max(quantile(abs(de$signed_sig), 0.99, na.rm = TRUE), 1.0)

  # Pre-sort so high-|signed_sig| points draw last (on top)
  de <- de |> arrange(abs(signed_sig))

  # Label significant markers only
  de_label <- de |> filter(is_marker, is_sig)

  p <- ggplot(de, aes(x = logFC, y = neglog_fdr, color = signed_sig)) +
    geom_point(size = 0.5, alpha = 0.75, shape = 16) +
    scale_color_gradient2(
      low = div_scale$low, mid = "grey90", high = div_scale$high,
      midpoint = 0, limits = c(-sig_cap, sig_cap),
      oob = scales::squish, guide = "none") +
    geom_hline(yintercept = -log10(args$fdr), linetype = "dashed",
               color = "grey50", linewidth = 0.25) +
    geom_vline(xintercept = c(-args$lfc, args$lfc), linetype = "dashed",
               color = "grey50", linewidth = 0.25) +
    ggrepel::geom_text_repel(data = de_label,
                              aes(label = symbol),
                              inherit.aes = TRUE, size = 2.0, color = "grey15",
                              segment.size = 0.15, segment.color = "grey50",
                              max.overlaps = Inf, min.segment.length = 0,
                              box.padding = 0.3, point.padding = 0.2,
                              force = 2, seed = 42, show.legend = FALSE) +
    labs(x = expression(log[2]~FC~(BR1 / AR)),
         y = expression(-log[10]~FDR)) +
    panel_theme +
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.text  = element_text(size = 6),
          axis.title = element_text(size = 7))

  out_pdf <- file.path(args$out_dir, paste0(cell, ".pdf"))
  ggsave(out_pdf, p, width = 3.5, height = 3.2, units = "in",
         device = cairo_pdf)
  cat(sprintf("  %s -> %s (%d sig markers labeled)\n",
              cell, out_pdf, nrow(de_label)))
  p
}

for (csv in de_csvs) build_volcano(csv)
cat("\nDone.\n")
