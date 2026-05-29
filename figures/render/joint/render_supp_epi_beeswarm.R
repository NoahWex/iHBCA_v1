#!/usr/bin/env Rscript
# Single-question Milo beeswarm for the C5 epi_content contrast.
#
# Question per panel: for the given scope (e.g. flex_Immune), which L1.5 cell
# types co-enrich (logFC > 0) or co-deplete (logFC < 0) with total per-library
# epithelial content?
#
# Plot structure:
#   Y-axis: L1.5 label (nhood_label_l1p5) ordered by median signed logFC
#           within the scope
#   X-axis: logFC from Milo edgeR GLM (epi_rich vs epi_sparse, blocked
#           by patient_id + p_level)
#   Color:  signed logFC (RdBu diverging) on significant nhoods (SpatialFDR < 0.05);
#           non-significant in grey
#   Shape:  none (single platform per panel)
#   Refs:   x = 0 (vertical line, neutral)
#
# Deliberately omitted (single-question discipline):
#   - panel title / methods text (lives in figure legend)
#   - panel letter (added in Illustrator)
#   - L1.0 collapse separators (one resolution per panel)
#   - per-platform shape overlay (FLEX-only panel needs no shape variation)

suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(ggplot2)
  library(ggbeeswarm)
})

# load_aesthetics() uses %||% internally; define here so it's in scope
`%||%` <- function(a, b) if (is.null(a)) b else a

parser <- ArgumentParser()
parser$add_argument("--project-root", required = TRUE)
parser$add_argument("--da-csv",       required = TRUE,
                    help = "Path to da_results.csv (one nhood per row)")
parser$add_argument("--label-col",    default = "nhood_label_l1p5",
                    help = "Column name in da-csv carrying the L1.5 label per nhood")
parser$add_argument("--out-pdf",      required = TRUE)
parser$add_argument("--sig-fdr",      type = "double", default = 0.05)
parser$add_argument("--min-nhoods",   type = "integer", default = 10,
                    help = "Drop L1.5 labels with fewer than this many nhoods")
parser$add_argument("--width-in",     type = "double", default = 4.7)
parser$add_argument("--height-in",    type = "double", default = 5.0)
args <- parser$parse_args()

# Aesthetics framework
root <- normalizePath(args$project_root)
source(file.path(root, "publication", "config", "load_aesthetics.R"))
aes_cfg <- load_aesthetics(file.path(root, "publication", "config"))
div <- aes_cfg$scales$expression_diverging

# Load DA
da <- read.csv(args$da_csv, stringsAsFactors = FALSE)
stopifnot(args$label_col %in% colnames(da))
stopifnot(all(c("logFC", "SpatialFDR") %in% colnames(da)))

# Filter: drop unlabeled nhoods + low-population labels
da <- da[!is.na(da[[args$label_col]]) & da[[args$label_col]] != "NA", ]
counts <- table(da[[args$label_col]])
keep_labels <- names(counts)[counts >= args$min_nhoods]
da <- da[da[[args$label_col]] %in% keep_labels, ]
cat("[beeswarm] nhoods=", nrow(da), " distinct labels=", length(keep_labels), "\n", sep = "")

# Order labels by median signed logFC of sig nhoods (fallback to all-nhood median)
da$is_sig  <- da$SpatialFDR < args$sig_fdr
da$plot_lfc <- ifelse(da$is_sig, da$logFC, NA_real_)
med_sig <- tapply(da$logFC[da$is_sig], da[[args$label_col]][da$is_sig], median, na.rm = TRUE)
med_all <- tapply(da$logFC,             da[[args$label_col]],            median, na.rm = TRUE)
med_for_order <- ifelse(is.na(med_sig[names(med_all)]), med_all, med_sig[names(med_all)])
label_order   <- names(sort(med_for_order))
da$label_f    <- factor(da[[args$label_col]], levels = label_order)

# Color limits: symmetric around 0, capped at the 98th percentile of |logFC| sig
lfc_cap <- quantile(abs(da$logFC[da$is_sig]), 0.98, na.rm = TRUE)
if (!is.finite(lfc_cap) || lfc_cap == 0) lfc_cap <- max(abs(da$logFC), na.rm = TRUE)

p <- ggplot(da, aes(x = logFC, y = label_f)) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey60") +
  # non-sig nhoods: grey background
  geom_quasirandom(data = subset(da, !is_sig),
                   groupOnX = FALSE, method = "tukeyDense",
                   size = 0.45, alpha = 0.45, color = "grey75", stroke = 0) +
  # sig nhoods: diverging color by signed logFC
  geom_quasirandom(data = subset(da, is_sig),
                   aes(color = plot_lfc),
                   groupOnX = FALSE, method = "tukeyDense",
                   size = 0.65, alpha = 0.85, stroke = 0) +
  scale_color_gradient2(low = div$low, mid = div$mid, high = div$high,
                        midpoint = 0, limits = c(-lfc_cap, lfc_cap),
                        oob = scales::squish,
                        guide = guide_colorbar(title = "logFC", barwidth = 0.4,
                                               barheight = 4)) +
  labs(x = "logFC (epi_rich vs epi_sparse)", y = NULL) +
  get_theme(aes_cfg) +
  theme(panel.grid.major.y = element_line(color = "grey92", linewidth = 0.2),
        panel.grid.minor   = element_blank(),
        plot.title         = element_blank(),
        plot.subtitle      = element_blank(),
        legend.position    = "right",
        axis.text.y        = element_text(size = 6))

# Save with explicit dims (single-column-narrow Nature spec, vector PDF, no rasterization)
dir.create(dirname(args$out_pdf), recursive = TRUE, showWarnings = FALSE)
ggsave(args$out_pdf, plot = p, width = args$width_in, height = args$height_in,
       units = "in", device = cairo_pdf)
cat("[beeswarm] wrote ", args$out_pdf, "\n", sep = "")
