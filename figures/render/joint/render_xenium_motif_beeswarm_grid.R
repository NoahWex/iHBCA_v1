#!/usr/bin/env Rscript
# Fig 3 / supp — Xenium motif × N-contrasts Milo DA beeswarm grid.
#
# Generic horizontal-beeswarm grid: motifs M0..M8 on Y, logFC on X, faceted
# across an arbitrary number of contrasts (one column per contrast).
# Reads a CSV of (contrast_label, da_csv_path) pairs.  Same publication_v2
# ggbeeswarm aesthetic (signed −log10 FDR diverging color + load_aesthetics
# theme).  Use case: 7-anatomic-contrast master grid for §A abundance supp.
#
# Substrate:
#   --contrasts-csv  CSV with columns: contrast_label, da_csv_path (rows are
#                    one contrast each).  Paths are relative to --root or
#                    absolute.
#   --nhood-motif    per_nhood_motif_annot (parquet OR csv) → nhood_idx,
#                    dominant_motif

suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(scales)
})
have_arrow    <- requireNamespace("arrow",      quietly = TRUE)
have_beeswarm <- requireNamespace("ggbeeswarm", quietly = TRUE)
have_rastr    <- requireNamespace("ggrastr",    quietly = TRUE)

parser <- ArgumentParser()
parser$add_argument("--contrasts-csv", required = TRUE)
parser$add_argument("--nhood-motif", required = TRUE)
parser$add_argument("--publication-root", required = TRUE)
parser$add_argument("--root", default = "",
                    help = "Prefix for relative da_csv_path entries; empty = paths used as-is")
parser$add_argument("--out-dir", required = TRUE)
parser$add_argument("--panel-id",
                    default = "s_xenium_motif_beeswarm_grid_all_contrasts")
parser$add_argument("--motifs", default = "M0,M1,M2,M3,M4,M5,M6,M7,M8")
parser$add_argument("--fdr-floor", type = "double", default = 1e-15)
parser$add_argument("--width-per-contrast-in", type = "double", default = 1.1)
parser$add_argument("--width-legend-in", type = "double", default = 1.2)
parser$add_argument("--per-row-height-in", type = "double", default = 0.30)
parser$add_argument("--height-overhead-in", type = "double", default = 1.4)
args <- parser$parse_args()

source(file.path(args$publication_root, "publication", "config",
                 "load_aesthetics.R"))
aes_cfg <- load_aesthetics(file.path(args$publication_root, "publication",
                                     "config"))
panel_theme <- get_theme(aes_cfg)
div_scale   <- aes_cfg$scales$expression_diverging

dir.create(args$out_dir, showWarnings = FALSE, recursive = TRUE)

motif_order <- strsplit(args$motifs, ",", fixed = TRUE)[[1]] |> trimws()

# --- Load nhood→motif annotation ---
nhood_motif <- if (grepl("\\.parquet$", args$nhood_motif)) {
  if (!have_arrow) stop("arrow package required for parquet")
  arrow::read_parquet(args$nhood_motif)
} else {
  read_csv(args$nhood_motif, show_col_types = FALSE)
} |>
  select(nhood_idx, dominant_motif)

# --- Load contrast spec + each DA CSV, build long format ---
contrasts <- read_csv(args$contrasts_csv, show_col_types = FALSE)
cat(sprintf("[1] %d contrasts to render\n", nrow(contrasts)))

resolve_path <- function(p) {
  if (substr(p, 1, 1) == "/") p else file.path(args$root, p)
}

da_all <- contrasts |>
  rowwise() |>
  do({
    p <- resolve_path(.$da_csv_path)
    if (!file.exists(p)) stop("missing DA CSV: ", p)
    df <- read_csv(p, show_col_types = FALSE) |>
      filter(!is.na(SpatialFDR), !is.na(logFC)) |>
      inner_join(nhood_motif, by = c("Nhood" = "nhood_idx")) |>
      filter(dominant_motif %in% motif_order) |>
      rename(motif = dominant_motif) |>
      mutate(contrast = .$contrast_label)
    df
  }) |>
  ungroup() |>
  as_tibble()
cat(sprintf("[2] loaded %d nhoods across all contrasts\n", nrow(da_all)))

# --- Color: signed −log10(SpatialFDR), shared scale ---
da_all <- da_all |>
  mutate(neglog_fdr = -log10(pmax(SpatialFDR, args$fdr_floor)),
         signed_sig = sign(logFC) * neglog_fdr)
sig_cap <- max(quantile(abs(da_all$signed_sig), 0.99, na.rm = TRUE), 1.0)

# --- Factor ordering ---
da_all$motif <- factor(da_all$motif, levels = rev(motif_order))
da_all$contrast <- factor(da_all$contrast, levels = contrasts$contrast_label)

# --- Subsample ns for density control ---
NS_SIG_THRESHOLD <- 0.10
NS_SAMPLE_CAP <- 1200
set.seed(42)
sub_sig <- da_all |> filter(SpatialFDR <  NS_SIG_THRESHOLD)
sub_ns  <- da_all |> filter(SpatialFDR >= NS_SIG_THRESHOLD) |>
  group_by(contrast, motif) |> slice_sample(n = NS_SAMPLE_CAP) |> ungroup()
da_plot <- bind_rows(sub_sig, sub_ns) |>
  arrange(contrast, motif, abs(signed_sig))
cat(sprintf("[3] subsample kept %d / %d\n", nrow(da_plot), nrow(da_all)))

# --- Shared x-range ---
x_lim <- max(max(abs(da_plot$logFC), na.rm = TRUE) * 1.05, 1)

swarm_geom <- if (have_beeswarm) {
  ggbeeswarm::geom_quasirandom(
    aes(x = logFC, y = motif, color = signed_sig),
    groupOnX = FALSE, size = 0.40, alpha = 0.75,
    bandwidth = 0.35, width = 0.45)
} else {
  geom_jitter(aes(x = logFC, y = motif, color = signed_sig),
              width = 0.0, height = 0.4, size = 0.4, alpha = 0.75)
}
swarm_layer <- if (have_rastr) {
  ggrastr::rasterise(swarm_geom, dpi = 600, dev = "cairo")
} else {
  swarm_geom
}

p <- ggplot(da_plot) +
  geom_vline(xintercept = 0, color = "grey60",
             linewidth = 0.3, linetype = "dashed") +
  swarm_layer +
  scale_color_gradient2(low = div_scale$low, mid = "grey90",
                        high = div_scale$high,
                        midpoint = 0, limits = c(-sig_cap, sig_cap),
                        oob = scales::squish,
                        name = "signed-log10(FDR)") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  coord_cartesian(xlim = c(-x_lim, x_lim), clip = "off") +
  facet_grid(. ~ contrast, switch = "y") +
  labs(x = "logFC (per nhood)", y = NULL) +
  panel_theme +
  theme(strip.background = element_blank(),
        strip.text.x     = element_text(face = "bold", size = 7),
        legend.position  = "right",
        legend.direction = "vertical",
        legend.title     = element_text(size = 6, angle = 90, hjust = 0.5),
        legend.title.position = "left",
        legend.key.height = unit(0.9, "cm"),
        legend.key.width  = unit(0.30, "cm"),
        panel.grid       = element_blank(),
        panel.background = element_blank(),
        axis.text.y      = element_text(size = 6, face = "bold"),
        axis.text.x      = element_text(size = 5),
        axis.title.x     = element_text(size = 7),
        axis.ticks.y     = element_blank(),
        panel.spacing.x  = unit(0.4, "lines"),
        plot.margin      = margin(4, 8, 4, 4))

n_motifs <- length(motif_order)
n_contrasts <- nrow(contrasts)
pdf_width  <- n_contrasts * args$width_per_contrast_in + args$width_legend_in
pdf_height <- n_motifs * args$per_row_height_in + args$height_overhead_in
pdf_path <- file.path(args$out_dir, paste0(args$panel_id, ".pdf"))
ggsave(pdf_path, plot = p, width = pdf_width, height = pdf_height,
       units = "in")
cat(sprintf("Wrote %s (%.2f x %.2f in; %d contrasts × %d motifs)\n",
            pdf_path, pdf_width, pdf_height, n_contrasts, n_motifs))
