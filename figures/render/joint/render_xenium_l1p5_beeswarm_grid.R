#!/usr/bin/env Rscript
# Fig 3 / supp — Xenium L1.5 × N-contrasts Milo DA beeswarm grid.
#
# Companion to render_xenium_motif_beeswarm_grid.R but grouping by L1.5
# cell type (nhood_label_l1p5) instead of motif.  Reads the same DA CSVs;
# L1.5 column is already present, so no motif annotation join required.
#
# Substrate:
#   --contrasts-csv  CSV with columns: contrast_label, da_csv_path
#   --l1p5-frac-min  Filter nhoods by minimum nhood_label_l1p5_frac (default
#                    0.7 to drop mixed-identity nhoods).

suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(scales)
})
have_beeswarm <- requireNamespace("ggbeeswarm", quietly = TRUE)
have_rastr    <- requireNamespace("ggrastr",    quietly = TRUE)

parser <- ArgumentParser()
parser$add_argument("--contrasts-csv", required = TRUE)
parser$add_argument("--publication-root", required = TRUE)
parser$add_argument("--root", default = "")
parser$add_argument("--out-dir", required = TRUE)
parser$add_argument("--panel-id",
                    default = "s_xenium_l1p5_beeswarm_grid_all_contrasts")
parser$add_argument("--l1p5-frac-min", type = "double", default = 0.7,
                    help = "Drop nhoods with mixed-identity (nhood_label_l1p5_frac < this)")
parser$add_argument("--fdr-floor", type = "double", default = 1e-15)
parser$add_argument("--width-per-contrast-in", type = "double", default = 1.1)
parser$add_argument("--width-legend-in", type = "double", default = 1.2)
parser$add_argument("--per-row-height-in", type = "double", default = 0.20,
                    help = "L1.5 has ~20 cell types so rows are shorter")
parser$add_argument("--height-overhead-in", type = "double", default = 1.4)
args <- parser$parse_args()

source(file.path(args$publication_root, "publication", "config",
                 "load_aesthetics.R"))
aes_cfg <- load_aesthetics(file.path(args$publication_root, "publication",
                                     "config"))
panel_theme <- get_theme(aes_cfg)
div_scale   <- aes_cfg$scales$expression_diverging

dir.create(args$out_dir, showWarnings = FALSE, recursive = TRUE)

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
      filter(!is.na(SpatialFDR), !is.na(logFC),
             !is.na(nhood_label_l1p5),
             nhood_label_l1p5_frac >= args$l1p5_frac_min) |>
      mutate(contrast = .$contrast_label,
             l1p5 = nhood_label_l1p5)
    df
  }) |>
  ungroup() |>
  as_tibble()
cat(sprintf("[2] loaded %d nhoods (l1p5_frac >= %.2f)\n",
            nrow(da_all), args$l1p5_frac_min))

# --- Order L1.5 cell types: epi → stromal → endothelial → immune ---
L1P5_ORDER <- c(
  # Epithelial
  "BMYO-myo", "LASP-basal", "LASP", "LHS",
  # Stromal
  "Fb", "Fb_Activated", "Fb_SFRP4", "Adipo",
  # Endothelial / vasculature
  "EC", "Vas-cap", "PV", "LEC",
  # T / NK
  "CD4 T cell", "CD4T", "CD8 T cell", "CD8T", "Treg", "T-NK",
  # Dendritic
  "cDC", "cDC1", "cDC2", "pDC",
  # B / Plasma
  "B cell", "B", "Plasma cell", "Plas",
  # Macrophage / mast / neutrophil
  "Macrophage", "Mac", "Mac_art", "Mast", "Neutrophil", "Neu"
)
present <- unique(da_all$l1p5)
ordered_present <- c(intersect(L1P5_ORDER, present),
                     setdiff(present, L1P5_ORDER))

da_all$l1p5 <- factor(da_all$l1p5, levels = rev(ordered_present))
da_all$contrast <- factor(da_all$contrast, levels = contrasts$contrast_label)

# --- Color: signed −log10(SpatialFDR), shared scale ---
da_all <- da_all |>
  mutate(neglog_fdr = -log10(pmax(SpatialFDR, args$fdr_floor)),
         signed_sig = sign(logFC) * neglog_fdr)
sig_cap <- max(quantile(abs(da_all$signed_sig), 0.99, na.rm = TRUE), 1.0)

# --- Subsample ns ---
NS_SIG_THRESHOLD <- 0.10
NS_SAMPLE_CAP <- 600  # per (contrast, l1p5)
set.seed(42)
sub_sig <- da_all |> filter(SpatialFDR <  NS_SIG_THRESHOLD)
sub_ns  <- da_all |> filter(SpatialFDR >= NS_SIG_THRESHOLD) |>
  group_by(contrast, l1p5) |> slice_sample(n = NS_SAMPLE_CAP) |> ungroup()
da_plot <- bind_rows(sub_sig, sub_ns) |>
  arrange(contrast, l1p5, abs(signed_sig))
cat(sprintf("[3] subsample kept %d / %d\n", nrow(da_plot), nrow(da_all)))

x_lim <- max(max(abs(da_plot$logFC), na.rm = TRUE) * 1.05, 1)

swarm_geom <- if (have_beeswarm) {
  ggbeeswarm::geom_quasirandom(
    aes(x = logFC, y = l1p5, color = signed_sig),
    groupOnX = FALSE, size = 0.35, alpha = 0.75,
    bandwidth = 0.35, width = 0.45)
} else {
  geom_jitter(aes(x = logFC, y = l1p5, color = signed_sig),
              width = 0.0, height = 0.4, size = 0.35, alpha = 0.75)
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
        axis.text.y      = element_text(size = 5.5),
        axis.text.x      = element_text(size = 5),
        axis.title.x     = element_text(size = 7),
        axis.ticks.y     = element_blank(),
        panel.spacing.x  = unit(0.4, "lines"),
        plot.margin      = margin(4, 8, 4, 4))

n_l1p5 <- length(ordered_present)
n_contrasts <- nrow(contrasts)
pdf_width  <- n_contrasts * args$width_per_contrast_in + args$width_legend_in
pdf_height <- n_l1p5 * args$per_row_height_in + args$height_overhead_in
pdf_path <- file.path(args$out_dir, paste0(args$panel_id, ".pdf"))
ggsave(pdf_path, plot = p, width = pdf_width, height = pdf_height,
       units = "in")
cat(sprintf("Wrote %s (%.2f x %.2f in; %d contrasts × %d L1.5 types)\n",
            pdf_path, pdf_width, pdf_height, n_contrasts, n_l1p5))
