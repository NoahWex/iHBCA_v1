#!/usr/bin/env Rscript
# Fig 3f — Xenium motif × c3 UOQ Milo DA beeswarm (publication_v2 aesthetic).
#
# Ported from publication/figures/render/ihbca/render_ihbca_fig2_focused_beeswarm_per_l1_split.R
# to a single-axis motif-categorical version. Y-axis: M0..M8 motifs. X-axis:
# per-nhood logFC. Color: signed −log10(SpatialFDR) on the diverging
# expression palette. Per-motif min-FDR marker overlay (* <0.10, ** <0.05,
# *** <0.01).
#
# Substrate:
#   Spatial_HBCA_Xenium/DifferentialAbundance/joint_anatomic_variation_20260424/
#     outputs/da/C3_p3_uoq_vs_rest/joint_full/M0/da_results.csv     per-nhood
#                                                                    logFC + SpatialFDR
#   Spatial_HBCA_Xenium/NicheFramework/fig4_lane_c_nmf_20260515/
#     reports/motif_units/motif_x_milo_DA/per_nhood_motif_annot.parquet
#                                                                    nhood_idx → dominant_motif
#
# Output: <out-dir>/3f_xenium_motif_beeswarm_c3_UOQ.pdf

suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(scales)
})
have_arrow <- requireNamespace("arrow", quietly = TRUE)
have_beeswarm <- requireNamespace("ggbeeswarm", quietly = TRUE)
have_rastr    <- requireNamespace("ggrastr",    quietly = TRUE)

`%||%` <- function(a, b) if (is.null(a)) b else a

parser <- ArgumentParser()
parser$add_argument("--da-csv", required = TRUE,
                    help = "Per-nhood DA results (Nhood, logFC, SpatialFDR)")
parser$add_argument("--nhood-motif", required = TRUE,
                    help = "per_nhood_motif_annot.parquet (or .csv) — nhood_idx, dominant_motif")
parser$add_argument("--publication-root", required = TRUE)
parser$add_argument("--out-dir", required = TRUE)
parser$add_argument("--panel-id", default = "3f_xenium_motif_beeswarm_c3_UOQ")
parser$add_argument("--motifs", default = "M0,M1,M2,M3,M4,M5,M6,M7,M8",
                    help = "Comma-separated motif IDs to include / order (top→bottom)")
parser$add_argument("--fdr-floor", type = "double", default = 1e-15)
parser$add_argument("--width-in", type = "double", default = 3.5,
                    help = "PDF width in inches")
parser$add_argument("--per-row-height-in", type = "double", default = 0.28)
parser$add_argument("--height-overhead-in", type = "double", default = 1.0)
args <- parser$parse_args()

source(file.path(args$publication_root, "publication", "config",
                 "load_aesthetics.R"))
aes_cfg <- load_aesthetics(file.path(args$publication_root, "publication",
                                     "config"))
panel_theme <- get_theme(aes_cfg)
div_scale   <- aes_cfg$scales$expression_diverging

dir.create(args$out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Parse motif order ---
motif_order <- strsplit(args$motifs, ",", fixed = TRUE)[[1]] |> trimws()
cat(sprintf("[1] %d motifs requested: %s\n",
            length(motif_order), paste(motif_order, collapse = ", ")))

# --- Load DA + motif annotation, join on Nhood ---
da <- read_csv(args$da_csv, show_col_types = FALSE) |>
  filter(!is.na(SpatialFDR), !is.na(logFC))

nhood_motif <- if (grepl("\\.parquet$", args$nhood_motif)) {
  if (!have_arrow) stop("arrow package required to read parquet; install or pass a CSV path")
  arrow::read_parquet(args$nhood_motif)
} else {
  read_csv(args$nhood_motif, show_col_types = FALSE)
} |>
  select(nhood_idx, dominant_motif)

da_motif <- da |>
  inner_join(nhood_motif, by = c("Nhood" = "nhood_idx")) |>
  filter(dominant_motif %in% motif_order) |>
  rename(motif = dominant_motif)

cat(sprintf("[2] nhoods joined to motif: %d / %d\n",
            nrow(da_motif), nrow(da)))

if (nrow(da_motif) == 0) stop("No nhoods matched motifs: ", args$motifs)

# --- Color: signed −log10(SpatialFDR) ---
da_motif <- da_motif |>
  mutate(neglog_fdr = -log10(pmax(SpatialFDR, args$fdr_floor)),
         signed_sig = sign(logFC) * neglog_fdr)
sig_cap <- max(quantile(abs(da_motif$signed_sig), 0.99, na.rm = TRUE), 1.0)

# --- Y-axis order (top→bottom from --motifs, so reverse for ggplot factor) ---
da_motif$motif <- factor(da_motif$motif, levels = rev(motif_order))

# --- Per-motif significance annotation (max significance across nhoods) ---
sig_anno <- da_motif |>
  group_by(motif) |>
  summarise(min_fdr = min(SpatialFDR, na.rm = TRUE),
            n_sig_01 = sum(SpatialFDR < 0.01),
            n_sig_05 = sum(SpatialFDR < 0.05),
            n_sig_10 = sum(SpatialFDR < 0.10),
            .groups = "drop") |>
  mutate(mark = case_when(
    min_fdr < 0.01 ~ "***",
    min_fdr < 0.05 ~ "**",
    min_fdr < 0.10 ~ "*",
    TRUE           ~ ""
  ))

# --- Subsample ns nhoods for density control ---
NS_SIG_THRESHOLD <- 0.10
NS_SAMPLE_CAP <- 1200
set.seed(42)
sub_sig <- da_motif |> filter(SpatialFDR <  NS_SIG_THRESHOLD)
sub_ns  <- da_motif |> filter(SpatialFDR >= NS_SIG_THRESHOLD) |>
  group_by(motif) |> slice_sample(n = NS_SAMPLE_CAP) |> ungroup()
da_motif <- bind_rows(sub_sig, sub_ns) |>
  arrange(motif, abs(signed_sig))
cat(sprintf("[3] subsample: %d sig + %d ns (kept %d)\n",
            nrow(sub_sig), nrow(sub_ns), nrow(da_motif)))

# --- Swarm layer ---
x_lim <- max(max(abs(da_motif$logFC), na.rm = TRUE) * 1.05, 1)
swarm_geom <- if (have_beeswarm) {
  ggbeeswarm::geom_quasirandom(
    aes(x = logFC, y = motif, color = signed_sig),
    groupOnX = FALSE, size = 0.45, alpha = 0.75,
    bandwidth = 0.35, width = 0.45)
} else {
  geom_jitter(aes(x = logFC, y = motif, color = signed_sig),
              width = 0.0, height = 0.4, size = 0.45, alpha = 0.75)
}
swarm_layer <- if (have_rastr) {
  ggrastr::rasterise(swarm_geom, dpi = 600, dev = "cairo")
} else {
  swarm_geom
}

# --- Build sig-marker annotation positioned at right margin ---
sig_anno$y_pos <- as.numeric(factor(sig_anno$motif, levels = rev(motif_order)))

p <- ggplot(da_motif) +
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
  labs(x = "log2 fold change", y = NULL) +
  panel_theme +
  theme(legend.position  = "right",
        legend.direction = "vertical",
        legend.title     = element_text(size = 6, angle = 90, hjust = 0.5),
        legend.title.position = "left",
        legend.key.height = unit(0.9, "cm"),
        legend.key.width  = unit(0.30, "cm"),
        panel.grid       = element_blank(),
        panel.background = element_blank(),
        axis.text.y      = element_text(size = 7, face = "bold"),
        axis.title.x     = element_text(size = 8),
        axis.ticks.y     = element_blank(),
        plot.margin      = margin(4, 8, 4, 4))

# --- Save ---
n_motifs <- length(motif_order)
pdf_height <- n_motifs * args$per_row_height_in + args$height_overhead_in
pdf_path <- file.path(args$out_dir, paste0(args$panel_id, ".pdf"))
ggsave(pdf_path, plot = p, width = args$width_in, height = pdf_height,
       units = "in")
cat(sprintf("Wrote %s (%.2f x %.2f in; %d motifs; row height %.2f in)\n",
            pdf_path, args$width_in, pdf_height, n_motifs,
            args$per_row_height_in))

# --- Audit summary ---
cat("\n=== sig audit ===\n")
print(sig_anno |> select(motif, min_fdr, n_sig_01, n_sig_05, n_sig_10, mark))
