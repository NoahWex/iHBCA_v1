#!/usr/bin/env Rscript
# Per-L1 focused Milo DA beeswarm renderer. Produces one PDF per L1 cell-type
# split (Fibroblast, T-cell, Stromal non-fibro) at the tested cohort
# (BR1_vs_AR_tested). The caller (render_ihbca_fig2_focused_beeswarm_per_l1_split.sh)
# invokes this three times with different --populations subsets and panel IDs.
#
# Outputs: <out-dir>/<panel-id>.pdf
#
# Plot conventions:
#   - 3-column compartment facet (epi / imm / str), drop empty compartments.
#   - Quasi-random swarm of per-nhood (logFC, signed -log10(SpatialFDR)).
#   - Continuous diverging color scale anchored on expression_diverging.
#   - Per-L2 max-significance marker overlay (* <0.10, ** <0.05, *** <0.01).
#   - 1.5-column panel size (120 mm width default).

suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(scales)
})
have_beeswarm <- requireNamespace("ggbeeswarm", quietly = TRUE)
have_rastr    <- requireNamespace("ggrastr",    quietly = TRUE)

`%||%` <- function(a, b) if (is.null(a)) b else a

parser <- ArgumentParser()
parser$add_argument("--da-csv", required = TRUE,
                    help = "Per-nhood DA results (compartment, label, logFC, SpatialFDR)")
parser$add_argument("--publication-root", required = TRUE)
parser$add_argument("--out-dir", required = TRUE)
parser$add_argument("--populations", required = TRUE,
                    help = "Comma-sep keys: 'compartment__label' (e.g., 'str__Fibro-major')")
parser$add_argument("--panel-id", default = "fig2_ihbca_focused_beeswarm_tested")
parser$add_argument("--fdr-floor", type = "double", default = 1e-15)
parser$add_argument("--width-in", type = "double", default = 5.0,
                    help = "Total PDF width in inches")
parser$add_argument("--per-row-height-in", type = "double", default = 0.28,
                    help = "Height per L2 swarm row in inches (default 0.28 ≈ 1:12 row aspect at width 5.0)")
parser$add_argument("--height-overhead-in", type = "double", default = 1.0,
                    help = "Fixed overhead for axis title + facet strips + plot margins")
args <- parser$parse_args()

source(file.path(args$publication_root, "publication", "config",
                 "load_aesthetics.R"))
aes_cfg <- load_aesthetics(file.path(args$publication_root, "publication",
                                     "config"))
panel_theme <- get_theme(aes_cfg)
div_scale   <- aes_cfg$scales$expression_diverging

dir.create(args$out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Parse population subset ---
pop_keys <- strsplit(args$populations, ",", fixed = TRUE)[[1]] |> trimws()
pop_df <- data.frame(
  key = pop_keys,
  compartment = sub("__.*$", "", pop_keys),
  label = sub("^[^_]*__", "", pop_keys),
  stringsAsFactors = FALSE
)
cat(sprintf("[1] %d populations requested: %s\n",
            nrow(pop_df), paste(pop_df$key, collapse = ", ")))

# --- Load DA + filter to populations ---
da <- read_csv(args$da_csv, show_col_types = FALSE) |>
  filter(!is.na(compartment), !is.na(label),
         !is.na(SpatialFDR), !is.na(logFC))

da_sub <- da |>
  inner_join(pop_df, by = c("compartment", "label"))
cat(sprintf("[2] nhoods after pop filter: %d / %d\n",
            nrow(da_sub), nrow(da)))

if (nrow(da_sub) == 0) stop("No nhoods matched populations: ", args$populations)

# --- Color: signed -log10(SpatialFDR) ---
da_sub <- da_sub |>
  mutate(neglog_fdr = -log10(pmax(SpatialFDR, args$fdr_floor)),
         signed_sig = sign(logFC) * neglog_fdr)
sig_cap <- max(quantile(abs(da_sub$signed_sig), 0.99, na.rm = TRUE), 1.0)

# --- Compartment facet ordering: epi / imm / str ---
COMP_LEVELS <- c("epi", "imm", "str")
COMP_LABELS <- c(epi = "Epithelial", imm = "Immune", str = "Stromal")
da_sub$compartment <- factor(da_sub$compartment, levels = COMP_LEVELS,
                              labels = COMP_LABELS)
pop_df$compartment_f <- factor(pop_df$compartment, levels = COMP_LEVELS,
                                labels = COMP_LABELS)

# --- Y-axis order per compartment (lock from --populations order) ---
y_order <- pop_df |> arrange(match(compartment, COMP_LEVELS)) |> pull(label)
da_sub$label <- factor(da_sub$label, levels = rev(y_order))

# --- Per-L2 significance annotation (max significance across nhoods) ---
sig_anno <- da_sub |>
  group_by(compartment, label) |>
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
# Keep all sig + larger ns subsample to give swarms visible density
NS_SIG_THRESHOLD <- 0.10
NS_SAMPLE_CAP <- 1200
set.seed(42)
sub_sig <- da_sub |> filter(SpatialFDR <  NS_SIG_THRESHOLD)
sub_ns  <- da_sub |> filter(SpatialFDR >= NS_SIG_THRESHOLD) |>
  group_by(compartment, label) |> slice_sample(n = NS_SAMPLE_CAP) |> ungroup()
da_sub <- bind_rows(sub_sig, sub_ns) |>
  arrange(compartment, label, abs(signed_sig))
cat(sprintf("[3] subsample: %d sig + %d ns (kept %d)\n",
            nrow(sub_sig), nrow(sub_ns), nrow(da_sub)))

# Drop empty compartments (e.g., epi placeholder) so they don't render
present_comps <- da_sub |> dplyr::distinct(compartment) |> dplyr::pull(compartment)
da_sub$compartment <- droplevels(factor(da_sub$compartment, levels = present_comps))

# --- Swarm layer (dense-swarm aesthetic ported from
#     coordination/staging/fig3_promotion_preview_20260520/figures/render/fig3/
#     render_xenium_beeswarm_grid_with_diagrams.R: tiny tight-packed dots) ---
x_lim <- max(max(abs(da_sub$logFC), na.rm = TRUE) * 1.05, 1)
swarm_geom <- if (have_beeswarm) {
  ggbeeswarm::geom_quasirandom(
    aes(x = logFC, y = label, color = signed_sig),
    groupOnX = FALSE, size = 0.10, alpha = 0.75,
    bandwidth = 0.35, width = 0.27)
} else {
  geom_jitter(aes(x = logFC, y = label, color = signed_sig),
              width = 0.0, height = 0.4, size = 0.10, alpha = 0.75)
}
# Rasterise swarm: facets are non-empty after droplevels (~5400 points)
swarm_layer <- if (have_rastr) {
  ggrastr::rasterise(swarm_geom, dpi = 600, dev = "cairo")
} else {
  swarm_geom
}

# --- Violin outline layer (dense-swarm aesthetic) ---
# Sharp black outline around each L2's density envelope. fill=NA so the swarm
# stays visible through the violin shape. Ported from fig3 grid renderer.
violin_layer <- geom_violin(
  aes(x = logFC, y = label, group = label),
  orientation = "y",
  fill = NA, color = "black", linewidth = 0.30,
  scale = "width", width = 0.55, linetype = "solid")

p <- ggplot(da_sub) +
  geom_vline(xintercept = 0, color = "grey60",
             linewidth = 0.3, linetype = "dashed") +
  swarm_layer +
  violin_layer +
  scale_color_gradient2(low = div_scale$low, mid = "grey90", high = div_scale$high,
                        midpoint = 0, limits = c(-sig_cap, sig_cap),
                        oob = scales::squish,
                        name = "sign(logFC) × −log10(SpatialFDR)") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  coord_cartesian(xlim = c(-x_lim, x_lim), clip = "off") +
  facet_grid(compartment ~ ., drop = TRUE, scales = "free_y", space = "free_y", switch = "y") +
  labs(x = "logFC (per nhood)", y = NULL) +
  panel_theme +
  theme(strip.background = element_blank(),
        strip.text.y.left = element_text(face = "bold", size = 8, angle = 90),
        legend.position  = "right",
        legend.direction = "vertical",
        legend.title     = element_text(size = 6, angle = 90, hjust = 0.5),
        legend.title.position = "left",
        legend.key.height = unit(0.9, "cm"),
        legend.key.width  = unit(0.30, "cm"),
        panel.grid       = element_blank(),
        panel.background = element_blank(),
        axis.text.y      = element_text(size = 6, face = "bold"),
        axis.title.x     = element_text(size = 8),
        axis.ticks.y     = element_blank(),
        panel.spacing.y  = unit(0.5, "lines"),
        plot.margin      = margin(4, 8, 4, 4),
        strip.placement  = "outside")

# Dynamic height: n_pops * per_row_height + overhead, locks 1:4 row aspect at width=5
n_pops <- nrow(pop_df)
pdf_height <- n_pops * args$per_row_height_in + args$height_overhead_in
pdf_path <- file.path(args$out_dir, paste0(args$panel_id, ".pdf"))
ggsave(pdf_path, plot = p, width = args$width_in, height = pdf_height, units = "in",
       device = cairo_pdf)
cat(sprintf("Wrote %s (%.2f x %.2f in; %d pops; row height %.2f in)\n",
            pdf_path, args$width_in, pdf_height, n_pops, args$per_row_height_in))

# --- Audit summary ---
cat("\n=== sig audit ===\n")
print(sig_anno |> select(compartment, label, min_fdr, n_sig_01, n_sig_05, n_sig_10, mark))
