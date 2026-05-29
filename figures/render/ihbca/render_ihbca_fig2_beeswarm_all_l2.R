#!/usr/bin/env Rscript
# All-L2 per-nhood Milo DA beeswarm. Used for the supplemental omnibus panels
# (S1-omnibus-tested and S1-omnibus-full) covering every L2 cell type in a
# single figure.
#
# Reads per-neighborhood DA results (one row per nhood) and plots one
# quasi-random swarm row per L2 cell type. Each point is one neighborhood at
# (logFC, signed -log10(SpatialFDR)).
#
# Aesthetic conventions (Fig 2 series):
#   - Diverging color scale by sign(logFC) × −log10(SpatialFDR), with the
#     expression_diverging.low / .high anchors from publication aesthetics.
#   - Vertical dashed line at logFC = 0.
#   - Horizontal dotted separators between L1 cell-type blocks (only when
#     L1-diagonal default ordering is in effect; suppressed when an external
#     L2 order file is supplied).
#   - cairo_pdf output, no on-plot title.
#
# Y-axis ordering: default = per-L1 median-logFC diagonal block ordering;
# override = pass --l2-order-file to lock the order to an external canonical
# list (used for the supp omnibus panels so the y-axis matches the Fig 1
# marker heatmap L2 order).

suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(scales)
})
have_beeswarm <- requireNamespace("ggbeeswarm", quietly = TRUE)
have_rastr <- requireNamespace("ggrastr", quietly = TRUE)

`%||%` <- function(a, b) if (is.null(a)) b else a

parser <- ArgumentParser()
parser$add_argument("--da-csv", required = TRUE,
                    help = "Per-nhood DA results CSV (compartment, label, logFC, SpatialFDR)")
parser$add_argument("--publication-root", required = TRUE)
parser$add_argument("--out-dir", required = TRUE)
parser$add_argument("--panel-id", default = "fig2_ihbca_beeswarm_tested")
parser$add_argument("--fdr-floor", type = "double", default = 1e-15)
parser$add_argument("--n-sig-min", type = "integer", default = 1,
                    help = "Keep L2s with at least this many sig nhoods (FDR<0.10)")
parser$add_argument("--l2-order-file", default = NULL,
                    help = paste("Optional path to a text file (one L2 per line) defining canonical y-axis order.",
                                 "Each entry is a compartment::label string (e.g., 'str::Fibro-major').",
                                 "When supplied, this order overrides the default L1-diagonal-block ordering",
                                 "so the supp omnibus panel can match the Fig 1 marker heatmap L2 order."))
args <- parser$parse_args()

source(file.path(args$publication_root, "publication", "config",
                 "load_aesthetics.R"))
aes_cfg <- load_aesthetics(file.path(args$publication_root, "publication",
                                     "config"))
panel_theme <- get_theme(aes_cfg)
div_scale   <- aes_cfg$scales$expression_diverging

dir.create(args$out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Load DA ---
cat(sprintf("[1] load DA: %s\n", args$da_csv))
da <- read_csv(args$da_csv, show_col_types = FALSE) |>
  filter(!is.na(compartment), !is.na(label), !is.na(SpatialFDR), !is.na(logFC))
cat(sprintf("  nhoods: %d\n", nrow(da)))

# --- L2 label (compartment::label) and artifact filter ---
labels_path <- file.path(args$publication_root, "publication", "analysis",
                          "annotation", "labels_full.csv")
artifact_l2s <- character(0)
if (file.exists(labels_path)) {
  art <- read_csv(labels_path, show_col_types = FALSE,
                  col_select = c("compartment", "label", "is_artifact")) |>
    mutate(is_artifact = as.logical(is_artifact)) |>
    group_by(compartment, label) |>
    summarise(pct = mean(is_artifact, na.rm = TRUE), .groups = "drop") |>
    filter(pct >= 0.5)
  artifact_l2s <- paste(art$compartment, art$label, sep = "::")
  cat(sprintf("  artifact L2s (>=50%% flagged): %d\n", length(artifact_l2s)))
}

da <- da |>
  mutate(L2 = paste(compartment, label, sep = "::")) |>
  filter(!L2 %in% artifact_l2s)

# --- L2 selection: any L2 with >= N_SIG_MIN sig nhoods (FDR < 0.10) ---
l2_summary <- da |>
  group_by(L2, compartment) |>
  summarise(n_sig = sum(SpatialFDR < 0.10, na.rm = TRUE),
            med_lfc = median(logFC, na.rm = TRUE),
            .groups = "drop")
l2_keep <- l2_summary |> filter(n_sig >= args$n_sig_min)
cat(sprintf("[2] L2s with >= %d sig nhoods (FDR<0.10): %d / %d\n",
            args$n_sig_min, nrow(l2_keep), nrow(l2_summary)))

da_sub <- da |> filter(L2 %in% l2_keep$L2)

# --- L1 mapping for diagonal block ordering ---
if (file.exists(labels_path)) {
  l2_to_l1 <- read_csv(labels_path, show_col_types = FALSE,
                       col_select = c("compartment", "label", "lineage")) |>
    group_by(compartment, label, lineage) |>
    summarise(n_cells = dplyr::n(), .groups = "drop_last") |>
    arrange(desc(n_cells)) |>
    slice_head(n = 1) |> ungroup() |>
    mutate(L2 = paste(compartment, label, sep = "::")) |>
    select(L2, L1 = lineage)
  da_sub <- da_sub |> left_join(l2_to_l1, by = "L2") |>
    mutate(L1 = ifelse(is.na(L1), "unknown", L1))
} else {
  da_sub$L1 <- "unknown"
}

# --- Continuous color: signed -log10(SpatialFDR) ---
da_sub <- da_sub |>
  mutate(neglog_fdr = -log10(pmax(SpatialFDR, args$fdr_floor)),
         signed_sig = sign(logFC) * neglog_fdr)
sig_cap <- max(quantile(abs(da_sub$signed_sig), 0.99, na.rm = TRUE), 1.0)

# --- Subsample ns nhoods for visual density (preserves all sig, caps ns) ---
NS_SIG_THRESHOLD <- 0.10
NS_SAMPLE_CAP <- 1500
set.seed(42)
n_before <- nrow(da_sub)
sub_sig <- da_sub |> filter(SpatialFDR < NS_SIG_THRESHOLD)
sub_ns  <- da_sub |> filter(SpatialFDR >= NS_SIG_THRESHOLD) |>
  group_by(L2) |> slice_sample(n = NS_SAMPLE_CAP) |> ungroup()
da_sub <- bind_rows(sub_sig, sub_ns) |>
  arrange(L2, abs(signed_sig))
cat(sprintf("[3] subsample: %d sig + %d ns (kept %d / %d total)\n",
            nrow(sub_sig), nrow(sub_ns), nrow(da_sub), n_before))

# --- Y-axis ordering ---
# If --l2-order-file is supplied, the y-axis follows the external canonical
# L2 order (intersected with the L2s actually present in the data). Otherwise
# fall back to L1-diagonal-block ordering by per-L1 median logFC.
if (!is.null(args$l2_order_file)) {
  cat(sprintf("[4] L2 order from external file: %s\n", args$l2_order_file))
  external_order <- readLines(args$l2_order_file)
  external_order <- trimws(external_order)
  external_order <- external_order[nzchar(external_order) & !startsWith(external_order, "#")]
  # Intersect with L2s present in the panel, preserving file order
  present <- intersect(external_order, l2_keep$L2)
  missing_from_file <- setdiff(l2_keep$L2, external_order)
  if (length(missing_from_file) > 0) {
    cat(sprintf("  WARNING: %d L2s in data not listed in order file; appended in alpha order: %s\n",
                length(missing_from_file),
                paste(missing_from_file, collapse = ", ")))
    present <- c(present, sort(missing_from_file))
  }
  l2_order_df <- l2_keep |>
    left_join(da_sub |> distinct(L2, L1), by = "L2") |>
    mutate(L2 = factor(L2, levels = present)) |>
    arrange(L2) |>
    mutate(L2 = as.character(L2),
           L2_short = sub("^[^:]*::", "", L2))
} else {
  # Default: L1 diagonal ordering by per-L1 median LFC
  l1_order_tbl <- da_sub |>
    group_by(L1) |>
    summarise(l1_med_lfc = median(logFC, na.rm = TRUE), .groups = "drop") |>
    arrange(l1_med_lfc)
  l2_order_df <- l2_keep |>
    left_join(da_sub |> distinct(L2, L1), by = "L2") |>
    mutate(L1 = factor(L1, levels = l1_order_tbl$L1)) |>
    arrange(L1, L2) |>
    mutate(L2_short = sub("^[^:]*::", "", L2))
}
da_sub$L2 <- factor(da_sub$L2, levels = rev(l2_order_df$L2))
y_labels_named <- setNames(l2_order_df$L2_short, l2_order_df$L2)
l2_order_df$y_pos <- rev(seq_len(nrow(l2_order_df)))

# --- L1 separator x-positions (now vertical separators since L2 is on x-axis) ---
sep_x <- if (is.null(args$l2_order_file)) {
  l1_breaks <- l2_order_df |>
    group_by(L1) |>
    summarise(y_min = min(y_pos), y_max = max(y_pos), .groups = "drop") |>
    arrange(y_min)
  l1_breaks |> arrange(desc(y_max)) |>
    mutate(sep = y_min - 0.5) |> filter(sep > 0.5) |> pull(sep)
} else {
  numeric(0)
}

# --- Swarm + violin layers (rotated horizontal: L2 on x-axis, logFC on y-axis;
#     tight-dot + violin-outline aesthetic ported from v3 2b panel) ---
y_lim <- max(max(abs(da_sub$logFC), na.rm = TRUE) * 1.05, 1)
swarm_geom <- if (have_beeswarm) {
  ggbeeswarm::geom_quasirandom(
    aes(x = L2, y = logFC, color = signed_sig),
    groupOnX = TRUE, size = 0.10, alpha = 0.75,
    bandwidth = 0.35, width = 0.27)
} else {
  geom_jitter(aes(x = L2, y = logFC, color = signed_sig),
              width = 0.4, height = 0.0, size = 0.10, alpha = 0.75)
}
swarm_layer <- if (have_rastr) {
  ggrastr::rasterise(swarm_geom, dpi = 600, dev = "cairo")
} else swarm_geom

violin_layer <- geom_violin(
  aes(x = L2, y = logFC, group = L2),
  orientation = "x",
  fill = NA, color = "black", linewidth = 0.30,
  scale = "width", width = 0.55, linetype = "solid")

p <- ggplot(da_sub) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3, linetype = "dashed") +
  geom_vline(xintercept = sep_x, color = "grey50",
             linewidth = 0.3, linetype = "dotted") +
  swarm_layer +
  violin_layer +
  scale_color_gradient2(low = div_scale$low, mid = "grey90", high = div_scale$high,
                        midpoint = 0, limits = c(-sig_cap, sig_cap),
                        oob = scales::squish,
                        name = "sign(logFC) × −log10(SpatialFDR)") +
  scale_x_discrete(labels = setNames(l2_order_df$L2_short, l2_order_df$L2),
                   limits = as.character(l2_order_df$L2),
                   expand = expansion(add = 0.5)) +
  coord_cartesian(ylim = c(-y_lim, y_lim), clip = "off") +
  labs(x = NULL, y = "logFC (per nhood)") +
  panel_theme +
  theme(strip.background = element_blank(),
        legend.position  = "right",
        legend.direction = "vertical",
        legend.title     = element_text(size = 6, angle = 90, hjust = 0.5),
        legend.title.position = "left",
        legend.key.height = unit(1.2, "cm"),
        legend.key.width  = unit(0.30, "cm"),
        panel.grid       = element_blank(),
        panel.background = element_blank(),
        axis.text.x      = element_text(size = 6, face = "bold", angle = 45,
                                        hjust = 1, vjust = 1),
        axis.text.y      = element_text(size = 7),
        axis.title.y     = element_text(size = 8),
        axis.ticks.x     = element_blank(),
        plot.margin      = margin(4, 8, 4, 4))

# Sized as a wide-short strip: width scales with L2 count, fixed-ish height
n_l2 <- nrow(l2_order_df)
panel_w <- max(11.0, 0.32 * n_l2 + 4.0)  # ~17.4" for 42 L2s
panel_h <- 2.8
pdf_path <- file.path(args$out_dir, paste0(args$panel_id, ".pdf"))
ggsave(pdf_path, plot = p, width = panel_w, height = panel_h, units = "in",
       device = cairo_pdf)
cat(sprintf("Wrote %s (%gx%g in; %d L2s)\n", pdf_path, panel_w, panel_h, n_l2))

# PNG companion for composite assembly
png_path <- file.path(args$out_dir, paste0(args$panel_id, ".png"))
ggsave(png_path, plot = p, width = panel_w, height = panel_h, units = "in",
       device = "png", dpi = 300)
cat(sprintf("Wrote %s\n", png_path))
