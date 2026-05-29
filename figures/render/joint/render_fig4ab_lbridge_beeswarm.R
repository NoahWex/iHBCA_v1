#!/usr/bin/env Rscript
# Render the LBridge differential-abundance beeswarm from a cached data RDS.
#
# Reads pre-computed plot_data_<contrast>.rds from prepare_fig4ab_lbridge_beeswarm.R
# and produces the PDF. Separating data loading from rendering keeps aesthetic
# iteration fast (~15 s) without re-reading the large upstream DA files.
#
# Plot structure
#   Y-axis: LBridge family groups ordered by compartment (Epithelial / Immune /
#     Stromal), then by median signed logFC of joint rows within each family.
#     Joint (FLEX+Xenium) rows appear first in each family group; FLEX-only L2S
#     rows (indented 3 spaces) appear below.
#   X-axis: logFC from Milo edgeR GLM with study as covariate.
#   Color: logFC diverging blue-white-red (canonical expression_diverging scale)
#     on significant nhoods (SpatialFDR < 0.1); non-significant nhoods in grey.
#   Shape: Chromium (circle) vs. Xenium (cross) in joint rows.
#   Separators: solid lines between LBridge families; heavier solid line between
#     compartments; dashed line between joint and FLEX L2S rows within a family.
#   Adaptive LFC vlines: LOESS-smoothed |logFC| vs. -log10(SpatialFDR) crossing
#     -log10(0.1); marks the approximate logFC where nhoods become significant.
#
# Outputs
#   c1  fig4a_lbridge_beeswarm[_joint_only].pdf   (Fig 4a / supplemental)
#   c3  fig4_uoq_lbridge_beeswarm[_joint_only].pdf    (Fig 4 UOQ axis / supplemental)

suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(ggplot2)
})

parser <- ArgumentParser()
parser$add_argument("--project-root", required = TRUE,
                    help = "Repository root (resolves to publication/config/)")
parser$add_argument("--plot-data",  required = TRUE,
                    help = "Path to plot_data_<contrast>.rds from prepare step")
parser$add_argument("--joint-only", action = "store_true", default = FALSE)
parser$add_argument("--out-dir",    default = NULL)
# composite-bound renders: suppress legend (so it appears once in shared legend
# panel) and override save dimensions.
parser$add_argument("--no-legend",  action = "store_true", default = FALSE)
parser$add_argument("--width-in",   type = "double", default = NULL,
                    help = "Override panel width in inches (bypasses save_panel default)")
parser$add_argument("--height-in",  type = "double", default = NULL,
                    help = "Override panel height in inches (bypasses save_panel default)")
args <- parser$parse_args()

# ── Aesthetics framework ───────────────────────────────────────────────────────
# Source the canonical aesthetics loader; pull divergent color scale, theme,
# dimensions, and validation/save helpers from publication/config/aesthetics.yaml.

root <- normalizePath(args$project_root)
source(file.path(root, "publication", "config", "load_aesthetics.R"))
aes_cfg <- load_aesthetics(file.path(root, "publication", "config"))

div_scale <- aes_cfg$scales$expression_diverging

# Adaptive LFC cutoff: |logFC| where LOESS-smoothed -log10(SpatialFDR) crosses
# -log10(fdr_threshold). Uses symmetric LOESS (degree 1) on all nhoods,
# interpolating the first crossing from the left of the sorted |logFC| axis.
compute_adaptive_lfc_cutoff <- function(da, fdr_threshold = 0.1,
                                        fallback = 0.5, span = 0.35) {
  eps      <- 1e-300
  target_y <- -log10(pmax(pmin(fdr_threshold, 1), eps))
  d <- da |>
    filter(is.finite(logFC), is.finite(SpatialFDR)) |>
    mutate(abs_lfc   = abs(logFC),
           neg_log_fdr = -log10(pmax(pmin(SpatialFDR, 1), eps))) |>
    filter(is.finite(abs_lfc), is.finite(neg_log_fdr))
  if (nrow(d) < 10) return(fallback)
  tryCatch({
    fit   <- loess(neg_log_fdr ~ abs_lfc, data = d, span = span,
                   degree = 1, family = "symmetric")
    gx    <- seq(min(d$abs_lfc), max(d$abs_lfc), length.out = 800)
    gy    <- predict(fit, newdata = data.frame(abs_lfc = gx))
    keep  <- is.finite(gy)
    gx    <- gx[keep]; gy <- gy[keep]
    i     <- which(gy >= target_y)[1]
    if (is.na(i))     return(max(d$abs_lfc))
    if (i == 1)       return(gx[1])
    w <- min(max((target_y - gy[i-1]) / (gy[i] - gy[i-1] + 1e-12), 0), 1)
    gx[i-1] + w * (gx[i] - gx[i-1])
  }, error = function(e) fallback)
}

cat("[render] Loading cached plot data...\n")
d <- readRDS(normalizePath(args$plot_data))
list2env(d, envir = environment())   # joint_df, flex_df, lbridge_stats, SIG_FDR, contrast

out_dir <- if (!is.null(args$out_dir)) args$out_dir else dirname(args$plot_data)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

adaptive_lfc <- compute_adaptive_lfc_cutoff(joint_df, fdr_threshold = SIG_FDR)
cat(sprintf("[render] Adaptive LFC cutoff (FDR<%.2f): %.3f\n", SIG_FDR, adaptive_lfc))

if (args$joint_only) {
  flex_df <- flex_df[0, ]
  cat("[render] --joint-only: FLEX L2S suppressed\n")
}

# ── Combine & order ────────────────────────────────────────────────────────────

plot_df <- bind_rows(joint_df, flex_df) |>
  mutate(sig      = SpatialFDR < SIG_FDR,
         platform = factor(platform, levels = c("FLEX", "Xenium")))

type_stats <- plot_df |>
  group_by(type_label, layer, lbridge_family) |>
  summarise(med_lfc = median(if (any(sig)) logFC[sig] else logFC, na.rm = TRUE),
            .groups = "drop")

# Row ordering: within each compartment, order all rows (joint + FLEX L2S
# combined) by signed median logFC ascending so that ggplot's bottom-up
# factor convention puts the highest logFC at the TOP of each compartment
# block (most-positive at top, most-negative at bottom). Compartment blocks
# themselves are ordered by compartment_rank ascending — smallest
# compartment_med at bottom, largest at top.
# Note: this flattens the family-block grouping; FLEX L2S indented subtypes
# now interleave with joint rows by their own med_lfc rather than being
# bound to their joint family anchor. The 3-space indent prefix remains the
# visual marker for FLEX L2S rows.
type_order_df <- type_stats |>
  left_join(lbridge_stats, by = "lbridge_family") |>
  filter(!is.na(compartment_rank)) |>
  arrange(compartment_rank, med_lfc) |>
  distinct(type_label, layer, .keep_all = TRUE)

type_order_display <- if_else(
  type_order_df$layer == "flex_l2s",
  paste0("   ", type_order_df$type_label),
  type_order_df$type_label
)

plot_df <- plot_df |>
  mutate(
    type_display = if_else(layer == "flex_l2s",
                           paste0("   ", type_label),
                           type_label),
    type_display = factor(type_display, levels = type_order_display)
  ) |>
  filter(!is.na(type_display)) |>
  mutate(y_num = as.integer(type_display))

# Manual y-dodge per platform. position_jitterdodge() dodges on x, but our x
# is continuous (logFC) — there are no x-intervals to dodge along, so platforms
# overlap when only jitterdodge is used. Compute a per-row y-offset by platform
# so FLEX and Xenium sit in adjacent vertical sub-bands within each row, then
# add a small y-jitter inside each sub-band.
dodge_half <- 0.18   # half the gap between FLEX and Xenium sub-bands
plot_df <- plot_df |>
  mutate(y_dodged = y_num + dplyr::case_when(
    platform == "FLEX"   ~ -dodge_half,
    platform == "Xenium" ~  dodge_half,
    TRUE                 ~  0
  ))

# Numeric y axis lets sidebar layers (rect, text) coexist on the same scale
# as the point layers. Discrete-style row labels are restored via scale_y_continuous.
y_breaks <- seq_along(type_order_display)
y_labels <- type_order_display

n_types <- n_distinct(plot_df$type_display)
cat(sprintf("[render] Y-axis: %d rows (%d joint + %d FLEX L2S)\n",
    n_types, sum(type_order_df$layer == "joint"),
    sum(type_order_df$layer == "flex_l2s")))

# ── Separators ─────────────────────────────────────────────────────────────────

type_meta <- plot_df |>
  distinct(type_display, lbridge_family, layer) |>
  mutate(y_pos = as.integer(type_display)) |>
  arrange(y_pos)

lbridge_solid_y <- type_meta |>
  group_by(lbridge_family) |>
  summarise(y_last = max(y_pos), .groups = "drop") |>
  arrange(y_last) |>
  filter(y_last < max(y_last)) |>
  pull(y_last) + 0.5

compartment_solid_y <- type_meta |>
  left_join(select(lbridge_stats, lbridge_family, lbridge_compartment),
            by = "lbridge_family") |>
  group_by(lbridge_compartment) |>
  summarise(y_last = max(y_pos), .groups = "drop") |>
  arrange(y_last) |>
  filter(y_last < max(y_last)) |>
  pull(y_last) + 0.5

joint_flex_dashed_y <- type_meta |>
  group_by(lbridge_family) |>
  filter(any(layer == "joint") & any(layer == "flex_l2s")) |>
  summarise(y_joint_bottom = min(y_pos[layer == "joint"]), .groups = "drop") |>
  pull(y_joint_bottom) - 0.5

# Restrict x-axis ticks to integer values within the data range. Compartment
# grouping is implicit in the row ordering (rows sorted by compartment, then
# by median signed logFC) — no separate compartment visual element is drawn.
data_lfc_max <- ceiling(max(abs(plot_df$logFC), na.rm = TRUE))
x_breaks <- seq(-data_lfc_max, data_lfc_max, by = 1)

# ── Plot ───────────────────────────────────────────────────────────────────────

p <- ggplot(plot_df, aes(y = y_num, x = logFC,
                         shape = platform, group = platform)) +
  geom_vline(xintercept = 0, color = "grey50", linewidth = 0.35) +
  geom_vline(xintercept = c(-adaptive_lfc, adaptive_lfc),
             color = "grey50", linewidth = 0.3, linetype = "dashed") +
  geom_point(
    data = . %>% filter(!sig),
    aes(y = y_dodged),
    color = "grey60",
    position = position_jitter(width = 0, height = 0.10, seed = 1L),
    size = 0.20, alpha = 0.18
  ) +
  geom_point(
    data = . %>% filter(sig),
    aes(y = y_dodged, color = logFC),
    position = position_jitter(width = 0, height = 0.10, seed = 1L),
    size = 0.61, alpha = 0.55
  ) +
  scale_color_gradient2(
    low = div_scale$low, mid = div_scale$mid, high = div_scale$high,
    midpoint = 0, name = "logFC",
    limits = c(-2.5, 2.5), oob = scales::squish,
    guide = guide_colorbar(order = 1)
  ) +
  scale_shape_manual(
    values = c(FLEX = 16, Xenium = 4),
    name   = "Platform\n(joint rows)",
    labels = c(FLEX = "Chromium", Xenium = "Xenium"),
    guide  = guide_legend(order = 2)
  ) +
  scale_y_continuous(breaks = y_breaks, labels = y_labels,
                     expand = expansion(add = 0.5)) +
  scale_x_continuous(breaks = x_breaks) +
  coord_cartesian(clip = "off") +
  get_theme(aes_cfg) +
  theme(
    axis.text.y = element_text(size = aes_cfg$typography$min_pt + 0.5, hjust = 1),
    panel.grid  = element_blank(),
    plot.margin = margin(t = 5, r = 10, b = 5, l = 35, unit = "mm")
  ) +
  labs(x = "logFC", y = NULL)

if (isTRUE(args$no_legend)) {
  p <- p + theme(legend.position = "none")
}

# Narrow-mode margin scaling: when explicit width is supplied AND it's narrower
# than the default panel width, shrink margins proportionally. The default 35mm
# left margin is sized for the 7.2-inch wide standard render; at 2.5 inches it
# would consume 55% of the panel. Scale linearly with width below 5 inches.
if (!is.null(args$width_in) && args$width_in < 5.0) {
  scale_factor <- max(0.20, args$width_in / 7.2)
  p <- p + theme(plot.margin = margin(
    t = 2,
    r = round(10 * scale_factor, 1),
    b = 2,
    l = round(35 * scale_factor, 1),
    unit = "mm"
  ))
  cat(sprintf("[render] narrow-mode margins: l=%.1fmm r=%.1fmm (scale=%.2f)\n",
              35 * scale_factor, 10 * scale_factor, scale_factor))
}

# ── Save ───────────────────────────────────────────────────────────────────────

fig_prefix <- switch(contrast,
  c1  = "fig4a",
  c2a = "fig4_p3_axis_a_vs_m",
  c2b = "fig4_p3_axis_p_vs_m",
  c3  = "fig4_uoq",
  c3a = "fig4_uiq",
  c3b = "fig4_loq",
  c3c = "fig4_liq",
  stop(sprintf("unknown contrast: %s", contrast))
)
suffix     <- if (args$joint_only) "_joint_only" else ""
panel_type <- if (args$joint_only) "beeswarm_compact" else "beeswarm_tall"
panel_id   <- sprintf("%s_lbridge_beeswarm%s", fig_prefix, suffix)

if (!is.null(args$width_in) && !is.null(args$height_in)) {
  # custom dimensions path: bypass save_panel, emit PDF + PNG directly
  pdf_path <- file.path(out_dir, paste0(panel_id, ".pdf"))
  png_path <- file.path(out_dir, paste0(panel_id, ".png"))
  ggsave(pdf_path, p, width = args$width_in, height = args$height_in,
          units = "in", device = cairo_pdf, dpi = 1200)
  ggsave(png_path, p, width = args$width_in, height = args$height_in,
          units = "in", dpi = 300)
  cat(sprintf("[save_panel] custom dims %s: PDF=%s PNG=%s (%.2fx%.2f in)\n",
              panel_id, pdf_path, png_path, args$width_in, args$height_in))
} else {

validate_panel(p, panel_type = panel_type, n_groups = n_types, config = aes_cfg)
save_panel(p, panel_id = panel_id, panel_type = panel_type,
           output_dir = out_dir, config = aes_cfg)
}
