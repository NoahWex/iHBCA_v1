#!/usr/bin/env Rscript
# AUCell L1 fibroblast aggregation, Welch's t-test on per-donor mean AUC.
#
# Methodology matches Nee et al. 2023 (Nat Genet 55:595-606): per-donor mean
# AUC, Welch's two-sample t-test, low-cell-count donor exclusion (<10% of
# group mean). The methods-QC claim is "atlas preserves prior-reported BRCA1-
# stromal fibroblast signatures from the source studies whose donors compose
# the atlas." Method-aligning with the source paper makes the comparison crisp.
#
# Inputs (already on disk from scripts/12_aucell_preneoplastic):
#   - 5 per_cell_scores.csv (cell_id, patientID, L2_label, risk_class, AUC)
#
# Operation:
#   1. Aggregate cells per (donor, panel) -> mean AUC + n_cells
#      (L1 = Fibroblast: pool the 4 fibro L2s)
#   2. Apply Nee cell-count filter: drop donors with n_cells < 0.10 * group_mean
#   3. Welch's two-sample t-test BR1 vs AR per panel
#   4. Render Fig 2-styled box+swarm: AR=expression_diverging.low,
#      BR1=expression_diverging.high; mean+SD overlay; significance bracket
#
# Outputs:
#   l1_per_donor_summary.csv  (donor x panel mean AUC + n_cells, kept flag)
#   l1_stats.csv              (panel x Welch + Wilcoxon + filter audit)
#   fibro_l1_beeswarm.pdf     (Fig 2-styled box+swarm)
#
# Code pattern provenance:
#   - ggbeeswarm idiom: publication/figures/render/ihbca/render_fig2_beeswarm_faceted.R:246-258
#   - AR/BR1 anchors: publication/figures/render/ihbca/render_fig2_split_forest.R:72-76
#   - Welch + cell filter: matches Nee et al. 2023 Nat Genet supp methods

suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(ggbeeswarm)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

parser <- ArgumentParser()
parser$add_argument("--score-dir", required = TRUE,
                    help = "Directory containing {panel}_per_cell_scores.csv from script 12")
parser$add_argument("--publication-root", required = TRUE)
parser$add_argument("--out-dir", required = TRUE)
parser$add_argument("--cell-count-frac", type = "double", default = 0.10,
                    help = "Drop donors with n_cells < frac * mean(n_cells) per panel (Nee 2023 rule = 0.10)")
parser$add_argument("--panels-file", default = "",
                    help = "Optional file with one panel key per line (overrides default 11-panel set)")
parser$add_argument("--out-stem", default = "fibro_l1_beeswarm",
                    help = "Output filename stem (without .pdf)")
args <- parser$parse_args()

source(file.path(args$publication_root, "publication", "config",
                 "load_aesthetics.R"))
aes_cfg <- load_aesthetics(file.path(args$publication_root, "publication",
                                     "config"))
panel_theme <- get_theme(aes_cfg)
div_scale   <- aes_cfg$scales$expression_diverging

AR_COLOR <- div_scale$low   # "#2166AC"
BR_COLOR <- div_scale$high  # "#B2182B"

PANELS <- list(
  # Tumor-context anchors
  wu2020              = list(label = "Wu 2020\n(breast tumor)",       group = "tumor_iCAF"),
  ohlund2017          = list(label = "Öhlund 2017\n(pancreas)",       group = "tumor_iCAF"),
  elyada2019_icaf     = list(label = "Elyada 2019\niCAF (PDAC)",      group = "tumor_iCAF"),
  galbo2021_pan_icaf  = list(label = "Galbo 2021\npan iCAF",          group = "tumor_iCAF"),
  elyada2019_mycaf    = list(label = "Elyada 2019\nmyCAF (PDAC)",     group = "tumor_myCAF"),
  galbo2021_pan_mycaf = list(label = "Galbo 2021\npan myCAF",         group = "tumor_myCAF"),
  # Internal reproduction (source studies in iHBCA)
  nee2023_precaf      = list(label = "Nee 2023\npreCAF",              group = "internal_repro"),
  reed2024_fb2        = list(label = "Reed 2024\nFB2",                group = "internal_repro"),
  mmp_convergent      = list(label = "MMP3/10/12",                    group = "internal_repro"),
  # Senescence (Suzuki 2024 SenNet biomarker DB)
  senescence_sasp     = list(label = "Senescence\nSASP",              group = "senescence"),
  senescence_cca      = list(label = "Senescence\nCCA",               group = "senescence")
)

# Optional subset of panels via --panels-file
if (nchar(args$panels_file) > 0L && file.exists(args$panels_file)) {
  keep <- readLines(args$panels_file) |> trimws()
  keep <- keep[nzchar(keep)]
  PANELS <- PANELS[keep[keep %in% names(PANELS)]]
  cat(sprintf("[subset] %d panels kept: %s\n",
              length(PANELS), paste(names(PANELS), collapse = ", ")))
}

if (!dir.exists(args$out_dir)) dir.create(args$out_dir, recursive = TRUE)

# --- Load per-cell AUC, collapse to L1 with cell-count filter ---
per_donor_long <- list()
filter_audit <- list()

for (pkey in names(PANELS)) {
  csv_path <- file.path(args$score_dir, paste0(pkey, "_per_cell_scores.csv"))
  if (!file.exists(csv_path)) stop("Missing per-cell scores: ", csv_path)
  per_cell <- read_csv(csv_path, show_col_types = FALSE)
  per_cell <- per_cell |> filter(risk_class %in% c("AR", "BR1"))

  # Per donor: mean AUC + n_cells across all fibroblast L2s (= L1=Fibroblast)
  per_donor <- per_cell |>
    group_by(patientID, risk_class) |>
    summarise(AUC     = mean(AUC, na.rm = TRUE),
              n_cells = dplyr::n(),
              .groups = "drop") |>
    mutate(panel = pkey)

  # Nee 2023 cell-count filter: drop donors with n_cells < frac * group_mean
  group_mean <- mean(per_donor$n_cells)
  cell_thr <- args$cell_count_frac * group_mean
  per_donor$kept <- per_donor$n_cells >= cell_thr
  n_drop <- sum(!per_donor$kept)
  cat(sprintf("[%s] %d donors total; mean cells/donor=%.0f; threshold=%.0f (%.0f%% of mean); dropped %d\n",
              pkey, nrow(per_donor), group_mean, cell_thr,
              100 * args$cell_count_frac, n_drop))

  filter_audit[[pkey]] <- list(
    panel = pkey,
    n_donors_total = nrow(per_donor),
    n_donors_kept  = sum(per_donor$kept),
    n_donors_dropped = n_drop,
    group_mean_cells = group_mean,
    cell_threshold   = cell_thr
  )

  per_donor_long[[pkey]] <- per_donor
}
df <- bind_rows(per_donor_long)
write_csv(df, file.path(args$out_dir, "l1_per_donor_summary.csv"))

# --- Statistics: Welch's t-test (primary, Nee 2023) + Wilcoxon (sensitivity) ---
stats_rows <- lapply(names(PANELS), function(pkey) {
  pan <- df |> filter(panel == pkey, kept)
  ar  <- pan |> filter(risk_class == "AR")  |> pull(AUC)
  br1 <- pan |> filter(risk_class == "BR1") |> pull(AUC)
  fa <- filter_audit[[pkey]]

  if (length(ar) < 2 || length(br1) < 2) {
    return(data.frame(panel = pkey, n_AR = length(ar), n_BR1 = length(br1),
                      mean_AR = NA, mean_BR1 = NA, sd_AR = NA, sd_BR1 = NA,
                      delta_mean = NA, p_welch = NA, p_wilcox = NA,
                      n_donors_dropped = fa$n_donors_dropped,
                      group_mean_cells = fa$group_mean_cells,
                      cell_threshold   = fa$cell_threshold))
  }

  welch  <- t.test(br1, ar, var.equal = FALSE)
  wilcox <- wilcox.test(br1, ar)

  data.frame(
    panel = pkey,
    n_AR  = length(ar), n_BR1 = length(br1),
    mean_AR = mean(ar), mean_BR1 = mean(br1),
    sd_AR   = sd(ar),   sd_BR1   = sd(br1),
    delta_mean = mean(br1) - mean(ar),
    p_welch  = welch$p.value,
    p_wilcox = wilcox$p.value,
    n_donors_dropped = fa$n_donors_dropped,
    group_mean_cells = fa$group_mean_cells,
    cell_threshold   = fa$cell_threshold
  )
})
stats <- bind_rows(stats_rows)
stats$p_welch_bh  <- p.adjust(stats$p_welch,  method = "BH")
stats$p_wilcox_bh <- p.adjust(stats$p_wilcox, method = "BH")
write_csv(stats, file.path(args$out_dir, "l1_stats.csv"))

cat("\nL1=Fibroblast statistics (Welch primary, Wilcoxon sensitivity):\n")
print(stats)

# --- Render: Fig 2-styled box + swarm + significance bracket ---
df_kept <- df |> filter(kept)
panel_order  <- names(PANELS)
panel_labels <- vapply(PANELS, function(x) x$label, character(1))
df_kept$panel <- factor(df_kept$panel, levels = panel_order, labels = panel_labels)
df_kept$risk_class <- factor(df_kept$risk_class, levels = c("AR", "BR1"))

# Significance bracket positions: top of each facet's data + small offset
brackets <- df_kept |>
  group_by(panel) |>
  summarise(y_top = max(AUC, na.rm = TRUE),
            y_range = diff(range(AUC, na.rm = TRUE)),
            .groups = "drop") |>
  mutate(y_bar  = y_top + 0.10 * y_range,
         y_tick = y_top + 0.08 * y_range,
         y_text = y_top + 0.16 * y_range)

stats_lbl <- stats |>
  mutate(panel = factor(panel, levels = panel_order, labels = panel_labels),
         label = ifelse(p_welch_bh < 0.001,
                        sprintf("p_adj = %.1e", p_welch_bh),
                        sprintf("p_adj = %.3f", p_welch_bh))) |>
  select(panel, label) |>
  left_join(brackets, by = "panel")

p <- ggplot(df_kept, aes(x = risk_class, y = AUC)) +
  # Box with AR/BR1 fills
  geom_boxplot(aes(fill = risk_class), width = 0.55, alpha = 0.45,
               outlier.shape = NA, linewidth = 0.3, color = "grey20") +
  # Per-donor points overlay
  geom_quasirandom(aes(color = risk_class), size = 1.4, alpha = 0.85,
                   width = 0.2) +
  # Mean + SD error bars (matches Nee 2023 "data as mean ± s.d.")
  stat_summary(fun = mean, geom = "point", shape = 18, size = 2.6,
               color = "black") +
  stat_summary(fun.data = function(y) {
                 m <- mean(y); s <- sd(y)
                 data.frame(y = m, ymin = m - s, ymax = m + s)
               },
               geom = "errorbar", width = 0.25, linewidth = 0.4,
               color = "black") +
  # Significance bracket
  geom_segment(data = stats_lbl,
               aes(x = 1, xend = 2, y = y_bar, yend = y_bar),
               inherit.aes = FALSE, linewidth = 0.3, color = "grey20") +
  geom_segment(data = stats_lbl,
               aes(x = 1, xend = 1, y = y_tick, yend = y_bar),
               inherit.aes = FALSE, linewidth = 0.3, color = "grey20") +
  geom_segment(data = stats_lbl,
               aes(x = 2, xend = 2, y = y_tick, yend = y_bar),
               inherit.aes = FALSE, linewidth = 0.3, color = "grey20") +
  geom_text(data = stats_lbl,
            aes(x = 1.5, y = y_text, label = label),
            inherit.aes = FALSE, size = 2.0, color = "grey20") +
  scale_fill_manual(values  = c(AR = AR_COLOR, BR1 = BR_COLOR), guide = "none") +
  scale_color_manual(values = c(AR = AR_COLOR, BR1 = BR_COLOR), guide = "none") +
  facet_wrap(~ panel, nrow = 2, scales = "free_y") +
  labs(x = NULL, y = "AUC (per donor)") +
  panel_theme +
  theme(strip.background = element_blank(),
        strip.text       = element_text(face = "bold", size = 7),
        panel.grid       = element_blank(),
        panel.background = element_blank(),
        axis.text.x      = element_text(size = 7),
        axis.text.y      = element_text(size = 6),
        axis.title.y     = element_text(size = 7),
        panel.spacing.x  = unit(0.4, "lines"))

validate_panel(p, panel_type = "beeswarm_compact",
               n_groups = length(PANELS) * 2, config = aes_cfg)

pdf_path <- file.path(args$out_dir, paste0(args$out_stem, ".pdf"))
ggsave(pdf_path, plot = p, width = 9.5, height = 5.0,
       units = "in", device = cairo_pdf)
cat(sprintf("\nWrote %s (7.2 x 2.8 in)\n", pdf_path))
