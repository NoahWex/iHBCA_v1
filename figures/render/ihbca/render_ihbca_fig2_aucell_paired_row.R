#!/usr/bin/env Rscript
# Single-row paired boxplots for Fig 2 main C.
#
# Layout: 4 facet panels in a single row (one per signature) × 4 fibro L2 groups
# on x-axis × paired AR|BR1 boxes per L2.
#
# Pattern source: 17_aucell_l2_welch.R (same per-(panel, L2) Welch + Nee filter
# + boxplot+swarm aesthetic; differs in faceting — single row, AR|BR1 as paired
# dodge instead of separate x positions).
#
# Output:
#   <out-stem>.pdf       7.2 × 1.7 in
#   <out-stem>_stats.csv per-(panel, L2) Welch + Wilcoxon + BH stats

suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggbeeswarm)
})

parser <- ArgumentParser()
parser$add_argument("--score-dir", required = TRUE)
parser$add_argument("--publication-root", required = TRUE)
parser$add_argument("--out-dir", required = TRUE)
parser$add_argument("--panels-file", required = TRUE,
                    help = "File with panel keys, one per line (4 expected for main C)")
parser$add_argument("--donor-meta",
                    help = "Optional donor metadata CSV (ihbca_donor_id, study, ...) used to map shape = study on the per-donor dots. If omitted, shape is uniform.")
parser$add_argument("--out-stem", default = "fig2_ihbca_aucell_l2_focused_paired_row")
parser$add_argument("--cell-count-frac", type = "double", default = 0.10)
parser$add_argument("--pdf-width",  type = "double", default = 7.2)
parser$add_argument("--pdf-height", type = "double", default = 1.7)
parser$add_argument("--nrow", type = "integer", default = 1,
                    help = "Facet rows. nrow=1 = single row of paired box panels; nrow>1 = facet_wrap")
parser$add_argument("--show-pvalues", action = "store_true", default = FALSE,
                    help = "If set, prints BH-adjusted p_welch above each box pair")
parser$add_argument("--show-asterisks", action = "store_true", default = FALSE,
                    help = "If set, prints * (<0.05) / ** (<0.01) / *** (<0.001) above each box pair")
args <- parser$parse_args()

source(file.path(args$publication_root, "publication", "config",
                 "load_aesthetics.R"))
aes_cfg <- load_aesthetics(file.path(args$publication_root, "publication",
                                     "config"))
panel_theme <- get_theme(aes_cfg)
div_scale   <- aes_cfg$scales$expression_diverging
AR_COLOR <- div_scale$low
BR_COLOR <- div_scale$high

PANELS_ALL <- list(
  wu2020              = "Wu 2020",
  ohlund2017          = "Öhlund 2017",
  elyada2019_icaf     = "Elyada iCAF",
  galbo2021_pan_icaf  = "Galbo iCAF",
  elyada2019_mycaf    = "Elyada myCAF",
  galbo2021_pan_mycaf = "Galbo myCAF",
  nee2023_precaf      = "Nee preCAF",
  reed2024_fb2        = "Reed FB2",
  mmp_convergent      = "MMP3/10/12",
  senescence_sasp     = "Sen SASP",
  senescence_cca      = "Sen CCA"
)
keep <- trimws(readLines(args$panels_file))
keep <- keep[nzchar(keep)]
miss <- setdiff(keep, names(PANELS_ALL))
if (length(miss)) stop("Unknown panels: ", paste(miss, collapse = ", "))
PANELS <- PANELS_ALL[keep]
cat(sprintf("[panels] %d active: %s\n", length(PANELS),
            paste(names(PANELS), collapse = ", ")))

FIBROS <- c("str__Fibro-major", "str__Fibro-IGF1",
            "str__Fibro-SFRP4", "str__Fibro-prematrix")

if (!dir.exists(args$out_dir)) dir.create(args$out_dir, recursive = TRUE)

# --- Load + filter + per-(panel, L2) Welch + Wilcoxon ---
all_rows <- list()
stats_rows <- list()
for (pkey in names(PANELS)) {
  csv_path <- file.path(args$score_dir, paste0(pkey, "_per_donor_scores.csv"))
  if (!file.exists(csv_path)) { cat(sprintf("[skip] %s\n", csv_path)); next }
  per_donor <- read_csv(csv_path, show_col_types = FALSE) |>
    filter(risk_class %in% c("AR", "BR1"),
           L2_label %in% FIBROS) |>
    mutate(panel = pkey) |>
    group_by(panel, L2_label) |>
    mutate(group_mean = mean(n_cells),
           cell_thr   = args$cell_count_frac * group_mean,
           kept       = n_cells >= cell_thr) |>
    ungroup()
  all_rows[[pkey]] <- per_donor
  for (l2 in FIBROS) {
    sub <- per_donor |> filter(L2_label == l2, kept)
    ar  <- sub |> filter(risk_class == "AR")  |> pull(mean_AUC)
    br1 <- sub |> filter(risk_class == "BR1") |> pull(mean_AUC)
    if (length(ar) < 2 || length(br1) < 2) {
      stats_rows[[length(stats_rows)+1L]] <- data.frame(
        panel = pkey, L2 = l2, n_AR = length(ar), n_BR1 = length(br1),
        delta = NA, p_welch = NA, p_wilcox = NA, p_glm_study = NA); next
    }
    w <- t.test(br1, ar, var.equal = FALSE)
    u <- wilcox.test(br1, ar)
    # GLM with study covariate (Wald test on risk_classBR1 coefficient).
    # Only fit if `study` is in the score CSV AND >1 study level present
    # in the kept donors (else lm drops the covariate and equals Welch).
    p_glm <- NA_real_
    if ("study" %in% names(sub) && length(unique(sub$study)) >= 2) {
      sub_lm <- sub |> mutate(risk_class = factor(risk_class, levels = c("AR", "BR1")))
      fit <- tryCatch(lm(mean_AUC ~ risk_class + study, data = sub_lm),
                      error = function(e) NULL)
      if (!is.null(fit)) {
        coefs <- summary(fit)$coefficients
        if ("risk_classBR1" %in% rownames(coefs)) {
          p_glm <- coefs["risk_classBR1", "Pr(>|t|)"]
        }
      }
    }
    stats_rows[[length(stats_rows)+1L]] <- data.frame(
      panel = pkey, L2 = l2, n_AR = length(ar), n_BR1 = length(br1),
      mean_AR  = mean(ar), mean_BR1 = mean(br1),
      delta = mean(br1) - mean(ar),
      p_welch = w$p.value, p_wilcox = u$p.value, p_glm_study = p_glm)
  }
}
df <- bind_rows(all_rows)
stats <- bind_rows(stats_rows) |>
  group_by(panel) |>
  mutate(p_welch_bh      = p.adjust(p_welch,      method = "BH"),
         p_wilcox_bh     = p.adjust(p_wilcox,     method = "BH"),
         p_glm_study_bh  = p.adjust(p_glm_study,  method = "BH")) |>
  ungroup()
write_csv(stats, file.path(args$out_dir, paste0(args$out_stem, "_stats.csv")))
cat("\nPer-(panel, L2) stats (BH within panel × 4-L2 family):\n")
print(stats |> select(panel, L2, n_AR, n_BR1, delta, p_welch, p_welch_bh))

# --- Plot prep ---
df_kept <- df |> filter(kept)
panel_keys   <- names(PANELS)
panel_labels <- unname(unlist(PANELS))
df_kept$panel <- factor(df_kept$panel, levels = panel_keys, labels = panel_labels)
df_kept$L2_short <- factor(sub("^str__", "", df_kept$L2_label),
                            levels = c("Fibro-major", "Fibro-IGF1",
                                       "Fibro-SFRP4", "Fibro-prematrix"),
                            labels = c("Major", "IGF1", "SFRP4", "Pre"))
df_kept$risk_class <- factor(df_kept$risk_class, levels = c("AR", "BR1"))

# --- v3: shape-by-study mapping ---
# Prefer the `study` column emitted directly by the AUCell builder (post-2026-05-22:
# build_ihbca_fig2_aucell_scores.py joins donor metadata at scoring time and
# emits study into per_donor_scores.csv). Fall back to legacy --donor-meta join
# if the score CSV is from an older build.
use_study_shape <- FALSE
STUDY_LEVELS <- c("gray", "kumar", "murrow", "nee", "pal", "reed", "twigger")

if ("study" %in% names(df_kept)) {
  n_with_study <- sum(!is.na(df_kept$study))
  cat(sprintf("[study-shape] study column already in score CSV: %d / %d rows have study\n",
              n_with_study, nrow(df_kept)))
  if (n_with_study > 0) {
    # Order by canonical 7-study sequence, then drop levels not present in this
    # panel's kept donors so the shape legend only lists represented studies.
    df_kept$study <- droplevels(factor(df_kept$study, levels = STUDY_LEVELS))
    cat(sprintf("[study-shape] %d studies represented after droplevels: %s\n",
                nlevels(df_kept$study),
                paste(levels(df_kept$study), collapse = ", ")))
    use_study_shape <- TRUE
  }
} else if (!is.null(args$donor_meta) && nzchar(args$donor_meta)) {
  # Legacy fallback: join donor metadata if score CSV doesn't carry study
  donor_meta <- read_csv(args$donor_meta, show_col_types = FALSE) |>
    select(ihbca_donor_id, study) |>
    distinct()
  df_kept <- df_kept |>
    left_join(donor_meta, by = c("patientID" = "ihbca_donor_id"))
  n_with_study <- sum(!is.na(df_kept$study))
  cat(sprintf("[study-shape] legacy join via --donor-meta: %d / %d rows have study\n",
              n_with_study, nrow(df_kept)))
  if (n_with_study > 0) {
    df_kept$study <- droplevels(factor(df_kept$study, levels = STUDY_LEVELS))
    cat(sprintf("[study-shape] %d studies represented after droplevels: %s\n",
                nlevels(df_kept$study),
                paste(levels(df_kept$study), collapse = ", ")))
    use_study_shape <- TRUE
  }
}

# --- Optional sig brackets (asterisks or numeric p_adj) ---
sig_layer <- list()
if (args$show_pvalues || args$show_asterisks) {
  brackets <- df_kept |>
    group_by(panel, L2_short) |>
    summarise(y_top = max(mean_AUC), y_range = diff(range(mean_AUC)),
              .groups = "drop") |>
    mutate(y_text = y_top + 0.12 * y_range)
  stats_lbl <- stats |>
    mutate(panel = factor(panel, levels = panel_keys, labels = panel_labels),
           L2_short = factor(sub("^str__", "", L2),
                              levels = c("Fibro-major", "Fibro-IGF1",
                                         "Fibro-SFRP4", "Fibro-prematrix"),
                              labels = c("Major", "IGF1", "SFRP4", "Pre")),
           label = if (args$show_asterisks) {
             dplyr::case_when(
               is.na(p_welch_bh)   ~ "",
               p_welch_bh < 0.001  ~ "***",
               p_welch_bh < 0.01   ~ "**",
               p_welch_bh < 0.05   ~ "*",
               TRUE                 ~ "")
           } else {
             sprintf("%.3f", p_welch_bh)
           }) |>
    left_join(brackets, by = c("panel", "L2_short"))
  txt_size <- if (args$show_asterisks) 2.4 else 1.6
  sig_layer <- list(
    geom_text(data = stats_lbl,
              aes(x = L2_short, y = y_text, label = label),
              inherit.aes = FALSE, size = txt_size, color = "grey15",
              fontface = "bold")
  )
}

# --- Single-row paired boxplots ---
DODGE <- 0.65

# v2: shape=study layered onto the quasirandom dots if donor_meta was supplied.
# Shape values are color-only (not fill-aware), so dots inherit color from
# risk_class color scale and shape from the study scale.
#
# Explicit group=risk_class avoids ggplot's auto-group inferring interaction(
# risk_class, study), which would fragment the dodge per-study and visually
# interleave AR/BR1 dots within a single L2 column.
quasirandom_layer <- if (use_study_shape) {
  geom_quasirandom(aes(shape = study, group = risk_class),
                   dodge.width = DODGE, size = 0.55,
                   alpha = 0.85, width = 0.16, stroke = 0.25)
} else {
  geom_quasirandom(dodge.width = DODGE, size = 0.35, alpha = 0.75,
                   width = 0.16)
}
STUDY_SHAPES <- c(gray = 16, kumar = 17, murrow = 15, nee = 5,
                  pal = 6, reed = 4, twigger = 3)
shape_scale_layer <- if (use_study_shape) {
  scale_shape_manual(values = STUDY_SHAPES, name = "Study", drop = FALSE)
} else {
  NULL
}

p <- ggplot(df_kept,
            aes(x = L2_short, y = mean_AUC, fill = risk_class, color = risk_class)) +
  geom_boxplot(width = 0.55, alpha = 0.40, outlier.shape = NA,
               linewidth = 0.25, position = position_dodge(width = DODGE)) +
  quasirandom_layer +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 1.2,
               color = "black", position = position_dodge(width = DODGE)) +
  sig_layer +
  scale_fill_manual(values = c(AR = AR_COLOR, BR1 = BR_COLOR),
                    name = NULL) +
  scale_color_manual(values = c(AR = AR_COLOR, BR1 = BR_COLOR), guide = "none") +
  shape_scale_layer +
  facet_wrap(~ panel, nrow = args$nrow, scales = "free_y") +
  labs(x = NULL, y = "AUC (per donor)") +
  panel_theme +
  theme(strip.background = element_blank(),
        strip.text       = element_text(face = "bold", size = 7),
        panel.grid       = element_blank(),
        panel.background = element_blank(),
        axis.text.x      = element_text(size = 6, angle = 35, hjust = 1),
        axis.text.y      = element_text(size = 6),
        axis.title.y     = element_text(size = 7),
        legend.position  = "top",
        legend.text      = element_text(size = 6),
        legend.key.size  = unit(0.30, "cm"),
        legend.margin    = margin(0, 0, 0, 0),
        panel.spacing.x  = unit(0.4, "lines"),
        plot.margin      = margin(2, 4, 2, 4))

pdf_path <- file.path(args$out_dir, paste0(args$out_stem, ".pdf"))
ggsave(pdf_path, p, width = args$pdf_width, height = args$pdf_height,
       units = "in", device = cairo_pdf)
cat(sprintf("\nWrote %s (%.2f x %.2f in; %d panels × 4 fibro L2 × paired AR|BR1)\n",
            pdf_path, args$pdf_width, args$pdf_height, length(PANELS)))
