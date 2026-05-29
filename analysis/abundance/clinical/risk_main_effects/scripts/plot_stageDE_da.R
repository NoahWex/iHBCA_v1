#!/usr/bin/env Rscript
# plot_stageDE_da.R
# Hero + supplemental beeswarms from Stage D outputs.
#
# Layout (rotated): y = L2 (faceted by compartment, rows); x = logFC.
# Within each L2 row, swarms are vertically dodged by contrast.
#
# Outputs to: reports/figures/
#   da_beeswarm_main.pdf            — AR / HR_sporadic / HR_BRCA1 (subgroup effects)
#   da_beeswarm_interactions.pdf    — parity_x_HR_sporadic / BRCA1 / BRCA2 (modifier tests)
#   da_beeswarm_supplemental.pdf    — HR_BRCA2 (sensitivity)
#   da_beeswarm_per_study.pdf       — main contrasts dodged by study (transparency)
#   da_volcano_<contrast>.pdf       — per-contrast volcano
#   per_l2_summary_heatmap.pdf      — L2 × contrast n_sig (Stage E)

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(yaml)
  library(scales); library(tidyr)
})
have_beeswarm <- requireNamespace("ggbeeswarm", quietly = TRUE)

`%||%` <- function(a, b) if (is.null(a)) b else a

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "plot_stageDE_da.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths
inquiry <- cfg$inquiry

fig_dir <- file.path(paths$outputs$reports, "figures")
ensure_dir(fig_dir)

cat(sprintf("=== Plot Stage D/E DA results for '%s' ===\n", inquiry$inquiry))
cat(sprintf("ggbeeswarm available: %s\n", have_beeswarm))

contrast_dirs <- list.dirs(paths$outputs$stageD, recursive = FALSE)
if (length(contrast_dirs) == 0) {
  cat("No Stage D outputs yet.\n"); quit(status = 0)
}

# Load all contrasts ----------------------------------------------------------
all_da <- list()
for (cdir in contrast_dirs) {
  cname <- basename(cdir)
  da_csv <- file.path(cdir, "da_results.csv")
  if (!file.exists(da_csv)) next
  da <- read_csv(da_csv, show_col_types = FALSE)
  fdr_t <- 0.05
  rs_path <- file.path(cdir, "run_summary.yaml")
  if (file.exists(rs_path)) fdr_t <- yaml::read_yaml(rs_path)$spatial_fdr %||% 0.05
  comp_col <- intersect(c("compartment", "L2_compartment"), colnames(da))[1]
  lab_col  <- intersect(c("label", "L2_label"), colnames(da))[1]
  if (is.na(comp_col) || is.na(lab_col)) {
    cat(sprintf("  %s: missing compartment/label cols, skipping\n", cname)); next
  }
  da$contrast <- cname
  da$fdr_t <- fdr_t
  da <- da %>%
    rename(L2_compartment = all_of(comp_col), L2_label = all_of(lab_col)) %>%
    filter(!is.na(L2_compartment), !is.na(L2_label)) %>%
    mutate(sig = !is.na(SpatialFDR) & SpatialFDR < fdr_t,
           neglog_fdr = -log10(pmax(SpatialFDR, 1e-300)))
  all_da[[cname]] <- da
}

if (length(all_da) == 0) {
  cat("No usable contrast outputs.\n"); quit(status = 0)
}

cat(sprintf("Loaded %d contrasts: %s\n",
            length(all_da), paste(names(all_da), collapse = ", ")))

# Build dodged beeswarm -------------------------------------------------------
# `contrast_set`: ordered vector of contrast names to render; top contrast plots
# at bottom row of the dodge (ggplot y is reversed for clarity).
# `drop_artifacts`: when TRUE, exclude L2 labels flagged is_artifact (Doublet, SS, DS, etc.)
build_beeswarm <- function(contrast_set, fig_path,
                            main_title, sub_title,
                            contrast_palette = NULL,
                            drop_artifacts = TRUE) {
  available <- intersect(contrast_set, names(all_da))
  if (length(available) < 2) {
    cat(sprintf("  skip %s: need ≥2 of %s, have %s\n",
                fig_path, paste(contrast_set, collapse = ","),
                paste(available, collapse = ",")))
    return(invisible(NULL))
  }

  d <- bind_rows(all_da[available])
  d$contrast <- factor(d$contrast, levels = available)
  d$lfc_capped <- pmin(pmax(d$logFC, -3), 3)

  # Artifact filter — operates on da_results' is_artifact column (canonical
  # publication source: annotation_v2_*.yaml flags Doublet / SS / DS / *_Doublets).
  if (drop_artifacts && "is_artifact" %in% colnames(d)) {
    art <- as.logical(d$is_artifact)
    art[is.na(art)] <- FALSE
    d <- d[!art, , drop = FALSE]
  }

  # Order L2 within each compartment by the FIRST contrast's median logFC
  base <- d %>% filter(contrast == available[1])
  l2_order <- base %>%
    group_by(L2_compartment, L2_label) %>%
    summarise(med = median(logFC, na.rm = TRUE), .groups = "drop") %>%
    arrange(L2_compartment, med) %>%
    mutate(L2_full = paste(L2_compartment, L2_label, sep = "::"))
  d$L2_full <- paste(d$L2_compartment, d$L2_label, sep = "::")
  d$L2_full <- factor(d$L2_full, levels = l2_order$L2_full)

  # Per-L2 alpha modulation. Direction: small L2s get HIGHER baseline alpha so
  # sparse swarms read; large L2s get LOWER baseline alpha so dense swarms don't
  # blob. Independent of, and multiplied with, the sig/ns alpha. log10(n_nhoods)
  # is the knob; modulate around median, clamped to [0.6, 1.5].
  l2_n <- d %>% group_by(L2_full) %>% summarise(n = dplyr::n(), .groups = "drop") %>%
    mutate(log_n = log10(pmax(n, 10)))
  mid_log_n <- median(l2_n$log_n, na.rm = TRUE)
  l2_n <- l2_n %>%
    mutate(alpha_mod = pmin(pmax(1 - 0.18 * (log_n - mid_log_n), 0.6), 1.5))
  d <- d %>% left_join(l2_n %>% select(L2_full, alpha_mod), by = "L2_full")
  # sig: 0.95 base; ns: 0.15 base. Multiply by per-L2 modulation.
  d$alpha_eff <- ifelse(d$sig, 0.95, 0.15) * d$alpha_mod
  d$alpha_eff <- pmin(d$alpha_eff, 1)

  # Subsample ns nhoods aggressively; keep all sig
  d_sig <- d %>% filter(sig)
  d_ns  <- d %>% filter(!sig)
  if (nrow(d_ns) > 80000) d_ns <- d_ns[sample.int(nrow(d_ns), 80000), ]

  # Dodge by contrast within each L2 row
  swarm_layer <- if (have_beeswarm) {
    list(
      ggbeeswarm::geom_quasirandom(data = d_ns,
                                    aes(x = logFC, y = L2_full, colour = contrast,
                                        alpha = alpha_eff),
                                    groupOnX = FALSE, dodge.width = 0.85,
                                    size = 0.18, bandwidth = 0.5,
                                    show.legend = FALSE),
      ggbeeswarm::geom_quasirandom(data = d_sig,
                                    aes(x = logFC, y = L2_full,
                                        colour = contrast, group = contrast,
                                        alpha = alpha_eff),
                                    groupOnX = FALSE, dodge.width = 0.85,
                                    size = 0.45, bandwidth = 0.5)
    )
  } else {
    list(
      geom_jitter(data = d_ns, aes(x = logFC, y = L2_full, colour = contrast,
                                    alpha = alpha_eff),
                  position = position_jitterdodge(jitter.width = 0.18,
                                                   dodge.width = 0.85),
                  size = 0.18, show.legend = FALSE),
      geom_jitter(data = d_sig, aes(x = logFC, y = L2_full, colour = contrast,
                                     alpha = alpha_eff),
                  position = position_jitterdodge(jitter.width = 0.18,
                                                   dodge.width = 0.85),
                  size = 0.45)
    )
  }

  pal <- if (!is.null(contrast_palette)) contrast_palette else
    setNames(c("#444444", "#E07B39", "#9B2226", "#3A86FF",
               "#E07B39", "#9B2226", "#3A86FF"),
             c("parity_in_AR", "parity_in_HR_sporadic",
               "parity_in_HR_BRCA1", "parity_in_HR_BRCA2",
               "parity_x_HR_sporadic", "parity_x_HR_BRCA1",
               "parity_x_HR_BRCA2"))

  p <- ggplot() +
    swarm_layer +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "grey40", linewidth = 0.3) +
    scale_colour_manual(values = pal, drop = FALSE,
                        guide = guide_legend(override.aes = list(size = 2.5,
                                                                  alpha = 1))) +
    scale_alpha_identity() +
    facet_grid(L2_compartment ~ ., scales = "free_y", space = "free_y") +
    coord_cartesian(xlim = c(-4, 4)) +
    theme_minimal(base_size = 9) +
    theme(axis.text.y = element_text(size = 7),
          axis.text.x = element_text(size = 8),
          panel.grid.major.y = element_blank(),
          strip.background = element_rect(fill = "grey95", colour = NA),
          strip.text = element_text(face = "bold", size = 9),
          legend.position = "right") +
    labs(title = main_title, subtitle = sub_title,
         x = "logFC", y = NULL, colour = "contrast")

  # Width scales modestly with #contrasts (more contrasts → wider dodge)
  # Height scales with #L2 rows (each row needs vertical room for dodged swarms)
  n_l2 <- nlevels(d$L2_full)
  n_c  <- length(available)
  fig_w <- 9 + 0.7 * n_c
  fig_h <- max(7, min(22, 0.42 * n_l2 + 0.5 * n_c + 3))
  ggsave(fig_path, p, width = fig_w, height = fig_h, limitsize = FALSE)
  cat(sprintf("  wrote: %s (%d L2s, %d contrasts, %.1fx%.1f in)\n",
              fig_path, n_l2, n_c, fig_w, fig_h))
}

# MAIN: AR baseline + BRCA1 modifier + Sporadic modifier (artifact-free) ------
# Layout: per L2 row, three contrasts dodged top→bot:
#   (1) parity_in_AR        — AR main effect baseline
#   (2) parity_x_HR_BRCA1   — BRCA1 carrier × parity interaction vs AR
#   (3) parity_x_HR_sporadic— Sporadic HR × parity interaction vs AR
# Stratified within-HR contrasts (parity_in_HR_*) are shown only as volcanos —
# they do not reach significance and are documented as null-by-power.
main_contrasts <- c("parity_in_AR", "parity_x_HR_BRCA1", "parity_x_HR_sporadic")
build_beeswarm(
  contrast_set = main_contrasts,
  fig_path = file.path(fig_dir, "da_beeswarm_main.pdf"),
  main_title = "Parity DA — AR baseline vs BRCA1 / sporadic interaction modifiers",
  sub_title  = "Artifact L2s excluded. Sig nhoods opaque; ns subsampled, alpha modulated by L2 density.",
  drop_artifacts = TRUE
)

# SUPP 1: same main panel but retain artifact L2s (transparency) -------------
build_beeswarm(
  contrast_set = main_contrasts,
  fig_path = file.path(fig_dir, "da_beeswarm_main_with_artifacts.pdf"),
  main_title = "Parity DA — main contrasts including artifact L2s (Doublet/SS/DS)",
  sub_title  = "Same statistics as main panel; artifact-flagged labels retained for transparency.",
  drop_artifacts = FALSE
)

# SUPP 2: extend main with BRCA2 + pooled-HR for transparency -----------------
build_beeswarm(
  contrast_set = c(main_contrasts, "parity_x_HR_BRCA2", "parity_in_HR_BRCA12"),
  fig_path = file.path(fig_dir, "da_beeswarm_main_supp_brca2_pooled.pdf"),
  main_title = "Parity DA — main + BRCA2 modifier + pooled HR_germline",
  sub_title  = "Adds parity_x_HR_BRCA2 (smaller cohort) and parity_in_HR_BRCA12 (BRCA1+2 pooled within-HR sensitivity).",
  drop_artifacts = TRUE
)

# AR BASELINE solo (context, narrow scope) -----------------------------------
build_beeswarm(
  contrast_set = c("parity_in_AR"),
  fig_path = file.path(fig_dir, "da_beeswarm_AR_baseline.pdf"),
  main_title = "Parity DA — AR baseline (within-stratum)",
  sub_title  = "Reference for modifier interpretation. n=125 AR donors.",
  drop_artifacts = TRUE
)

# PER-STUDY TRANSPARENCY: main contrasts faceted by contrast, dodged by study
build_per_study_beeswarm <- function(contrast_set, fig_path, main_title, sub_title) {
  available <- intersect(contrast_set, names(all_da))
  if (length(available) < 1) return(invisible(NULL))
  d <- bind_rows(all_da[available])
  d$contrast <- factor(d$contrast, levels = available)
  if (!"index_study" %in% colnames(d)) {
    cat(sprintf("  skip %s: no index_study column\n", fig_path)); return(invisible(NULL))
  }
  d <- d %>% filter(!is.na(index_study))

  base <- all_da[[available[1]]]
  l2_order <- base %>%
    group_by(L2_compartment, L2_label) %>%
    summarise(med = median(logFC, na.rm = TRUE), .groups = "drop") %>%
    arrange(L2_compartment, med) %>%
    mutate(L2_full = paste(L2_compartment, L2_label, sep = "::"))
  d$L2_full <- paste(d$L2_compartment, d$L2_label, sep = "::")
  d$L2_full <- factor(d$L2_full, levels = l2_order$L2_full)

  d_sig <- d %>% filter(sig)
  d_ns  <- d %>% filter(!sig)
  if (nrow(d_ns) > 60000) d_ns <- d_ns[sample.int(nrow(d_ns), 60000), ]

  swarm_layer <- if (have_beeswarm) {
    list(
      ggbeeswarm::geom_quasirandom(data = d_ns,
                                    aes(x = logFC, y = L2_full, colour = index_study),
                                    groupOnX = FALSE, dodge.width = 0.85,
                                    size = 0.15, alpha = 0.15, bandwidth = 0.5,
                                    show.legend = FALSE),
      ggbeeswarm::geom_quasirandom(data = d_sig,
                                    aes(x = logFC, y = L2_full,
                                        colour = index_study, group = index_study),
                                    groupOnX = FALSE, dodge.width = 0.85,
                                    size = 0.4, alpha = 0.85, bandwidth = 0.5)
    )
  } else {
    list(
      geom_jitter(data = d_ns, aes(x = logFC, y = L2_full, colour = index_study),
                  position = position_jitterdodge(0.18, dodge.width = 0.85),
                  size = 0.15, alpha = 0.15, show.legend = FALSE),
      geom_jitter(data = d_sig, aes(x = logFC, y = L2_full, colour = index_study),
                  position = position_jitterdodge(0.18, dodge.width = 0.85),
                  size = 0.4, alpha = 0.85)
    )
  }

  # Per-(L2 × study × contrast × compartment) mean ± SD over all nhoods (sig + ns).
  stat_df <- d %>%
    group_by(L2_full, L2_compartment, contrast, index_study) %>%
    summarise(mean_lfc = mean(logFC, na.rm = TRUE),
              sd_lfc   = stats::sd(logFC, na.rm = TRUE),
              n        = dplyr::n(),
              .groups  = "drop") %>%
    filter(is.finite(mean_lfc) & is.finite(sd_lfc) & n >= 30) %>%
    mutate(xmin = mean_lfc - sd_lfc, xmax = mean_lfc + sd_lfc)

  p <- ggplot() +
    swarm_layer +
    geom_errorbarh(data = stat_df,
                   aes(y = L2_full, xmin = xmin, xmax = xmax,
                       group = index_study),
                   position = position_dodge(width = 0.85),
                   height = 0.5, linewidth = 0.35,
                   colour = "grey25", alpha = 0.9, inherit.aes = FALSE) +
    geom_point(data = stat_df,
               aes(x = mean_lfc, y = L2_full, fill = index_study,
                   group = index_study),
               position = position_dodge(width = 0.85),
               shape = 23, size = 1.6, colour = "grey15",
               stroke = 0.3, inherit.aes = FALSE) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "grey40", linewidth = 0.3) +
    facet_grid(L2_compartment ~ contrast, scales = "free_y", space = "free_y") +
    coord_cartesian(xlim = c(-4, 4)) +
    scale_fill_discrete(guide = "none") +
    theme_minimal(base_size = 9) +
    theme(axis.text.y = element_text(size = 6),
          axis.text.x = element_text(size = 8),
          panel.grid.major.y = element_blank(),
          strip.background = element_rect(fill = "grey95", colour = NA),
          strip.text = element_text(face = "bold", size = 9),
          legend.position = "right") +
    labs(title = main_title, subtitle = paste0(sub_title,
         " Diamond = mean logFC, error bar = ±1 SD (per study, n≥30 nhoods)."),
         x = "logFC", y = NULL, colour = "study")

  n_l2 <- nlevels(d$L2_full)
  n_c  <- length(available)
  fig_w <- 6 + 4 * n_c
  fig_h <- max(8, min(22, 0.42 * n_l2 + 4))
  ggsave(fig_path, p, width = fig_w, height = fig_h, limitsize = FALSE)
  cat(sprintf("  wrote: %s (%d L2s, %d contrasts, study-dodged, %.1fx%.1f in)\n",
              fig_path, n_l2, n_c, fig_w, fig_h))
}

build_per_study_beeswarm(
  contrast_set = c("parity_in_AR", "parity_x_HR_sporadic",
                    "parity_x_HR_BRCA1", "parity_x_HR_BRCA2"),
  fig_path = file.path(fig_dir, "da_beeswarm_per_study.pdf"),
  main_title = "Parity DA per study (transparency)",
  sub_title  = "AR baseline + 3 modifier tests; dodged + summarised by source study."
)

# Per-contrast volcanos -------------------------------------------------------
for (cname in names(all_da)) {
  d <- all_da[[cname]]
  fdr_t <- d$fdr_t[1]
  d$dir <- with(d, ifelse(sig & logFC > 0, "up",
                    ifelse(sig & logFC < 0, "dn", "ns")))
  pv <- ggplot(d, aes(x = logFC, y = neglog_fdr, colour = dir)) +
    geom_point(size = 0.4, alpha = 0.5) +
    scale_colour_manual(values = c(ns = "grey80", up = "#B2182B", dn = "#2166AC")) +
    geom_hline(yintercept = -log10(fdr_t), linetype = "dashed", colour = "grey30") +
    theme_minimal(base_size = 9) +
    labs(title = sprintf("DA volcano — %s", cname),
         subtitle = sprintf("nhoods=%d, sig=%d (FDR<%.2f)",
                            nrow(d), sum(d$sig), fdr_t),
         x = "logFC", y = "-log10(SpatialFDR)", colour = NULL)
  ggsave(file.path(fig_dir, sprintf("da_volcano_%s.pdf", cname)),
         pv, width = 6, height = 5)
}

# Stage E heatmap (if present) -----------------------------------------------
e_path <- paths$outputs$stageE
if (file.exists(e_path)) {
  e <- read_csv(e_path, show_col_types = FALSE)
  ns_cols <- grep("__n_sig$", colnames(e), value = TRUE)
  if (length(ns_cols) > 0) {
    long <- e %>%
      select(L2_joint, all_of(ns_cols)) %>%
      pivot_longer(-L2_joint, names_to = "contrast", values_to = "n_sig") %>%
      mutate(contrast = sub("__n_sig$", "", contrast))
    ord <- long %>% group_by(L2_joint) %>%
      summarise(max_n = max(n_sig, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(max_n)) %>% head(40) %>% pull(L2_joint)
    long_top <- long %>% filter(L2_joint %in% ord) %>%
      mutate(L2_joint = factor(L2_joint, levels = rev(ord)))
    p <- ggplot(long_top, aes(x = contrast, y = L2_joint, fill = n_sig)) +
      geom_tile(colour = "white") +
      geom_text(aes(label = ifelse(n_sig > 0, n_sig, "")),
                size = 2.6, colour = "white") +
      scale_fill_gradient(low = "grey90", high = "firebrick",
                           na.value = "grey95", trans = "log1p") +
      theme_minimal(base_size = 9) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
      labs(title = "Per-L2 sig nhood count by contrast (top 40 L2s)",
           x = NULL, y = NULL, fill = "n_sig")
    ggsave(file.path(fig_dir, "per_l2_summary_heatmap.pdf"),
           p, width = 8, height = 12)
    cat("  per_l2_summary_heatmap.pdf rendered\n")
  }
}

cat("=== plotting done ===\n")
