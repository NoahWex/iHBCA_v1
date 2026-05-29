#!/usr/bin/env Rscript
# fig2_mockup.R
# Renders a fully-composed Fig 2 mockup for the parity-by-risk DA inquiry.
# 6 panels demonstrating the multi-study atlas DA framework:
#   A. Cohort substrate (donor metadata grid)
#   B. Per-nhood DA hero beeswarm (3 strata)
#   C. Stratification ranking heatmap (top-N populations × 4 contrasts)
#   D. Illustrative L2 depth: epi::LASP-major (per-nhood + top pathways)
#   E. Pathway summary dot heatmap (top L2s × top pathways × 3 contrasts)
#   F. Cross-stratum overlap of significant populations
# Output: reports/figures/fig2_mockup.pdf

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(tidyr)
  library(patchwork); library(stringr)
})
have_beeswarm <- requireNamespace("ggbeeswarm", quietly = TRUE)

`%||%` <- function(a, b) if (is.null(a)) b else a

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "fig2_mockup.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths
fig_dir <- file.path(paths$outputs$reports, "figures")
ensure_dir(fig_dir)

theme_fig2 <- function() {
  theme_minimal(base_size = 7) +
    theme(panel.grid = element_blank(),
          axis.line = element_line(linewidth = 0.25),
          axis.ticks = element_line(linewidth = 0.25),
          plot.title = element_text(face = "bold", size = 8),
          legend.key.size = unit(0.3, "cm"),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6, face = "bold"),
          strip.text = element_text(face = "bold", size = 7))
}
theme_set(theme_fig2())

CONTRAST_LABEL <- c(
  parity_in_AR = "AR baseline",
  parity_x_HR_BRCA1 = "BRCA1 modifier",
  parity_x_HR_sporadic = "sporadic modifier",
  parity_x_HR_BRCA2 = "BRCA2 modifier"
)
CONTRAST_ORDER <- c("parity_in_AR", "parity_x_HR_BRCA1",
                     "parity_x_HR_sporadic", "parity_x_HR_BRCA2")

# ============================================================
# Panel A — Cohort substrate
# ============================================================
cat("Panel A: cohort substrate ...\n")
cohort_files <- list(
  AR_baseline   = "cohort_AR_baseline_design.csv",
  HR_sporadic   = "cohort_HR_sporadic_design.csv",
  HR_BRCA1      = "cohort_HR_BRCA1_design.csv",
  HR_BRCA2      = "cohort_HR_BRCA2_design.csv"
)
cohort_dfs <- list()
for (cn in names(cohort_files)) {
  fp <- file.path(paths$outputs$stageA, cohort_files[[cn]])
  if (!file.exists(fp)) next
  d <- read_csv(fp, show_col_types = FALSE) %>% mutate(cohort = cn)
  cohort_dfs[[cn]] <- d
}
ck <- bind_rows(cohort_dfs)

# Donor identifier
donor_col <- intersect(c("ihbca_donor_id", "donor_id", "patientID"),
                        colnames(ck))[1]
ck$donor_id <- ck[[donor_col]]

ck$cohort <- factor(ck$cohort, levels = c("AR_baseline", "HR_sporadic",
                                            "HR_BRCA1", "HR_BRCA2"))
ck$parity_binary <- factor(ck$parity_binary,
                            levels = c("nulliparous", "parous"))

# Order donors: by cohort, then parity, then age
ck <- ck %>% arrange(cohort, parity_binary, age_continuous)
ck$donor_idx <- seq_len(nrow(ck))

# Long-form metadata strip
strip_vars <- c("cohort", "parity_binary", "study", "menopausal_status_binary",
                  "facs_status", "brca_genotype", "cancer_history")
strip_vars <- intersect(strip_vars, colnames(ck))
strip_long <- ck %>%
  select(donor_idx, all_of(strip_vars)) %>%
  pivot_longer(-donor_idx, names_to = "var", values_to = "val") %>%
  mutate(var = factor(var, levels = strip_vars))

panel_A <- ggplot(strip_long, aes(x = donor_idx, y = var, fill = val)) +
  geom_tile(linewidth = 0) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(limits = rev(strip_vars), expand = c(0, 0)) +
  scale_fill_viridis_d(option = "turbo", na.value = "grey90", guide = "none") +
  labs(title = "A. Cohort substrate",
        subtitle = sprintf("%d donors across 5 within-stratum parity contrasts (BRCA12 sensitivity overlaps shown twice)",
                            n_distinct(ck$donor_idx)),
        x = "Donor (sorted: cohort, parity, age)", y = NULL)

# ============================================================
# Panel B — Hero beeswarm
# ============================================================
cat("Panel B: hero beeswarm ...\n")
B_contrasts <- c("parity_in_AR", "parity_x_HR_sporadic", "parity_x_HR_BRCA1")
da_rows <- list()
for (cn in B_contrasts) {
  fp <- file.path(paths$outputs$stageD, cn, "da_results.csv")
  if (!file.exists(fp)) next
  d <- read_csv(fp, show_col_types = FALSE) %>%
    mutate(contrast = cn,
           sig_dir = ifelse(SpatialFDR < 0.05 & logFC > 0, "up",
                       ifelse(SpatialFDR < 0.05 & logFC < 0, "dn", "ns")))
  if ("is_artifact" %in% colnames(d)) d <- d %>% filter(!is_artifact)
  da_rows[[cn]] <- d
}
da <- bind_rows(da_rows)
da$contrast <- factor(da$contrast, levels = B_contrasts,
                       labels = CONTRAST_LABEL[B_contrasts])

# Restrict L2s shown: those with any sig nhood in any contrast
l2_keep <- da %>% group_by(label, compartment) %>%
  summarise(any_sig = any(sig_dir != "ns"), .groups = "drop") %>%
  filter(any_sig) %>%
  pull(label)
da_plot <- da %>% filter(label %in% l2_keep)
da_plot$compartment <- factor(da_plot$compartment, levels = c("epi", "imm", "str"))
# Sample down ns dots for legibility
set.seed(1)
da_sig <- da_plot %>% filter(sig_dir != "ns")
da_ns  <- da_plot %>% filter(sig_dir == "ns")
if (nrow(da_ns) > 30000) da_ns <- da_ns[sample.int(nrow(da_ns), 30000), ]
da_plot2 <- bind_rows(da_sig, da_ns)

swarm <- if (have_beeswarm) {
  ggbeeswarm::geom_quasirandom(aes(x = logFC, y = label, colour = sig_dir),
                                groupOnX = FALSE, size = 0.15, alpha = 0.5,
                                bandwidth = 0.5)
} else {
  geom_jitter(aes(x = logFC, y = label, colour = sig_dir),
              size = 0.15, alpha = 0.5, width = 0, height = 0.25)
}

panel_B <- ggplot(da_plot2) + swarm +
  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.3) +
  facet_grid(compartment ~ contrast, scales = "free_y", space = "free_y") +
  scale_colour_manual(values = c(ns = "grey75", up = "#B2182B", dn = "#2166AC"),
                       guide = guide_legend(override.aes = list(size = 1.5,
                                                                  alpha = 1))) +
  coord_cartesian(xlim = c(-4, 4)) +
  labs(title = "B. Per-nhood DA across strata",
        x = "logFC", y = NULL, colour = "FDR<0.05")

# ============================================================
# Panel C — Stratification ranking heatmap
# ============================================================
cat("Panel C: stratification heatmap ...\n")
strat <- read_csv(file.path(paths$inquiry_root, "stratification_two_tier.csv"),
                   show_col_types = FALSE)
# Restrict: top by abs_lfc within each contrast, union across contrasts
top_units <- strat %>%
  group_by(contrast) %>%
  slice_max(abs_lfc, n = 12, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(unit_id = ifelse(level == "L2", L2_joint,
                            paste0(L2_joint, ":", NhoodGroup_renamed)))
unit_order <- unique(top_units$unit_id)

# Re-pull all (contrast × unit) entries from the full table for these units
full_C <- strat %>%
  mutate(unit_id = ifelse(level == "L2", L2_joint,
                            paste0(L2_joint, ":", NhoodGroup_renamed))) %>%
  filter(unit_id %in% unit_order)

# Wide for ordering, then long for plot
full_C$contrast <- factor(full_C$contrast, levels = CONTRAST_ORDER,
                           labels = CONTRAST_LABEL[CONTRAST_ORDER])

# Use mean lfc across BRCA1+sporadic for ordering rows
order_score <- full_C %>%
  filter(contrast %in% c("BRCA1 modifier", "sporadic modifier")) %>%
  group_by(unit_id) %>%
  summarise(score = mean(abs_lfc, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(score))
full_C$unit_id <- factor(full_C$unit_id, levels = order_score$unit_id)

panel_C <- ggplot(full_C, aes(x = contrast, y = unit_id, fill = lfc)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.1f", lfc)), size = 1.7) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                        midpoint = 0, limits = c(-5, 5),
                        oob = scales::squish, na.value = "grey90",
                        name = "logFC") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev(order_score$unit_id)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 0, size = 6.5),
        axis.text.y = element_text(size = 6)) +
  labs(title = "C. Stratification ranking",
        subtitle = "Top units across non-AR contrasts (mean |logFC| sort)",
        x = NULL, y = NULL)

# ============================================================
# Panel D — Illustrative L2: epi::LASP-major
# ============================================================
cat("Panel D: LASP-major depth ...\n")
ILLUSTRATIVE_L2 <- "epi::LASP-major"

# 4-row beeswarm of nhood-level logFC across all 4 contrasts
da_all_rows <- list()
for (cn in CONTRAST_ORDER) {
  fp <- file.path(paths$outputs$stageD, cn, "da_results.csv")
  if (!file.exists(fp)) next
  d <- read_csv(fp, show_col_types = FALSE) %>%
    mutate(contrast = cn,
           sig_dir = ifelse(SpatialFDR < 0.05 & logFC > 0, "up",
                       ifelse(SpatialFDR < 0.05 & logFC < 0, "dn", "ns")))
  if ("is_artifact" %in% colnames(d)) d <- d %>% filter(!is_artifact)
  da_all_rows[[cn]] <- d
}
da_all <- bind_rows(da_all_rows)

l2_da <- da_all %>%
  mutate(L2_joint = paste0(compartment, "::", label)) %>%
  filter(L2_joint == ILLUSTRATIVE_L2)
l2_da$contrast <- factor(l2_da$contrast, levels = CONTRAST_ORDER,
                          labels = CONTRAST_LABEL[CONTRAST_ORDER])

swarm_d <- if (have_beeswarm) {
  ggbeeswarm::geom_quasirandom(aes(x = logFC, y = contrast, colour = sig_dir),
                                groupOnX = FALSE, size = 0.4, alpha = 0.6,
                                bandwidth = 0.5)
} else {
  geom_jitter(aes(x = logFC, y = contrast, colour = sig_dir),
              size = 0.4, alpha = 0.6, height = 0.2, width = 0)
}
panel_D_da <- ggplot(l2_da) + swarm_d +
  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.3) +
  scale_colour_manual(values = c(ns = "grey75", up = "#B2182B", dn = "#2166AC"),
                       guide = "none") +
  coord_cartesian(xlim = c(-5, 5)) +
  labs(title = sprintf("D. %s — per-nhood DA", ILLUSTRATIVE_L2),
        x = "logFC", y = NULL)

# Top pathways for this L2 from per_L2_nes
nes_path <- file.path(paths$inquiry_root, "outputs", "stageI3_gsea",
                       "per_L2_nes.csv")
if (file.exists(nes_path)) {
  nes <- read_csv(nes_path, show_col_types = FALSE)
  l2_nes <- nes %>%
    filter(parent_L2_joint == ILLUSTRATIVE_L2,
            !is.na(padj), padj < 0.01,
            collection == "hallmark" |
              (collection == "reactome" & abs(NES) >= 2.5)) %>%
    group_by(contrast, pathway) %>% slice_max(abs(NES), n = 1) %>%
    ungroup()
  # Pick top 12 pathways by max |NES| across contrasts shown
  pathway_keep <- l2_nes %>% group_by(pathway) %>%
    summarise(maxabs = max(abs(NES)), .groups = "drop") %>%
    arrange(desc(maxabs)) %>% head(12) %>% pull(pathway)
  l2_nes_plot <- l2_nes %>% filter(pathway %in% pathway_keep) %>%
    mutate(contrast = factor(contrast,
                               levels = c("parity_in_AR", "parity_x_HR_BRCA1",
                                          "parity_x_HR_sporadic"),
                               labels = c("AR", "BRCA1", "sporadic")),
           # short label
           pathway_short = str_replace_all(pathway,
                                            c("HALLMARK_" = "H_",
                                              "REACTOME_" = "R_",
                                              "_" = " ")) %>%
                             str_trunc(45))

  pathway_order <- l2_nes_plot %>% group_by(pathway, pathway_short) %>%
    summarise(maxabs = max(abs(NES)), .groups = "drop") %>%
    arrange(desc(maxabs))
  l2_nes_plot$pathway_short <- factor(l2_nes_plot$pathway_short,
                                       levels = pathway_order$pathway_short)

  panel_D_path <- ggplot(l2_nes_plot,
                          aes(x = contrast, y = pathway_short)) +
    geom_point(aes(size = -log10(padj), fill = NES),
                shape = 21, stroke = 0.2, colour = "black") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                          midpoint = 0, limits = c(-3.5, 3.5),
                          oob = scales::squish, name = "NES") +
    scale_size_continuous(range = c(1, 4), name = "-log10(padj)") +
    scale_y_discrete(limits = rev(levels(l2_nes_plot$pathway_short))) +
    labs(title = sprintf("Top pathways · %s", ILLUSTRATIVE_L2),
          x = NULL, y = NULL) +
    theme(axis.text.y = element_text(size = 5.5))

  panel_D <- panel_D_da / panel_D_path + plot_layout(heights = c(1, 2.5))
} else {
  panel_D <- panel_D_da
}

# ============================================================
# Panel E — Pathway summary dot heatmap
# ============================================================
cat("Panel E: pathway dot heatmap ...\n")
nes_path <- file.path(paths$inquiry_root, "outputs", "stageI3_gsea",
                       "per_L2_nes.csv")
if (file.exists(nes_path)) {
  nes <- read_csv(nes_path, show_col_types = FALSE)
  # Top-6 hallmark pathways per (contrast × L2) subset for the matrix
  E_data <- nes %>%
    filter(collection == "hallmark", padj < 0.05,
            contrast %in% c("parity_in_AR", "parity_x_HR_BRCA1",
                            "parity_x_HR_sporadic")) %>%
    mutate(contrast_lbl = recode(contrast,
                                   parity_in_AR = "AR",
                                   parity_x_HR_BRCA1 = "BRCA1",
                                   parity_x_HR_sporadic = "sporadic"))
  # Top L2s by recurrence
  l2_count <- E_data %>% group_by(parent_L2_joint) %>%
    summarise(n_sig_paths = n(), .groups = "drop") %>%
    arrange(desc(n_sig_paths)) %>% head(15)
  pw_count <- E_data %>%
    filter(parent_L2_joint %in% l2_count$parent_L2_joint) %>%
    group_by(pathway) %>%
    summarise(n_groups = n(), max_abs = max(abs(NES)), .groups = "drop") %>%
    arrange(desc(n_groups), desc(max_abs)) %>% head(20)
  E_plot <- E_data %>%
    filter(parent_L2_joint %in% l2_count$parent_L2_joint,
            pathway %in% pw_count$pathway) %>%
    mutate(pathway_short = str_replace(pathway, "HALLMARK_", "") %>%
                            str_replace_all("_", " ") %>% str_trunc(28),
           contrast_lbl = factor(contrast_lbl,
                                   levels = c("AR", "sporadic", "BRCA1")))
  E_plot$parent_L2_joint <- factor(E_plot$parent_L2_joint,
                                     levels = l2_count$parent_L2_joint)
  E_plot$pathway_short <- factor(E_plot$pathway_short,
                                   levels = unique(E_plot$pathway_short[
                                     order(match(E_plot$pathway,
                                                   pw_count$pathway))]))

  panel_E <- ggplot(E_plot,
                     aes(x = parent_L2_joint, y = pathway_short)) +
    geom_point(aes(size = -log10(padj), fill = NES),
                shape = 21, stroke = 0.2, colour = "black") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                          midpoint = 0, limits = c(-3, 3),
                          oob = scales::squish, name = "NES") +
    scale_size_continuous(range = c(0.5, 3), name = "-log10(padj)") +
    facet_wrap(~ contrast_lbl, ncol = 1) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5),
          axis.text.y = element_text(size = 5.5)) +
    labs(title = "E. Pathway summary",
          subtitle = "Top hallmark pathways × top L2s × 3 contrasts",
          x = NULL, y = NULL)
} else {
  panel_E <- ggplot() + theme_void() +
    labs(title = "E. Pathway summary (data missing)")
}

# ============================================================
# Panel F — Cross-stratum overlap of significant populations
# ============================================================
cat("Panel F: cross-stratum overlap ...\n")
sig_overlap <- strat %>%
  mutate(unit_id = ifelse(level == "L2", L2_joint,
                            paste0(L2_joint, ":", NhoodGroup_renamed)),
          is_strong = abs_lfc >= 1.5) %>%
  filter(is_strong) %>%
  group_by(unit_id) %>%
  summarise(strata = paste(sort(unique(contrast)), collapse = "+"),
             n_strata = n_distinct(contrast),
             max_abs = max(abs_lfc),
             .groups = "drop")

# For UpSet-style: count units in each combination
combo_counts <- sig_overlap %>%
  count(strata, n_strata) %>%
  arrange(desc(n_strata), desc(n))
combo_counts$strata_lbl <- str_replace_all(combo_counts$strata,
                                            c("parity_in_AR" = "AR",
                                              "parity_x_HR_BRCA1" = "BRCA1",
                                              "parity_x_HR_sporadic" = "spor",
                                              "parity_x_HR_BRCA2" = "BRCA2"))

panel_F <- ggplot(combo_counts %>% head(15),
                   aes(x = reorder(strata_lbl, n), y = n,
                       fill = factor(n_strata))) +
  geom_col() +
  geom_text(aes(label = n), hjust = -0.2, size = 2.2) +
  scale_fill_manual(values = c(`1` = "#FB6A4A", `2` = "#9ECAE1",
                                 `3` = "#3182BD", `4` = "#08306B"),
                     name = "# strata") +
  coord_flip(clip = "off") +
  labs(title = "F. Cross-stratum overlap",
        subtitle = "Populations with |logFC| ≥ 1.5 by which contrasts they appear in",
        x = NULL, y = "# units") +
  theme(plot.margin = margin(5, 30, 5, 5))

# ============================================================
# Compose
# ============================================================
cat("Composing fig2_mockup.pdf ...\n")

top_row <- panel_A + plot_layout(widths = 1)
mid_row <- (panel_B | panel_C) + plot_layout(widths = c(1.5, 1))
bot_row <- (panel_D | (panel_E / panel_F + plot_layout(heights = c(2, 1)))) +
            plot_layout(widths = c(1.4, 1))

fig2 <- top_row / mid_row / bot_row +
          plot_layout(heights = c(0.6, 1.6, 2.5)) +
          plot_annotation(title = "Fig 2 (mockup) — Demonstration of multi-study atlas DA framework: parity × risk-stratum",
                            theme = theme(plot.title = element_text(face = "bold",
                                                                       size = 10)))

out_pdf <- file.path(fig_dir, "fig2_mockup.pdf")
ggsave(out_pdf, fig2, width = 14, height = 18, limitsize = FALSE)
cat(sprintf("Wrote: %s\n", out_pdf))

cat("\n=== fig2 mockup done ===\n")
