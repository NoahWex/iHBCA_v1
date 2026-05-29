#!/usr/bin/env Rscript
# stratification_two_tier.R
# Two-tier stratification: chooses L2 vs NhoodGroup resolution per parent_L2
# based on whether L2-level signed signal collapses below the nhoodgroup-level
# magnitude (i.e., L2 is heterogeneous and the action is at sub-L2).
#
# Filter: viable nhoodgroups only (n_nhoods_in_group >= 10, F.3 limma-viable).
#
# Per (parent_L2 × modifier contrast):
#   l2_signed_lfc   = mean(nhood logFC) over nhoods in viable nhoodgroups (w/in L2)
#   nhood_mean_abs  = mean(|group_med_lfc|) over viable nhoodgroups in L2
#   heterogeneous   = |l2_signed_lfc| < nhood_mean_abs
#
# Resolution choice per L2:
#   coherent       -> emit one row at L2 resolution
#   heterogeneous  -> emit one row per viable NhoodGroup
#
# Score: sign(beta_AR) * lfc  (same direction as AR baseline -> AMPLIFY-like).
#
# Outputs:
#   stratification_two_tier.csv
#   reports/figures/stratification_two_tier.pdf  (per-contrast bar chart)

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(ggplot2); library(stringr)
})
`%||%` <- function(a, b) if (is.null(a)) b else a

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "stratification_two_tier.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths
fig_dir <- file.path(paths$outputs$reports, "figures")
ensure_dir(fig_dir)
out_dir <- paths$inquiry_root

AR <- "parity_in_AR"
MODIFIERS <- c("parity_x_HR_BRCA1", "parity_x_HR_sporadic", "parity_x_HR_BRCA2")

VIABLE_MIN_NHOODS <- 10        # F.3 limma-viable threshold
MIN_N_SIG_AR      <- 30        # below this, AR sign is unreliable -> drop
PCT_AR_NEAR_ZERO  <- 0.02      # |beta_AR| below this -> sign unreliable

# ------------------------------------------------------------------------
# 1. AR baseline per L2 (from Stage E)
# ------------------------------------------------------------------------
e_wide <- read_csv(paths$outputs$stageE, show_col_types = FALSE)

ar_tbl <- e_wide %>%
  transmute(L2_joint, compartment, label,
            beta_AR = .data[[paste0(AR, "__med_lfc")]],
            n_sig_AR = .data[[paste0(AR, "__n_sig")]],
            n_total_AR = .data[[paste0(AR, "__n")]])

# AR sign — usable only if enough sig AR nhoods AND |beta_AR| not near zero
ar_tbl <- ar_tbl %>%
  mutate(ar_sign = case_when(
           is.na(beta_AR) | is.na(n_sig_AR) ~ NA_real_,
           n_sig_AR < MIN_N_SIG_AR          ~ NA_real_,
           abs(beta_AR) < PCT_AR_NEAR_ZERO  ~ NA_real_,
           beta_AR > 0                       ~ 1,
           beta_AR < 0                       ~ -1,
           TRUE                              ~ NA_real_))

cat(sprintf("AR baseline: %d L2s with usable sign (of %d total)\n",
            sum(!is.na(ar_tbl$ar_sign)), nrow(ar_tbl)))

# ------------------------------------------------------------------------
# 2. For each modifier contrast: load F.1 nhood-level + group-level
# ------------------------------------------------------------------------
process_contrast <- function(cname) {
  f1_dir <- file.path(paths$outputs$stageF1, cname)
  nhood_csv  <- file.path(f1_dir, "nhood_groups.csv")
  group_csv  <- file.path(f1_dir, "nhood_groups_summary.csv")
  if (!file.exists(nhood_csv) || !file.exists(group_csv)) {
    cat(sprintf("  [%s] F.1 outputs missing, skipping\n", cname)); return(NULL)
  }

  nh <- read_csv(nhood_csv, show_col_types = FALSE) %>%
    filter(!is.na(NhoodGroup_renamed))
  gr <- read_csv(group_csv, show_col_types = FALSE)

  # viable groups
  gr_viable <- gr %>% filter(n_nhoods_in_group >= VIABLE_MIN_NHOODS)
  viable_set <- gr_viable$NhoodGroup_renamed
  nh_viable <- nh %>% filter(NhoodGroup_renamed %in% viable_set)

  cat(sprintf("  [%s] viable groups: %d / %d (n_nhoods >= %d), nhoods covered: %d\n",
              cname, nrow(gr_viable), nrow(gr), VIABLE_MIN_NHOODS, nrow(nh_viable)))

  # 2a. L2-level signed LFC: mean over nhoods in viable groups within parent_L2
  l2_signed <- nh_viable %>%
    group_by(parent_L2_joint, parent_compartment, parent_label) %>%
    summarise(l2_signed_lfc = mean(logFC, na.rm = TRUE),
              n_nhoods_used = n(),
              .groups = "drop") %>%
    rename(L2_joint = parent_L2_joint,
           compartment = parent_compartment,
           label = parent_label)

  # 2b. Nhoodgroup aggregate per L2: mean(|group_med_lfc|) over viable groups
  l2_agg <- gr_viable %>%
    group_by(parent_L2_joint) %>%
    summarise(nhood_mean_abs = mean(abs(group_med_lfc), na.rm = TRUE),
              n_groups_viable = n(),
              .groups = "drop") %>%
    rename(L2_joint = parent_L2_joint)

  l2 <- l2_signed %>% inner_join(l2_agg, by = "L2_joint") %>%
    mutate(contrast = cname,
           heterogeneous = abs(l2_signed_lfc) < nhood_mean_abs)

  # 2c. NhoodGroup-level rows (used when heterogeneous)
  ng <- gr_viable %>%
    transmute(L2_joint = parent_L2_joint,
              compartment = parent_compartment,
              label = parent_label,
              NhoodGroup_renamed,
              group_med_lfc, n_nhoods_in_group, n_sig, group_pct_up,
              contrast = cname)

  list(l2 = l2, ng = ng)
}

ALL_CONTRASTS <- c(AR, MODIFIERS)
contrast_data <- lapply(ALL_CONTRASTS, process_contrast)
names(contrast_data) <- ALL_CONTRASTS

# ------------------------------------------------------------------------
# 3. Resolution choice + signed score
# AR contrast: signed_score = l2_signed_lfc directly (AR is its own reference).
# Modifier contrasts: signed_score = sign(β_AR) × lfc; rows without usable
#   AR sign are dropped from modifier output.
# ------------------------------------------------------------------------
# Classification helper: compare modifier lfc to AR baseline lfc per L2.
# Returns one of: AMPLIFY, ATTENUATE, REVERSE, ATTENUATE_REVERSE, no_baseline.
# AR rows always classified NA (AR is the baseline, not a modifier).
classify_vs_ar <- function(lfc, ar_lfc, n_sig_ar) {
  is_baseline_usable <- !is.na(ar_lfc) & !is.na(n_sig_ar) &
                          n_sig_ar >= MIN_N_SIG_AR &
                          abs(ar_lfc) >= PCT_AR_NEAR_ZERO
  out <- character(length(lfc))
  out[!is_baseline_usable] <- "no_baseline"
  has_base <- which(is_baseline_usable)
  if (length(has_base) > 0) {
    same_sign <- sign(lfc[has_base]) == sign(ar_lfc[has_base])
    larger    <- abs(lfc[has_base]) >= abs(ar_lfc[has_base])
    cls <- ifelse(same_sign & larger,  "AMPLIFY",
            ifelse(same_sign & !larger, "ATTENUATE",
              ifelse(!same_sign & larger, "REVERSE", "ATTENUATE_REVERSE")))
    out[has_base] <- cls
  }
  out
}

two_tier_rows <- list()
for (cname in ALL_CONTRASTS) {
  cd <- contrast_data[[cname]]
  if (is.null(cd)) next

  # Always join AR-baseline lookup; never filter out by usability (it's a
  # contextual annotation, not a ranking gate).
  l2 <- cd$l2 %>% left_join(
    ar_tbl %>% select(L2_joint, beta_AR, n_sig_AR), by = "L2_joint")
  ng <- cd$ng %>% left_join(
    ar_tbl %>% select(L2_joint, beta_AR, n_sig_AR), by = "L2_joint")

  # Coherent L2s: one row at L2 resolution
  coherent <- l2 %>% filter(!heterogeneous) %>%
    transmute(contrast,
              level = "L2",
              L2_joint, compartment, label,
              NhoodGroup_renamed = NA_character_,
              beta_AR, n_sig_AR,
              lfc = l2_signed_lfc,
              n_units = n_nhoods_used,
              n_sig = NA_integer_,
              n_groups_viable, nhood_mean_abs,
              heterogeneous = FALSE)

  # Heterogeneous L2s: drop to nhoodgroup rows
  het_l2s <- l2 %>% filter(heterogeneous) %>% select(L2_joint, n_groups_viable,
                                                       nhood_mean_abs)
  hetero <- ng %>% inner_join(het_l2s, by = "L2_joint") %>%
    transmute(contrast,
              level = "NhoodGroup",
              L2_joint, compartment, label,
              NhoodGroup_renamed,
              beta_AR, n_sig_AR,
              lfc = group_med_lfc,
              n_units = n_nhoods_in_group,
              n_sig,
              n_groups_viable, nhood_mean_abs,
              heterogeneous = TRUE)

  rows <- bind_rows(coherent, hetero) %>%
    mutate(abs_lfc   = abs(lfc),
           direction = ifelse(lfc > 0, "up",
                         ifelse(lfc < 0, "dn", "ns")),
           modifier_class = if (cname == AR) NA_character_
                              else classify_vs_ar(lfc, beta_AR, n_sig_AR))

  two_tier_rows[[cname]] <- rows
}

two_tier <- bind_rows(two_tier_rows) %>%
  arrange(contrast, desc(abs_lfc))

cat(sprintf("\nTwo-tier rows: %d\n", nrow(two_tier)))
cat("  by contrast:\n")
two_tier %>% count(contrast, level) %>% print()

# Output table
out_csv <- file.path(out_dir, "stratification_two_tier.csv")
write_csv(two_tier, out_csv)
cat(sprintf("Wrote: %s\n", out_csv))

# ------------------------------------------------------------------------
# 4. Plot: per-contrast bar chart
# ------------------------------------------------------------------------
contrast_label <- c(
  parity_in_AR = "AR baseline",
  parity_x_HR_BRCA1 = "BRCA1 modifier",
  parity_x_HR_sporadic = "sporadic modifier",
  parity_x_HR_BRCA2 = "BRCA2 modifier"
)

# Top 40 per contrast by abs_lfc (primary ranking). modifier_class encoded
# as suffix in row label (compact tag).
class_short <- c(
  AMPLIFY            = "[AMP]",
  ATTENUATE          = "[att]",
  ATTENUATE_REVERSE  = "[~rev]",
  REVERSE            = "[REV]",
  no_baseline        = "[nob]"
)

top_per_contrast <- two_tier %>%
  group_by(contrast) %>%
  slice_max(abs_lfc, n = 40, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(base_id = ifelse(level == "L2", L2_joint,
                            paste0(L2_joint, " :: ", NhoodGroup_renamed)),
         class_tag = ifelse(contrast == AR | is.na(modifier_class), "",
                              paste0(" ", class_short[modifier_class])),
         row_id = paste0(base_id, class_tag),
         contrast_label = factor(contrast_label[contrast],
                                  levels = contrast_label[ALL_CONTRASTS]))

p <- ggplot(top_per_contrast,
             aes(x = lfc,
                 y = reorder(row_id, lfc),
                 fill = direction)) +
  geom_col() +
  geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.4) +
  facet_wrap(~ contrast_label, scales = "free_y") +
  scale_fill_manual(values = c(up = "#B2182B", dn = "#2166AC", ns = "grey70"),
                     name = "Direction") +
  theme_minimal(base_size = 8) +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 6),
        legend.position = "top") +
  labs(x = "logFC (signed; ranked by |logFC|)",
       y = NULL,
       title = "Two-tier stratification: ranked by |logFC|. modifier_class shown as suffix tag.",
       subtitle = sprintf("Filter: n_nhoods_in_group ≥ %d. Top 40 per contrast. Tags: AMP/att/~rev/REV/nob (vs AR baseline).",
                            VIABLE_MIN_NHOODS))

out_pdf <- file.path(fig_dir, "stratification_two_tier.pdf")
ggsave(out_pdf, p, width = 14, height = 12)
cat(sprintf("Wrote: %s\n", out_pdf))

cat("\n=== two-tier stratification done ===\n")
