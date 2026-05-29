#!/usr/bin/env Rscript
# 38_stageG_per_finding_context.R
# Stage G — per-finding context (descriptive, NOT disqualifying).
#
# For each row of the per-L2 / nhoodgroup table, attach:
#   1. per-study contribution (n libraries, n donors, n_sig per study)
#   2. within-Kumar dissoc-adjustment context (CORRECTED FRAMING):
#        - is dissoc independent of parity within the dissoc-stratified study?
#        - if collinear: does coef magnitude survive adjustment?
#        - separately: stress-marker signature from Stage F.3
#      Verdict matrix (collinearity, stress) → verdict label
#   3. cohort-imbalance descriptors (per-study share of HR/proph/contra etc.)
#
# Caveats are CONTEXTUAL columns, NOT filters. Reader judges each row.
#
# Outputs:
#   stageG/per_finding_context.csv (one row per L2_joint × contrast,
#                                   plus one row per nhoodgroup × contrast)
#   stageG/dissoc_collinearity_summary.csv (per L2 × study collinearity)

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(yaml); library(broom)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "38_stageG_per_finding_context.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths
inquiry <- cfg$inquiry
out_csv <- paths$outputs$stageG
ensure_dir(dirname(out_csv))

cat(sprintf("=== Stage G: per-finding context for '%s' ===\n", inquiry$inquiry))

# -- 0. Composition-weighted collinearity (per L2 × flagged pair) ------------
# For each cohort:
#   - read cohort_<name>_l2_composition.csv (per L2 × study cell counts + share)
#   - read collinearity_matrix.csv (per cohort × within_study × pair × assoc)
#   - for each flagged (term1, term2) pair:
#       weighted_assoc[L2] = Σ_study (study_share_of_L2) × within-study assoc
#                            (using 0 when study not present in cohort)
# Flagged = pairs where any per-study verdict is COLLINEAR or ALIASED.
collin_path <- file.path(paths$outputs$stageB, "collinearity_matrix.csv")
weighted_collin <- list()
if (file.exists(collin_path)) {
  collin <- read_csv(collin_path, show_col_types = FALSE)
  flagged_pairs <- collin %>%
    filter(!is.na(within_study), !is.na(association),
           grepl("COLLINEAR|ALIASED", verdict)) %>%
    distinct(cohort, term1, term2)

  for (cn in names(inquiry$cohorts)) {
    comp_csv <- file.path(paths$outputs$stageA,
                           sprintf("cohort_%s_l2_composition.csv", cn))
    if (!file.exists(comp_csv)) next
    comp <- read_csv(comp_csv, show_col_types = FALSE)
    pairs_this <- flagged_pairs %>% filter(cohort == cn)
    for (k in seq_len(nrow(pairs_this))) {
      pr <- pairs_this[k, ]
      ws <- collin %>%
        filter(cohort == cn, term1 == pr$term1, term2 == pr$term2,
               !is.na(within_study), !is.na(association)) %>%
        select(study = within_study, study_assoc = association)
      if (nrow(ws) == 0) next
      # join study-level assoc onto each L2's per-study composition
      l2w <- comp %>%
        left_join(ws, by = "study") %>%
        mutate(study_assoc = ifelse(is.na(study_assoc), 0, study_assoc)) %>%
        group_by(L2_joint) %>%
        summarise(
          weighted_assoc = sum((n_cells / total_cells) * study_assoc),
          dominant_study = study[which.max(n_cells)],
          dominant_study_pct = max(pct_of_l2_cells),
          .groups = "drop"
        ) %>%
        mutate(cohort = cn,
               term_pair = sprintf("%s__%s", pr$term1, pr$term2))
      weighted_collin[[length(weighted_collin) + 1]] <- l2w
    }
  }
}
weighted_df <- bind_rows(weighted_collin)
if (nrow(weighted_df) > 0) {
  write_csv(weighted_df,
            file.path(dirname(out_csv), "stageG_l2_weighted_collinearity.csv"))
  cat(sprintf("Wrote: stageG_l2_weighted_collinearity.csv (%d rows)\n",
              nrow(weighted_df)))
}

# -- 1. Per-study contribution per L2 -----------------------------------------
# (already produced by Stage E as stageE_per_l2_per_study_long.csv)
study_long_path <- file.path(dirname(paths$outputs$stageE),
                             "stageE_per_l2_per_study_long.csv")
if (!file.exists(study_long_path)) {
  stop("Missing Stage E per-study long table. Run Stage E first.")
}
study_long <- read_csv(study_long_path, show_col_types = FALSE)

# Wide pivot: per L2_joint × contrast, n_sig & direction split per study
study_wide <- study_long %>%
  mutate(study = as.character(study)) %>%
  pivot_wider(
    id_cols     = c(L2_joint, contrast),
    names_from  = study,
    values_from = c(n, n_sig, n_sig_up, n_sig_dn),
    names_glue  = "study_{study}__{.value}"
  )

# Concordance: fraction of contributing studies whose direction matches the
# pooled direction.
direction_concordance <- study_long %>%
  group_by(L2_joint, contrast) %>%
  mutate(pooled_dir = ifelse(sum(n_sig_up) >= sum(n_sig_dn), "up", "dn")) %>%
  group_by(L2_joint, contrast, pooled_dir) %>%
  summarise(
    n_studies_contrib = sum(n_sig > 0),
    n_studies_match   = sum((pooled_dir == "up" & n_sig_up >= n_sig_dn) |
                            (pooled_dir == "dn" & n_sig_dn >  n_sig_up)),
    direction_concordance = ifelse(n_studies_contrib > 0,
                                   n_studies_match / n_studies_contrib, NA),
    .groups = "drop"
  ) %>%
  select(-pooled_dir)

# -- 2. Within-Kumar dissoc context (CORRECTED FRAMING) ----------------------
dissoc_study <- inquiry$context$dissoc_study %||% NA_character_
dissoc_table_path <- paths$inputs$dissoc_table
dissoc_ctx <- NULL
collin_summary <- NULL
if (!is.na(dissoc_study) && !is.null(dissoc_table_path) &&
    file.exists(dissoc_table_path)) {

  cat(sprintf("Dissoc context (study=%s, table=%s)\n",
              dissoc_study, dissoc_table_path))
  dissoc_tab <- read_csv(dissoc_table_path, show_col_types = FALSE)
  donor_col_dt <- intersect(c("donor_id", "ihbca_donor_id"),
                            colnames(dissoc_tab))[1]

  # Load harmonized donor metadata to get parity per donor in dissoc study
  hm <- read_csv(paths$inputs$donor_metadata, show_col_types = FALSE)
  donor_col_hm <- intersect(c("ihbca_donor_id", "donor_id"), colnames(hm))[1]
  hm_sub <- hm %>%
    filter(study == dissoc_study, !is.na(parity_binary)) %>%
    select(all_of(donor_col_hm), parity_binary)

  joined <- hm_sub %>%
    left_join(dissoc_tab %>% rename(.dt_donor = !!donor_col_dt),
              by = setNames(".dt_donor", donor_col_hm)) %>%
    filter(!is.na(dissoc_min), parity_binary %in% c("parous", "nulliparous"))

  # ---- KEY collinearity test: parity vs dissoc within Kumar ----
  # Wilcoxon (continuous dissoc by binary parity) + |Spearman| of parity vs dissoc
  if (nrow(joined) >= 8) {
    w <- suppressWarnings(
      wilcox.test(dissoc_min ~ parity_binary, data = joined))
    s <- suppressWarnings(cor(as.integer(joined$parity_binary == "parous"),
                              joined$dissoc_min, method = "spearman"))
    median_par  <- median(joined$dissoc_min[joined$parity_binary == "parous"])
    median_null <- median(joined$dissoc_min[joined$parity_binary == "nulliparous"])
    collin_summary <- tibble(
      study              = dissoc_study,
      n_donors           = nrow(joined),
      n_parous           = sum(joined$parity_binary == "parous"),
      n_nulliparous      = sum(joined$parity_binary == "nulliparous"),
      median_dissoc_parous     = median_par,
      median_dissoc_nulliparous = median_null,
      wilcox_p           = signif(w$p.value, 3),
      spearman_r         = round(s, 3),
      collinear          = !is.na(w$p.value) && w$p.value < 0.05
    )
    cat(sprintf("Kumar parity-dissoc collinearity: wilcox p=%g  spearman=%.3f  collinear=%s\n",
                collin_summary$wilcox_p, collin_summary$spearman_r,
                collin_summary$collinear))
  }

  # ---- Per-L2 dissoc-adjustment within dissoc study ----
  # Compute per-L2 abundance from milo? No — too expensive at this stage.
  # Instead, derive from per-study DA results: for each contrast × L2,
  # report per-study n_sig (already in study_wide). Dissoc adjustment lives
  # at script 19 / 29 level; we attach the COEF-RETENTION per L2 if a precomputed
  # systematic_l2_dissoc CSV exists in the inquiry's config dir or context.
  presaved_dissoc <- inquiry$context$systematic_dissoc_csv %||% NA_character_
  if (!is.na(presaved_dissoc) && file.exists(presaved_dissoc)) {
    sd <- read_csv(presaved_dissoc, show_col_types = FALSE)
    dissoc_ctx <- sd %>%
      transmute(
        L2_joint               = L2_joint,
        kumar_dissoc_n_donors  = n_donors,
        kumar_parity_naive_p   = parity_naive_p,
        kumar_parity_adj_p     = parity_adj_p,
        kumar_dissoc_p         = dissoc_p,
        kumar_coef_retention_pct = coef_retention_pct
      )
  }
}

# -- 3. Stress-marker signature per nhoodgroup (Stage F.3) -------------------
stress_rows <- list()
contrast_dirs <- list.dirs(paths$outputs$stageF3, recursive = FALSE)
for (cdir in contrast_dirs) {
  cname <- basename(cdir)
  csvs <- list.files(cdir, pattern = "_vs_parent\\.csv$", full.names = TRUE)
  for (f in csvs) {
    m <- read_csv(f, show_col_types = FALSE)
    if (!"FDR" %in% colnames(m) || !"is_dissoc_stress" %in% colnames(m)) next
    grp <- sub("_vs_parent\\.csv$", "", basename(f))
    sig <- m %>% filter(!is.na(FDR), FDR < 0.05)
    n_sig <- nrow(sig)
    n_sig_up <- sum(sig$logFC > 0, na.rm = TRUE)
    n_sig_stress <- sum(sig$is_dissoc_stress, na.rm = TRUE)
    n_sig_stress_up <- sum(sig$is_dissoc_stress & sig$logFC > 0, na.rm = TRUE)
    stress_rows[[length(stress_rows) + 1]] <- tibble(
      contrast      = cname,
      group_renamed = grp,
      markers_total_sig    = n_sig,
      markers_total_sig_up = n_sig_up,
      markers_stress_sig   = n_sig_stress,
      markers_stress_sig_up = n_sig_stress_up,
      stress_fraction_up   = ifelse(n_sig_up > 0, n_sig_stress_up / n_sig_up, NA)
    )
  }
}
stress_df <- bind_rows(stress_rows)

# -- 4. Verdict synthesis -----------------------------------------------------
# verdict per (L2 × contrast):
#   - if dissoc study isn't in the cohort: NA
#   - if Kumar parity-dissoc independent: PARITY_DISSOC_INDEPENDENT
#   - if collinear and coef retention high (>=50%) and stress-fraction <= 0.10: PARITY_RECOVERY_RATE_LIKELY
#   - if collinear and stress-fraction > 0.10: PARITY_STRESS_ARTIFACT_RISK
#   - if no Kumar contribution to this L2: KUMAR_NOT_CONTRIBUTING
verdict_for_row <- function(retention, stress_frac, collinear,
                             kumar_n_sig) {
  if (is.na(collinear)) return(NA_character_)
  if (!isTRUE(collinear)) return("PARITY_DISSOC_INDEPENDENT")
  if (is.na(kumar_n_sig) || kumar_n_sig == 0) return("KUMAR_NOT_CONTRIBUTING")
  if (!is.na(stress_frac) && stress_frac > 0.10) return("PARITY_STRESS_ARTIFACT_RISK")
  if (!is.na(retention) && abs(retention) >= 50)  return("PARITY_RECOVERY_RATE_LIKELY")
  "PARITY_DISSOC_AMBIGUOUS"
}

# -- 5. Assemble Stage G table ------------------------------------------------
out <- study_wide %>%
  left_join(direction_concordance, by = c("L2_joint", "contrast"))

# attach dissoc context if available
if (!is.null(dissoc_ctx)) {
  out <- out %>% left_join(dissoc_ctx, by = "L2_joint")
} else {
  out$kumar_dissoc_n_donors <- NA_integer_
  out$kumar_parity_naive_p  <- NA_real_
  out$kumar_parity_adj_p    <- NA_real_
  out$kumar_dissoc_p        <- NA_real_
  out$kumar_coef_retention_pct <- NA_real_
}

# attach collinearity flag (study-level constant)
collinear_flag <- if (!is.null(collin_summary)) collin_summary$collinear else NA
out$dissoc_collinear_with_parity <- collinear_flag

# attach Kumar n_sig per L2 (already in study_wide as
# study_kumar__n_sig if Kumar is named "kumar")
kumar_col <- grep("^study_kumar__n_sig$", colnames(out), value = TRUE)[1]
if (length(kumar_col) > 0) {
  out$kumar_n_sig <- out[[kumar_col]]
} else {
  out$kumar_n_sig <- NA_integer_
}

out$verdict_dissoc <- mapply(verdict_for_row,
                              retention   = out$kumar_coef_retention_pct,
                              stress_frac = NA, # filled in nhoodgroup-level table
                              collinear   = collinear_flag,
                              kumar_n_sig = out$kumar_n_sig)

write_csv(out, out_csv)

# Per-nhoodgroup table with stress info
if (nrow(stress_df) > 0) {
  ng_csv <- file.path(dirname(out_csv), "stageG_per_nhoodgroup_context.csv")
  ng_out <- stress_df %>%
    mutate(verdict_dissoc = mapply(verdict_for_row,
                                    retention = NA,
                                    stress_frac = stress_fraction_up,
                                    collinear = collinear_flag,
                                    kumar_n_sig = NA))
  write_csv(ng_out, ng_csv)
  cat(sprintf("Wrote: %s  (%d nhoodgroups)\n", ng_csv, nrow(ng_out)))
}

if (!is.null(collin_summary)) {
  write_csv(collin_summary, file.path(dirname(out_csv),
                                       "stageG_dissoc_collinearity_summary.csv"))
}

cat(sprintf("Wrote: %s  (%d rows)\n", out_csv, nrow(out)))
cat("=== Stage G done ===\n")
