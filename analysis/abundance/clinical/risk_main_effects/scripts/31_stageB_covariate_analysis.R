#!/usr/bin/env Rscript
# 31_stageB_covariate_analysis.R
# Stage B — covariate analysis (Phase B per da_framework.md).
# Defines what each cohort CAN and CANNOT answer before any DA test runs.
#
# Inputs:
#   stageA cohort_<name>_design.csv (per cohort)
#   inquiry$terms_of_interest
#
# Outputs:
#   stageB/cohort_landscape.csv
#       per cohort × covariate × level: n libraries, n donors
#   stageB/collinearity_matrix.csv
#       per cohort × pair-of-terms: Cramér's V (categorical) or |Spearman| (numeric)
#       + identifiability flag (V > 0.95 -> aliased)
#   stageB/per_study_testability.csv
#       per cohort × study × term: levels present, estimable per-study?, estimable pooled?
#   stageB/STAGEB_covariate_landscape.md
#       written summary with red/yellow/green flags per term per cohort
#
# Usage:
#   Rscript 31_stageB_covariate_analysis.R --inquiry-dir /path/to/inquiry

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(purrr)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "31_stageB_covariate_analysis.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths
inquiry <- cfg$inquiry
in_dir <- paths$outputs$stageA
out_dir <- paths$outputs$stageB
ensure_dir(out_dir)

cat(sprintf("=== Stage B: covariate analysis for inquiry '%s' ===\n", inquiry$inquiry))

# -- helpers -------------------------------------------------------------------

cramers_v <- function(x, y) {
  ok <- complete.cases(x, y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 2) return(NA_real_)
  if (length(unique(x)) < 2 || length(unique(y)) < 2) return(NA_real_)
  tab <- suppressWarnings(table(x, y))
  chi <- suppressWarnings(chisq.test(tab, correct = FALSE))$statistic
  n <- sum(tab)
  k <- min(nrow(tab), ncol(tab))
  unname(sqrt(chi / (n * (k - 1))))
}

# assoc(): unified pairwise association in [0, 1]-ish.
# - numeric × numeric:        |Spearman r|
# - categorical × numeric:    sqrt(eta² from Kruskal-Wallis), in [0,1]
# - categorical × categorical: Cramér's V, in [0,1]
# Also returns the test name + p-value so we can flag significance, not just
# magnitude. Returns a list(assoc, p, test).
assoc <- function(x, y) {
  ok <- complete.cases(x, y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 3) return(list(assoc = NA_real_, p = NA_real_, test = "n<3"))
  xn <- is.numeric(x); yn <- is.numeric(y)
  if (xn && yn) {
    ct <- suppressWarnings(cor.test(x, y, method = "spearman"))
    return(list(assoc = abs(unname(ct$estimate)), p = ct$p.value, test = "spearman"))
  }
  # one numeric + one categorical
  if (xn != yn) {
    num <- if (xn) x else y
    cat <- as.character(if (xn) y else x)
    if (length(unique(cat)) < 2) return(list(assoc = NA_real_, p = NA_real_, test = "1-level"))
    kw <- suppressWarnings(kruskal.test(num, cat))
    n <- length(num); k <- length(unique(cat))
    # eta² approx for KW: (H - k + 1) / (n - k); clamp to [0,1]
    eta2 <- (unname(kw$statistic) - k + 1) / (n - k)
    if (is.na(eta2) || eta2 < 0) eta2 <- 0
    return(list(assoc = sqrt(min(eta2, 1)), p = kw$p.value, test = "kruskal"))
  }
  # both categorical
  xs <- as.character(x); ys <- as.character(y)
  if (length(unique(xs)) < 2 || length(unique(ys)) < 2) {
    return(list(assoc = NA_real_, p = NA_real_, test = "1-level"))
  }
  v <- cramers_v(xs, ys)
  tab <- suppressWarnings(table(xs, ys))
  if (any(dim(tab) < 2) || sum(tab) < 2) {
    return(list(assoc = v, p = NA_real_, test = "chisq-deg"))
  }
  ch <- tryCatch(suppressWarnings(chisq.test(tab, correct = FALSE)),
                 error = function(e) NULL)
  list(assoc = v, p = if (is.null(ch)) NA_real_ else ch$p.value, test = "chisq")
}

# -- main loop over cohorts ----------------------------------------------------

landscape_rows  <- list()
collin_rows     <- list()
testability_rows <- list()
md_lines <- c(sprintf("# Stage B covariate landscape — inquiry: %s",
                       inquiry$inquiry), "")

for (cohort_name in names(inquiry$cohorts)) {
  spec <- inquiry$cohorts[[cohort_name]]
  csv <- file.path(in_dir, sprintf("cohort_%s_design.csv", cohort_name))
  if (!file.exists(csv)) {
    cat(sprintf("WARN: cohort design missing for %s, skipping\n", cohort_name))
    next
  }
  d <- read_csv(csv, show_col_types = FALSE)
  cat(sprintf("\n--- cohort %s: %d libraries ---\n", cohort_name, nrow(d)))

  donor_col <- intersect(c("ihbca_donor_id", "donor_id", "donor"), colnames(d))[1]
  study_col <- intersect(c("study", "dataset", "source"), colnames(d))[1]
  terms <- intersect(inquiry$terms_of_interest, colnames(d))

  # Landscape: per term × level counts.
  # For continuous terms, "level" is bucketed to quartiles to keep the table
  # tractable; otherwise we cast level to character for type-stable bind.
  for (term in terms) {
    vals <- d[[term]]
    is_num <- is.numeric(vals)
    if (is_num) {
      brk <- tryCatch(
        unique(quantile(vals, probs = c(0, .25, .5, .75, 1), na.rm = TRUE)),
        error = function(e) NULL)
      if (!is.null(brk) && length(brk) >= 3) {
        bucket <- cut(vals, breaks = brk, include.lowest = TRUE)
        d2 <- d
        d2$.lvl <- as.character(bucket)
      } else {
        d2 <- d
        d2$.lvl <- as.character(vals)
      }
    } else {
      d2 <- d
      d2$.lvl <- as.character(vals)
    }
    tab <- d2 %>%
      group_by(.lvl) %>%
      summarise(n_libraries = n(),
                n_donors    = n_distinct(.data[[donor_col]]),
                .groups = "drop") %>%
      mutate(cohort = cohort_name, term = term, type = if (is_num) "numeric" else "categorical") %>%
      rename(level = .lvl)
    landscape_rows[[length(landscape_rows) + 1]] <- tab
  }

  # Collinearity: pairwise association + p-value, both pooled and per-study.
  # 'within_study' = NA → pooled across whole cohort.
  studies_in_cohort <- if (!is.na(study_col))
    sort(unique(as.character(d[[study_col]]))) else character(0)
  scopes <- c(NA_character_, studies_in_cohort)

  for (sc in scopes) {
    d_sc <- if (is.na(sc)) d else d %>% filter(.data[[study_col]] == sc)
    if (nrow(d_sc) < 4) next
    for (i in seq_along(terms)) {
      for (j in seq_along(terms)) {
        if (j <= i) next
        # Skip pairs where one side is the study column when we're within-study
        if (!is.na(sc) && (terms[i] == study_col || terms[j] == study_col)) next
        a <- tryCatch(assoc(d_sc[[terms[i]]], d_sc[[terms[j]]]),
                       error = function(e) list(assoc = NA_real_, p = NA_real_,
                                                 test = paste0("err:", conditionMessage(e))))
        collin_rows[[length(collin_rows) + 1]] <- tibble(
          cohort       = cohort_name,
          within_study = sc,
          n_obs        = sum(complete.cases(d_sc[[terms[i]]], d_sc[[terms[j]]])),
          term1        = terms[i],
          term2        = terms[j],
          test         = a$test,
          association  = if (is.null(a$assoc)) NA_real_ else round(a$assoc, 3),
          pvalue       = if (is.null(a$p))     NA_real_ else signif(a$p, 3),
          identifiable = !is.na(a$assoc) && a$assoc < 0.95
        )
      }
    }
  }

  # Per-study testability
  if (!is.na(study_col)) {
    for (term in terms) {
      if (term == study_col) next
      d2 <- d
      d2$.term_chr <- as.character(d2[[term]])
      tt <- d2 %>%
        group_by(.data[[study_col]]) %>%
        summarise(
          n_levels  = n_distinct(.term_chr),
          n_donors  = n_distinct(.data[[donor_col]]),
          levels    = paste(sort(unique(.term_chr)), collapse = "|"),
          .groups = "drop"
        ) %>%
        mutate(cohort = cohort_name, term = term,
               estimable_per_study = n_levels >= 2,
               n_donors_with_var = n_donors)
      testability_rows[[length(testability_rows) + 1]] <- tt %>%
        rename(study = !!study_col)
    }
  }

  # Markdown summary section
  md_lines <- c(md_lines, sprintf("## cohort_%s", cohort_name), "",
    sprintf("- libraries: %d", nrow(d)),
    sprintf("- donors:    %d", n_distinct(d[[donor_col]])),
    "")

  md_lines <- c(md_lines, "### term landscape", "", "| term | n levels | n libraries with each level |", "|---|---|---|")
  for (term in terms) {
    tab <- table(d[[term]], useNA = "ifany")
    levs <- paste(sprintf("%s=%d", names(tab), as.integer(tab)), collapse = "; ")
    md_lines <- c(md_lines, sprintf("| %s | %d | %s |", term, length(tab), levs))
  }
  md_lines <- c(md_lines, "")

  md_lines <- c(md_lines, "### collinearity (top |association| pairs, |V|>=0.5)", "", "| term1 | term2 | assoc | identifiable |", "|---|---|---|---|")
  collin_this <- bind_rows(collin_rows[vapply(collin_rows,
                                              function(r) isTRUE(r$cohort[1] == cohort_name),
                                              logical(1))]) %>%
    arrange(desc(association))
  for (k in seq_len(min(nrow(collin_this), 30))) {
    r <- collin_this[k, ]
    if (is.na(r$association) || r$association < 0.5) next
    md_lines <- c(md_lines, sprintf("| %s | %s | %.3f | %s |",
                                    r$term1, r$term2, r$association,
                                    ifelse(r$identifiable, "yes", "**NO (aliased)**")))
  }
  md_lines <- c(md_lines, "")
}

# -- write outputs ------------------------------------------------------------
landscape_df <- bind_rows(landscape_rows)
collin_df    <- bind_rows(collin_rows)
testab_df    <- bind_rows(testability_rows)

# Flag pairs whose pooled or any per-study association is high+significant.
# This becomes the "decision-relevant" view: which pairs would alter the
# interpretation of a contrast if not handled? It's a derived filter on the
# pairwise matrix — no separate machinery.
flag_alpha <- inquiry$context$collinearity_alpha %||% 0.05
flag_eff   <- inquiry$context$collinearity_min_assoc %||% 0.30
collin_df  <- bind_rows(collin_rows) %>%
  mutate(decision_flag = !is.na(association) & !is.na(pvalue) &
                          association >= flag_eff & pvalue < flag_alpha,
         verdict = case_when(
           is.na(pvalue) ~ "INDETERMINATE",
           association >= 0.95 ~ "ALIASED (not separable)",
           decision_flag ~ "COLLINEAR (adj is meaningful sensitivity)",
           TRUE ~ "ORTHOGONAL (adj would over-correct; do NOT include)"
         ))

# Markdown summary: surface flagged pairs prominently
flagged <- collin_df %>%
  filter(decision_flag | grepl("ALIASED", verdict)) %>%
  arrange(cohort, !is.na(within_study), within_study, desc(association))
if (nrow(flagged) > 0) {
  md_lines <- c(md_lines, "## decision-relevant collinearity (flagged pairs)", "",
    sprintf("Flag rule: |assoc| >= %.2f AND p < %.2f  (per inquiry config).",
            flag_eff, flag_alpha),
    "Verdicts:",
    "- ORTHOGONAL → adjusting for the second term would OVER-CORRECT; do NOT include in formula.",
    "- COLLINEAR  → adjusting IS a meaningful sensitivity contrast.",
    "- ALIASED    → terms not separable; one must be dropped.",
    "",
    "| cohort | within_study | term1 | term2 | test | n | assoc | p | verdict |",
    "|---|---|---|---|---|---|---|---|---|")
  for (i in seq_len(nrow(flagged))) {
    r <- flagged[i, ]
    md_lines <- c(md_lines, sprintf("| %s | %s | %s | %s | %s | %d | %.3f | %g | %s |",
      r$cohort, ifelse(is.na(r$within_study), "(pooled)", r$within_study),
      r$term1, r$term2, r$test, r$n_obs, r$association, r$pvalue, r$verdict))
  }
  md_lines <- c(md_lines, "")
}

# -- write outputs ------------------------------------------------------------
write_csv(landscape_df, file.path(out_dir, "cohort_landscape.csv"))
write_csv(collin_df,    file.path(out_dir, "collinearity_matrix.csv"))
write_csv(testab_df,    file.path(out_dir, "per_study_testability.csv"))
writeLines(md_lines, file.path(out_dir, "STAGEB_covariate_landscape.md"))

cat(sprintf("\nWrote:\n  %s\n  %s\n  %s\n  %s\n",
            file.path(out_dir, "cohort_landscape.csv"),
            file.path(out_dir, "collinearity_matrix.csv"),
            file.path(out_dir, "per_study_testability.csv"),
            file.path(out_dir, "STAGEB_covariate_landscape.md")))
cat("=== Stage B done ===\n")
