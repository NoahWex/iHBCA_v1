#!/usr/bin/env Rscript
# 39_stageH_assembly.R
# Stage H — assembly + reporting.
# Joins all upstream stage outputs into the final `findings_table.csv`.
# Two row-level keys:
#   - per-L2 row:           one row per (L2_joint × contrast)
#   - per-nhoodgroup row:   one row per (NhoodGroup_renamed × contrast)
#
# Output:
#   stageH_findings_table.csv
#   reports/findings_report.md  (auto-generated narrative skeleton)
#   reports/methods.md           (parameterized methods text)

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(yaml); library(purrr)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "39_stageH_assembly.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths
inquiry <- cfg$inquiry

out_csv <- paths$outputs$stageH
ensure_dir(dirname(out_csv))
ensure_dir(paths$outputs$reports)

cat(sprintf("=== Stage H: assembly for '%s' ===\n", inquiry$inquiry))

# -- load upstream outputs ---------------------------------------------------
e_wide <- read_csv(paths$outputs$stageE, show_col_types = FALSE)
g_l2   <- read_csv(paths$outputs$stageG, show_col_types = FALSE)

ng_ctx_path <- file.path(dirname(paths$outputs$stageG),
                          "stageG_per_nhoodgroup_context.csv")
g_ng <- if (file.exists(ng_ctx_path)) {
          read_csv(ng_ctx_path, show_col_types = FALSE)
        } else {
          tibble()
        }

# Aggregate Stage F.1 nhoodgroup summaries across all contrasts
ng_rows <- list()
contrast_dirs <- list.dirs(paths$outputs$stageF1, recursive = FALSE)
for (cdir in contrast_dirs) {
  cname <- basename(cdir)
  gs_csv <- file.path(cdir, "nhood_groups_summary.csv")
  if (!file.exists(gs_csv)) next
  gs <- read_csv(gs_csv, show_col_types = FALSE) %>%
    mutate(contrast = cname)
  # Harmonize column types — readr infers per-file and can disagree across
  # contrasts (e.g., when one contrast has all-NA / no rows for some column).
  char_cols <- c("NhoodGroup", "compartment", "label", "L2_joint",
                 "NhoodGroup_renamed", "parent_L2_joint",
                 "parent_compartment", "parent_label")
  num_cols  <- c("n_nhoods_in_group", "n_sig", "n_sig_up", "n_sig_dn",
                 "group_med_lfc", "group_pct_up")
  for (cc in intersect(char_cols, colnames(gs))) gs[[cc]] <- as.character(gs[[cc]])
  for (cc in intersect(num_cols,  colnames(gs))) gs[[cc]] <- as.numeric(gs[[cc]])
  ng_rows[[length(ng_rows) + 1]] <- gs
}
ng_summary <- bind_rows(ng_rows)

# -- assemble per-L2 findings ------------------------------------------------
# Wide stage E + Stage G per-L2 context
findings_l2 <- e_wide %>%
  left_join(g_l2, by = c("L2_joint")) %>%
  mutate(row_type = "L2")
# Note: L2 rows may collapse contrasts; the user can pivot back if desired.

# Per-nhoodgroup findings: nhood_groups_summary × Stage G nhoodgroup context
findings_ng <- ng_summary %>%
  left_join(g_ng, by = c("contrast", "NhoodGroup_renamed" = "group_renamed")) %>%
  mutate(row_type = "nhoodgroup")

# -- write outputs -----------------------------------------------------------
write_csv(findings_l2, file.path(dirname(out_csv), "stageH_findings_per_L2.csv"))
write_csv(findings_ng, file.path(dirname(out_csv), "stageH_findings_per_nhoodgroup.csv"))

# Combined wide (long along contrast dimension)
findings_combined <- bind_rows(
  findings_l2 %>% mutate(group_renamed = NA_character_),
  findings_ng %>% rename(L2_joint = parent_L2_joint) %>%
                   mutate(L2_joint = L2_joint)
)
write_csv(findings_combined, out_csv)
cat(sprintf("Wrote: %s  (rows=%d cols=%d)\n",
            out_csv, nrow(findings_combined), ncol(findings_combined)))

# -- reports/findings_report.md (narrative skeleton) -------------------------
md <- c(
  sprintf("# Findings — inquiry %s", inquiry$inquiry),
  sprintf("Generated %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  sprintf("- Description: %s", inquiry$description %||% ""),
  sprintf("- Cohorts: %s", paste(names(inquiry$cohorts), collapse = ", ")),
  sprintf("- Contrasts run: %d", n_distinct(findings_l2 %>% select(contrast) %>% unlist())),
  "",
  "## Top L2 by main effect",
  ""
)
# findings_l2 is wide-format (per-contrast columns); n_sig as a single column
# only exists if Stage E was pivoted long earlier. Guard accordingly so the
# narrative top-table is best-effort, never blocking on schema.
top_main <- if (all(c("contrast", "n_sig") %in% colnames(findings_l2))) {
  findings_l2 %>%
    filter(grepl("main", contrast, ignore.case = TRUE)) %>%
    arrange(desc(n_sig)) %>%
    head(15)
} else {
  tibble()
}
if (nrow(top_main) > 0) {
  md <- c(md, "| L2 | n_sig | med_lfc | pct_up | direction concordance |",
              "|---|---|---|---|---|")
  for (k in seq_len(nrow(top_main))) {
    r <- top_main[k, ]
    md <- c(md, sprintf("| %s | %d | %.3f | %.1f%% | %s |",
                        r$L2_joint %||% "?", r$n_sig %||% NA,
                        r$med_lfc %||% NA, r$pct_up %||% NA,
                        if (!is.null(r$direction_concordance))
                          sprintf("%.2f", r$direction_concordance)
                        else ""))
  }
}
md <- c(md, "", "## Dissoc-orthogonality verdict (Stage B)", "",
  "See `outputs/stageB_exploration/parity_dissoc_orthogonality.csv` for the per-cohort,",
  "per-study orthogonality test that determined whether dissoc-adjustment would over-correct.",
  "")

writeLines(md, file.path(paths$outputs$reports, "findings_report.md"))

# -- reports/methods.md (parameterized) --------------------------------------
methods <- c(
  sprintf("# Methods — inquiry %s", inquiry$inquiry),
  "",
  "## Cohort definitions"
)
for (cname in names(inquiry$cohorts)) {
  spec <- inquiry$cohorts[[cname]]
  methods <- c(methods,
    sprintf("### %s", cname),
    "",
    sprintf("- description: %s", spec$description %||% ""),
    sprintf("- formula: `%s`", spec$formula),
    sprintf("- contrasts: %s",
            paste(sapply(spec$contrasts, function(c) c$name), collapse = ", ")),
    sprintf("- spatial FDR threshold: %s",
            spec$spatial_fdr %||% 0.05),
    "")
}
methods <- c(methods,
  "## Pipeline stages",
  "Per `da_framework.md`, this inquiry passes through Stages A–H:",
  "",
  "- A: cohort design (`30_stageA_cohort.R`)",
  "- B: covariate analysis + parity-dissoc orthogonality (`31_stageB_covariate_analysis.R`)",
  "- C: contrast specification (`32_stageC_contrasts.R`)",
  "- D: testNhoods array (`33_stageD_da_array.R`)",
  "- E: per-L2 effect summary (`34_stageE_per_l2_summary.R`)",
  "- F.1: groupNhoods array (`35_stageF1_groupnhoods_array.R`)",
  "- F.2: pseudobulk array (`36_stageF2_pseudobulk_array.R`)",
  "- F.3: within-parent markers array (`37_stageF3_within_parent_markers_array.R`)",
  "- G: per-finding context (caveats are descriptive columns) (`38_stageG_per_finding_context.R`)",
  "- H: assembly + reports (`39_stageH_assembly.R`)",
  "",
  "## Caveat framework",
  "",
  "Caveats are **descriptive context columns**, not pre-filters. Each row of `findings_table.csv`",
  "carries: per-study contribution + direction concordance, parity-dissoc orthogonality verdict",
  "(study-level), within-parent stress-marker fraction (group-level), cohort-imbalance descriptors.",
  "",
  "Dissoc-correlation does NOT invalidate a parity finding **unless** dissoc covaries with parity",
  "in the dissoc-stratified study. The Stage B orthogonality test is the gating check; Stage F.3",
  "stress-marker fraction distinguishes recovery-rate (clean markers) from state-distortion",
  "(stress-flagged markers) when collinearity exists."
)
writeLines(methods, file.path(paths$outputs$reports, "methods.md"))

cat("=== Stage H done ===\n")
