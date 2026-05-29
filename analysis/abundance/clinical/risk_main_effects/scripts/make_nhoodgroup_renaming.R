#!/usr/bin/env Rscript
# make_nhoodgroup_renaming.R
# Build a NhoodGroup renaming lookup with size_rank counted SEPARATELY per parent_L2.
#
# Per Kai/Noah convention 2026-05-07:
#   "counts separate per L2 (Fibro-1, Fibro-2, basal 1, basal 2, etc.)"
#
# For each (contrast x parent_L2_joint):
#   - Filter to viable NhoodGroups (n_nhoods_in_group >= 10; matches F.3 limma-voom floor)
#   - Rank descending by n_nhoods_in_group
#   - Assign size_rank_within_L2 = 1, 2, 3, ...
#
# Emits:
#   outputs/nhoodgroup_renaming/lookup.csv
#
# Columns:
#   contrast, stratum, parent_L2_joint, parent_compartment, parent_label,
#   NhoodGroup, NhoodGroup_renamed (existing F.1 form: epi::BMYO-basal_11),
#   n_nhoods_in_group, n_sig, group_med_lfc, group_pct_up,
#   viable, size_rank_within_L2 (NA for non-viable),
#   name_short ([parent_label]_[rank], e.g. BMYO-basal_1),
#   name_full  ([parent_label]_[rank]_[stratum], e.g. BMYO-basal_1_BR1)

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr)
})

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "make_nhoodgroup_renaming.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths
inq <- paths$inquiry_root

VIABLE_MIN <- 10L

# Stratum codes. AR/HR-stratum within-cohort parity main-effect tests (the 5 we care
# about) get short codes; older `parity_in_HR_*` redundancies get an `m` suffix so
# they don't collide.
stratum_map <- c(
  "parity_in_AR"          = "AR",
  "parity_x_HR_BRCA1"     = "BR1",
  "parity_x_HR_BRCA2"     = "BR2",
  "parity_x_HR_BRCA12"    = "BR12",
  "parity_x_HR_sporadic"  = "HRS",
  "parity_in_HR_BRCA1"    = "BR1m",
  "parity_in_HR_BRCA2"    = "BR2m",
  "parity_in_HR_sporadic" = "HRSm"
)

nhg_dir <- file.path(inq, "outputs", "stageF1_nhoodgroups")
contrasts <- list.dirs(nhg_dir, full.names = FALSE, recursive = FALSE)
contrasts <- contrasts[contrasts != ""]
cat(sprintf("Contrasts found (%d): %s\n",
            length(contrasts), paste(contrasts, collapse = ", ")))

all_rows <- list()
for (cn in contrasts) {
  fp <- file.path(nhg_dir, cn, "nhood_groups_summary.csv")
  if (!file.exists(fp)) {
    cat(sprintf("  MISSING: %s\n", fp)); next
  }
  d <- read_csv(fp, show_col_types = FALSE)
  if (nrow(d) == 0) {
    cat(sprintf("  EMPTY (header-only, skipped): %s\n", fp)); next
  }
  d <- d %>%
    mutate(
      NhoodGroup          = as.character(NhoodGroup),
      n_nhoods_in_group   = as.integer(n_nhoods_in_group),
      n_sig               = as.integer(n_sig),
      group_med_lfc       = as.numeric(group_med_lfc),
      group_pct_up        = as.numeric(group_pct_up),
      parent_L2_joint     = as.character(parent_L2_joint),
      parent_compartment  = as.character(parent_compartment),
      parent_label        = as.character(parent_label),
      NhoodGroup_renamed  = as.character(NhoodGroup_renamed),
      contrast            = cn
    )
  all_rows[[cn]] <- d
}
if (length(all_rows) == 0) stop("No non-empty F.1 summaries found.")
df <- bind_rows(all_rows)
cat(sprintf("Read %d total NhoodGroup-summary rows from %d contrasts.\n",
            nrow(df), length(all_rows)))

# Resolve unmapped contrasts loudly (rather than silently producing NA stratum)
unmapped <- setdiff(unique(df$contrast), names(stratum_map))
if (length(unmapped) > 0) {
  warning("Unmapped contrasts (no stratum code): ",
          paste(unmapped, collapse = ", "))
}

df <- df %>%
  mutate(
    stratum = unname(stratum_map[contrast]),
    viable  = n_nhoods_in_group >= VIABLE_MIN
  )

# Rank within (contrast x parent_L2_joint). Viable rows take consecutive ranks 1..k;
# non-viable rows keep NA. cumsum(viable) gives the rank only on viable rows because
# we mask the non-viable cells to NA afterwards.
ranked <- df %>%
  group_by(contrast, parent_L2_joint) %>%
  arrange(desc(n_nhoods_in_group), NhoodGroup, .by_group = TRUE) %>%
  mutate(
    rank_running        = cumsum(as.integer(viable)),
    size_rank_within_L2 = if_else(viable, rank_running, NA_integer_)
  ) %>%
  ungroup() %>%
  select(-rank_running) %>%
  mutate(
    name_short = if_else(!is.na(size_rank_within_L2),
                          sprintf("%s_%d", parent_label, size_rank_within_L2),
                          NA_character_),
    name_full  = if_else(!is.na(size_rank_within_L2),
                          sprintf("%s_%d_%s", parent_label,
                                  size_rank_within_L2, stratum),
                          NA_character_)
  )

ranked_out <- ranked %>%
  select(contrast, stratum, parent_L2_joint, parent_compartment, parent_label,
         NhoodGroup, NhoodGroup_renamed,
         n_nhoods_in_group, n_sig, group_med_lfc, group_pct_up,
         viable, size_rank_within_L2, name_short, name_full)

out_dir <- file.path(inq, "outputs", "nhoodgroup_renaming")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(out_dir, "lookup.csv")
write_csv(ranked_out, out_path)
cat(sprintf("Wrote: %s\n  rows=%d  viable=%d  non-viable=%d\n",
            out_path, nrow(ranked_out),
            sum(ranked_out$viable), sum(!ranked_out$viable)))

# ---- Sanity prints ----
cat("\n=== Spot check: BMYO-basal in parity_x_HR_BRCA1 ===\n")
print(ranked_out %>%
        filter(contrast == "parity_x_HR_BRCA1",
               parent_L2_joint == "epi::BMYO-basal") %>%
        select(NhoodGroup, n_nhoods_in_group, group_med_lfc,
               viable, size_rank_within_L2, name_short, name_full),
      n = Inf)

cat("\n=== Spot check: Fibro-SFRP4 in parity_in_AR ===\n")
print(ranked_out %>%
        filter(contrast == "parity_in_AR",
               parent_L2_joint == "str::Fibro-SFRP4") %>%
        select(NhoodGroup, n_nhoods_in_group, group_med_lfc,
               viable, size_rank_within_L2, name_short, name_full),
      n = Inf)

cat("\n=== Per-(contrast x parent_L2) viable group counts (top 20 by viable n) ===\n")
print(ranked_out %>%
        filter(viable) %>%
        count(contrast, parent_L2_joint, name = "n_viable") %>%
        arrange(desc(n_viable)) %>%
        head(20))

cat("\n=== make_nhoodgroup_renaming done ===\n")
