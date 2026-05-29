#!/usr/bin/env Rscript
# annotate_gsea_nes.R
# Post-hoc annotation of per_nhoodgroup_nes.csv with parent_L2_joint,
# parent_compartment, n_nhoods_in_group via NhoodGroup_renamed lookup.
#
# Reason: F.3 outputs use filesystem-safe NhoodGroup IDs (epi__BMYO_basal_23)
# while F.1 uses the canonical form (epi::BMYO-basal_23). The GSEA reader
# emitted the F.3 form; this annotates against F.1 by normalizing both sides
# (gsub "::" -> "__", gsub "-" -> "_").
#
# Output: outputs/stageI3_gsea/per_nhoodgroup_nes_annotated.csv

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr)
})

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "annotate_gsea_nes.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths

CONTRASTS <- c("parity_in_AR", "parity_x_HR_BRCA1",
                "parity_x_HR_sporadic", "parity_x_HR_BRCA2")

nes_path <- file.path(paths$inquiry_root, "outputs", "stageI3_gsea",
                       "per_nhoodgroup_nes.csv")
if (!file.exists(nes_path)) {
  stop("Missing GSEA NES at: ", nes_path)
}
nes <- read_csv(nes_path, show_col_types = FALSE)

parent_lookup_rows <- list()
for (cname in CONTRASTS) {
  gs_csv <- file.path(paths$outputs$stageF1, cname, "nhood_groups_summary.csv")
  if (!file.exists(gs_csv)) next
  gs <- read_csv(gs_csv, show_col_types = FALSE) %>%
    transmute(NhoodGroup_renamed = as.character(NhoodGroup_renamed),
              parent_L2_joint = as.character(parent_L2_joint),
              parent_compartment = as.character(parent_compartment),
              n_nhoods_in_group = suppressWarnings(as.numeric(n_nhoods_in_group))) %>%
    mutate(NhoodGroup_renamed = gsub("-", "_",
                                       gsub("::", "__", NhoodGroup_renamed)),
           contrast = cname)
  parent_lookup_rows[[cname]] <- gs
}
parent_lookup <- bind_rows(parent_lookup_rows)

nes_ann <- nes %>%
  select(-any_of(c("parent_L2_joint", "parent_compartment", "n_nhoods_in_group"))) %>%
  left_join(parent_lookup, by = c("contrast", "NhoodGroup_renamed"))

cat(sprintf("Annotated rows: %d / %d (%.1f%% had parent L2 found)\n",
            sum(!is.na(nes_ann$parent_L2_joint)), nrow(nes_ann),
            100 * sum(!is.na(nes_ann$parent_L2_joint)) / nrow(nes_ann)))

out_nes <- file.path(paths$inquiry_root, "outputs", "stageI3_gsea",
                      "per_nhoodgroup_nes_annotated.csv")
write_csv(nes_ann, out_nes)
cat(sprintf("Wrote: %s\n", out_nes))
