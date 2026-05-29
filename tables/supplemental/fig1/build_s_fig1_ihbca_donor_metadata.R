#!/usr/bin/env Rscript
# Build Supplementary Table 1: harmonized per-donor metadata.
#
# 287 donors x 23 covariates harmonized across 7 contributing iHBCA cohorts.
# Companion table to Supplementary Fig. 1 (donor covariate heatmap).
# Visualized via render_s_ihbca_donor_covariate_heatmap.R (publication/figures/render/ihbca/).
#
# Upstream assembly:
#   per-study entry sheets
#     -> iHBCAv1_upload/external_studies/harmonization/scripts/build_harmonized_donors.py
#     -> iHBCAv1_upload/external_studies/outputs/harmonized_metadata/harmonized_donor_metadata.csv
#     -> iHBCAv1_upload/publication/scripts/stage_l1_metadata.py
#     -> iHBCAv1_upload/publication/config/metadata_stages/L1_harmonized_donor.csv  (consumed here)
#
# This script is a thin pass-through: load the staged canonical, validate,
# and write to publication/tables/supplemental/fig1/.

suppressPackageStartupMessages({
  library(argparse)
  library(yaml)
})

parser <- ArgumentParser(description = "Build Supp Table 1: harmonized donor metadata.")
parser$add_argument("--project-root", required = TRUE,
                    help = "Repository root for path resolution via publication/config/load_paths.R")
parser$add_argument("--out-path", default = NULL,
                    help = "Override default output path (default: publication/tables/supplemental/fig1/s_fig1_ihbca_donor_metadata.csv)")
args <- parser$parse_args()

# Resolve paths via repo path resolver.
source(file.path(args$project_root, "publication/config/load_paths.R"))
paths <- load_paths(file.path(args$project_root, "publication/config"))

# Upstream staged canonical (paths.yaml sources.donor_metadata).
source_csv <- resolve_path("donor_metadata", paths)
if (!file.exists(source_csv)) {
  stop("Source CSV not found at resolved path: ", source_csv)
}

donors <- read.csv(source_csv, stringsAsFactors = FALSE)

# Validation gates.
required_cols <- c(
  "ihbca_donor_id", "study",
  "age_continuous", "age_binary",
  "parity_count", "parity_binary", "age_at_first_birth",
  "brca_genotype", "cancer_history", "tissue_indication",
  "risk_status_binary", "risk_genotype_only",
  "menopausal_status_detailed", "menopausal_status_binary",
  "ethnicity_verbatim", "ethnicity_grouped",
  "bmi_continuous", "bmi_category",
  "sample_preservation", "sample_type", "facs_status",
  "dissociation_minutes", "metadata_notes"
)
missing_cols <- setdiff(required_cols, colnames(donors))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

if (nrow(donors) != 287) {
  stop("Expected 287 donors; got ", nrow(donors))
}
if (any(is.na(donors$ihbca_donor_id))) stop("ihbca_donor_id has NAs")
if (any(duplicated(donors$ihbca_donor_id))) stop("ihbca_donor_id has duplicates")
if (any(is.na(donors$study))) stop("study has NAs")
if (any(is.na(donors$risk_status_binary))) stop("risk_status_binary has NAs (prose claim: 100% coverage)")

expected_studies <- c("gray", "kumar", "murrow", "nee", "pal", "reed", "twigger")
extra_studies <- setdiff(unique(donors$study), expected_studies)
if (length(extra_studies) > 0) {
  warning("Unexpected study values: ", paste(extra_studies, collapse = ", "))
}

# Coverage report (per §1 prose claims).
coverage <- function(col) sum(!is.na(donors[[col]])) / nrow(donors) * 100
cat(sprintf("Coverage report (n = %d donors):\n", nrow(donors)))
cat(sprintf("  age_continuous OR age_binary: %5.1f%%\n",
            sum(!is.na(donors$age_continuous) | !is.na(donors$age_binary)) / nrow(donors) * 100))
cat(sprintf("  parity_binary:                %5.1f%%\n", coverage("parity_binary")))
cat(sprintf("  menopausal_status_binary:     %5.1f%%\n", coverage("menopausal_status_binary")))
cat(sprintf("  brca_genotype:                %5.1f%%\n", coverage("brca_genotype")))
cat(sprintf("  risk_status_binary:           %5.1f%%\n", coverage("risk_status_binary")))

# Reorder columns to schema order.
donors <- donors[, required_cols]

# Write.
out_path <- if (!is.null(args$out_path)) args$out_path else file.path(
  args$project_root,
  "publication/tables/supplemental/fig1/s_fig1_ihbca_donor_metadata.csv"
)
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
write.csv(donors, out_path, row.names = FALSE)
cat("Wrote", nrow(donors), "donors x", ncol(donors), "covariates to:\n  ", out_path, "\n")
