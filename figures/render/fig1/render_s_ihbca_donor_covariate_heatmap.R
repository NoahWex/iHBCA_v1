#!/usr/bin/env Rscript
# Render Fig 1 supplemental: donor covariate landscape heatmap.
#
# 287-donor matrix showing harmonized clinical, risk-detail, and demographic
# covariates, grouped by contributing study. Each row is a donor; each column
# is a covariate. Color encodes the donor's category for that covariate;
# missing entries are gray. A study row-annotation strip and a continuous-age
# barplot annotation accompany the discrete-covariate matrix.
#
# Display logic
#   Donors are sorted by study, then by descending continuous age (where
#   available). Donors with age_binary but no age_continuous are imputed at
#   {young = 35, old = 60} for the age barplot and shaded in a distinct color
#   to flag the imputation. The heatmap rows are split by study to make
#   per-cohort coverage and aliasing visible at a glance.
#
# Color encoding (as_built constraint per manifest fig1_supps_promotion_20260513)
#   Inline per-variable semantic palettes are retained for this panel because
#   the binary contrasts (BR1 vs AR, parous vs nulliparous, pre vs post,
#   high-risk vs average-risk) carry conventional clinical color mappings
#   that downstream readers are accustomed to. Aesthetics-framework refactor
#   into aesthetics.yaml is a deferred follow-up. The companion alluvial
#   panel (s_ihbca_study_alluvials.R) uses the framework Tol pool for its
#   single-color cycling pattern; this heatmap uses semantic palettes.
#
# Inputs (resolved via argparse)
#   --metadata: harmonized_donor_metadata.csv (one row per donor, harmonized covariates)
#   --aesthetics: optional aesthetics.yaml path (currently unused; reserved
#                 for future framework refactor)
#   --output-dir: destination directory
#
# Output
#   <output-dir>/s1_1_covariate_heatmap.pdf
#   <output-dir>/covariate_coverage.csv (per-covariate n_known + pct_known + n_levels)

suppressPackageStartupMessages({
  library(argparse)
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(yaml)
})

parser <- ArgumentParser(description = "Donor covariate landscape heatmap")
parser$add_argument("--metadata", required = TRUE,
                    help = "harmonized_donor_metadata.csv")
parser$add_argument("--aesthetics", default = NULL,
                    help = "Optional aesthetics.yaml path (reserved for future framework refactor)")
parser$add_argument("--output-dir", required = TRUE)
args <- parser$parse_args()

dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loading metadata:", args$metadata, "\n")
donors <- read.csv(args$metadata, stringsAsFactors = FALSE)
cat("  ", nrow(donors), "donors\n")

# Sort donors by study, then age_binary (young -> old, missing last), then by
# age_continuous (ascending). Donors with no age_binary value sort to the end
# of their study block; within each binary stratum, ascending continuous age.
donors <- donors %>%
  arrange(
    study,
    factor(age_binary, levels = c("young", "old")),
    as.numeric(age_continuous)
  )

# CSV stores parity_binary and menopausal_status_binary as integer codes (0/1);
# the palettes below key on semantic strings. Recode before palette lookup so
# values map to a color instead of falling through to "missing".
if ("parity_binary" %in% colnames(donors)) {
  donors$parity_binary <- dplyr::recode(as.character(donors$parity_binary),
                                        "0" = "nulliparous",
                                        "1" = "parous",
                                        .default = NA_character_)
}
if ("menopausal_status_binary" %in% colnames(donors)) {
  donors$menopausal_status_binary <- dplyr::recode(as.character(donors$menopausal_status_binary),
                                                   "0" = "pre",
                                                   "1" = "post",
                                                   .default = NA_character_)
}

# Covariate blocks, displayed left-to-right in the order: clinical contrasts,
# risk detail, demographic context.
clinical    <- c("age_binary", "parity_binary", "risk_status_binary",
                 "menopausal_status_binary")
risk_detail <- c("brca_genotype", "cancer_history", "tissue_indication")
demographic <- c("bmi_category", "ethnicity_grouped")

all_covs <- c(clinical, risk_detail, demographic)
all_covs <- all_covs[all_covs %in% colnames(donors)]

# Per-variable semantic color palettes (as_built constraint; see header note).
palettes <- list(
  age_binary               = c("young" = "#4DBBD5", "old" = "#E64B35"),
  parity_binary            = c("nulliparous" = "#91D1C2", "parous" = "#7E6148"),
  risk_status_binary       = c("AR" = "#3C5488", "HR" = "#DC0000"),
  menopausal_status_binary = c("pre" = "#F39B7F", "post" = "#8491B4"),
  brca_genotype            = c("negative" = "#00A087", "BRCA1" = "#DC0000",
                               "BRCA2" = "#E64B35", "RAD51C" = "#B09C85"),
  cancer_history           = c("no" = "#00A087", "yes" = "#DC0000"),
  tissue_indication        = c("reduction" = "#3C5488", "prophylactic" = "#DC0000",
                               "contralateral" = "#E64B35"),
  bmi_category             = c("underweight" = "#91D1C2", "normal" = "#00A087",
                               "overweight" = "#F39B7F", "obese" = "#E64B35"),
  ethnicity_grouped        = c("white" = "#4DBBD5", "black" = "#E64B35",
                               "hispanic" = "#00A087", "asian" = "#F39B7F",
                               "jewish" = "#8491B4", "other" = "#B09C85")
)

study_colors <- c(
  "gray" = "#E64B35", "kumar" = "#4DBBD5", "murrow" = "#00A087",
  "nee" = "#3C5488", "pal" = "#F39B7F", "reed" = "#8491B4",
  "twigger" = "#7E6148"
)

# Coerce missing/unknown to "missing" so the gray-fill encoding reads cleanly.
anno_df     <- data.frame(row.names = donors$ihbca_donor_id)
anno_colors <- list()
for (cov in all_covs) {
  vals <- donors[[cov]]
  vals[vals %in% c("unknown", "nan", "", NA)] <- "missing"
  anno_df[[cov]] <- vals
  pal <- palettes[[cov]]
  if (!is.null(pal)) {
    anno_colors[[cov]] <- c(pal, "missing" = "#CCCCCC")
  }
}

# Continuous age barplot: actual values where available, imputed from
# age_binary at {young = 35, old = 60} otherwise. Bar color flags the source:
# blue = actual, orange = imputed-from-binary, gray = missing.
age_cont    <- as.numeric(donors$age_continuous)
age_imputed <- rep(FALSE, nrow(donors))
for (i in seq_len(nrow(donors))) {
  if (is.na(age_cont[i]) && donors$age_binary[i] %in% c("young", "old")) {
    age_cont[i]    <- ifelse(donors$age_binary[i] == "young", 35, 60)
    age_imputed[i] <- TRUE
  }
}
age_bar_color <- ifelse(is.na(age_cont), "#CCCCCC",
                        ifelse(age_imputed, "#F39B7F", "#4DBBD5"))

cat("Generating heatmap...\n")
cat("  Age: ", sum(!is.na(age_cont) & !age_imputed), " actual, ",
    sum(age_imputed), " imputed from binary, ",
    sum(is.na(age_cont)), " missing\n", sep = "")

pdf(file.path(args$output_dir, "s1_1_covariate_heatmap.pdf"),
    width = 12, height = 10)

# Row annotation: study category + continuous age barplot (blue/orange/gray).
row_ha <- rowAnnotation(
  study = donors$study,
  age = anno_barplot(ifelse(is.na(age_cont), 0, age_cont),
                     width = unit(2, "cm"),
                     gp = gpar(fill = age_bar_color, col = NA)),
  col = list(study = study_colors),
  annotation_name_gp = gpar(fontsize = 7),
  annotation_legend_param = list(study = list(title_gp = gpar(fontsize = 7),
                                              labels_gp = gpar(fontsize = 6)))
)

study_split <- factor(donors$study, levels = names(study_colors))

# Build one Heatmap per covariate column; concatenate into the row layout.
ht_list <- NULL
first <- TRUE
for (cov in all_covs) {
  vals <- anno_df[[cov]]
  pal  <- anno_colors[[cov]]
  if (is.null(pal)) next

  col_fun <- pal[unique(vals)]
  col_fun <- col_fun[!is.na(names(col_fun))]

  ht_args <- list(
    matrix = matrix(vals, ncol = 1, dimnames = list(NULL, cov)),
    name = cov,
    col = col_fun,
    width = unit(8, "mm"),
    show_row_names = FALSE,
    show_column_names = TRUE,
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 7),
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 7),
      labels_gp = gpar(fontsize = 6)
    ),
    border = TRUE
  )
  if (first) {
    ht_args$row_split    <- study_split
    ht_args$row_title_rot <- 0
    ht_args$row_title_gp  <- gpar(fontsize = 7)
    ht_args$row_gap       <- unit(1, "mm")
  }
  ht <- do.call(Heatmap, ht_args)
  first <- FALSE

  if (is.null(ht_list)) {
    ht_list <- ht
  } else {
    ht_list <- ht_list + ht
  }
}

draw(row_ha + ht_list, heatmap_legend_side = "bottom")
dev.off()
cat("Saved: s1_1_covariate_heatmap.pdf\n")

# Per-covariate coverage summary (n + % donors with non-missing values + n_levels).
coverage <- data.frame(
  covariate  = all_covs,
  n_known    = sapply(all_covs, function(c) sum(anno_df[[c]] != "missing")),
  pct_known  = sapply(all_covs, function(c) round(mean(anno_df[[c]] != "missing") * 100, 1)),
  n_levels   = sapply(all_covs, function(c) length(unique(anno_df[[c]][anno_df[[c]] != "missing"])))
)
write.csv(coverage, file.path(args$output_dir, "covariate_coverage.csv"),
          row.names = FALSE)
cat("Saved: covariate_coverage.csv\n")
cat("Done.\n")
