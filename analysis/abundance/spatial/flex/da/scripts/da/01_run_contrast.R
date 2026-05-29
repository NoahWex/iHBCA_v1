#!/usr/bin/env Rscript
# Phase 1, step 01 — single-contrast Milo DA runner.
#
# Tier-aware YAML-driven runner. Reads experiment_metadata.yaml for:
#   contrast definition (factor_col, levels, libraries_filter)
#   tier formula (M0/M1/M2)
#   scope -> Milo path
#   testNhoods.reducedDim (embedding key in the Milo's reducedDims slot)
#
# Locked model skeleton (FLEX-only): ~ patient_id + log_n_cells + <terms> + group
# Sign convention: levels[2] = positive logFC direction (smaller / minority group).
#
# Direct lineage: joint_anatomic_variation_20260424/scripts/da/01_run_contrast.R
#   Adaptations for FLEX-only:
#     1. reducedDim is read from YAML (scVI_100 here, Concord in joint)
#     2. Scope names use flex_* prefix (flex_full, flex_Epithelial, ...)
#     3. No platform covariate (single-platform substrate)
#
# Pattern source:
#   anatomic_variation_analysis_20260323/scripts/2_position_da.Rmd  (testNhoods + design)
#
# Usage:
#   Rscript scripts/da/01_run_contrast.R \
#     --metadata experiment_metadata.yaml \
#     --bridge   outputs/shared/flex_l1p5.csv \
#     --covars   outputs/shared/flex_library_covariates.csv \
#     --scope    flex_full \
#     --contrast C1_p1p2_vs_p3 \
#     --tier     M0 \
#     --out-dir  outputs/da/C1_p1p2_vs_p3/flex_full/M0

suppressPackageStartupMessages({
  library(argparse)
  library(miloR)
  library(SingleCellExperiment)
  library(edgeR)
  library(dplyr)
  library(yaml)
  library(data.table)
})

# --- helpers --------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

derive_vertical <- function(pid) {
  ifelse(grepl("UI|UO|Upper", pid), "Upper",
  ifelse(grepl("LI|LO|Lower", pid), "Lower", NA_character_))
}

derive_axis_amp <- function(pid) {
  ifelse(grepl("P3_A_", pid), "A",
  ifelse(grepl("P3_M_", pid), "M",
  ifelse(grepl("P3_P_", pid), "P", NA_character_)))
}

build_factor_col <- function(covars, factor_col) {
  if (factor_col == "p1p2_vs_p3") {
    ifelse(covars$p_level %in% c("P1", "P2"), "P1P2",
    ifelse(covars$p_level == "P3", "P3", NA))
  } else if (factor_col == "p_level") {
    as.character(covars$p_level)
  } else if (factor_col == "vertical") {
    derive_vertical(covars$library_id)
  } else if (factor_col == "axis_amp") {
    derive_axis_amp(covars$library_id)
  } else if (factor_col == "is_uoq") {
    tolower(as.character(covars$is_uoq))   # "true" / "false"
  } else if (factor_col == "epi_class") {
    # Externally-classified factor — assignments loaded from class_assignments_yaml
    # and joined into covars before this function is called. Reading from the
    # YAML keeps the binarization choice (median split, cut value) auditable
    # outside the runner.
    if (!"epi_class" %in% colnames(covars)) {
      stop("factor_col='epi_class' requires class_assignments_yaml in contrast spec")
    }
    as.character(covars$epi_class)
  } else {
    stop("Unknown factor_col: ", factor_col)
  }
}

apply_libraries_filter <- function(covars, expr) {
  if (is.null(expr) || identical(expr, "null")) return(covars)
  expr <- gsub("AND", "&", expr, fixed = TRUE)
  # Quote bareword literals after `==`: "p_level == P3" -> "p_level == 'P3'"
  expr <- gsub("==\\s*([A-Za-z_][A-Za-z0-9_]*)",
               "== '\\1'", expr, perl = TRUE)
  # Translate `in [A, M, ...]` -> `%in% c('A','M', ...)` one match at a time
  while (grepl("\\bin\\s*\\[", expr, perl = TRUE)) {
    m <- regmatches(expr, regexpr("\\bin\\s*\\[[^]]*\\]", expr, perl = TRUE))
    vals <- regmatches(m, regexpr("\\[[^]]*\\]", m))
    vals <- gsub("\\[|\\]", "", vals, perl = TRUE)
    vals <- trimws(strsplit(vals, ",")[[1]])
    replacement <- paste0("%in% c(", paste0("'", vals, "'", collapse = ","), ")")
    expr <- sub("\\bin\\s*\\[[^]]*\\]", replacement, expr, perl = TRUE)
  }
  keep <- with(covars, eval(parse(text = expr)))
  covars[keep, , drop = FALSE]
}

annotate_nhoods <- function(milo, cell_labels) {
  # cell_labels: named vector keyed by cell_id -> string label (or NA)
  nhood_mat <- nhoods(milo)
  n <- ncol(nhood_mat)
  labs <- character(n); fracs <- numeric(n)
  cl <- cell_labels[match(colnames(milo), names(cell_labels))]
  for (i in seq_len(n)) {
    idx <- which(nhood_mat[, i] == 1)
    if (length(idx) == 0) { labs[i] <- NA; fracs[i] <- NA; next }
    tab <- sort(table(cl[idx]), decreasing = TRUE)
    if (length(tab) == 0) { labs[i] <- NA; fracs[i] <- NA; next }
    labs[i] <- names(tab)[1]
    fracs[i] <- as.integer(tab[1]) / sum(tab)
  }
  list(label = labs, frac = round(fracs, 3))
}

# --- main -----------------------------------------------------------------

parser <- ArgumentParser()
parser$add_argument("--metadata", required = TRUE)
parser$add_argument("--bridge",   required = TRUE)
parser$add_argument("--covars",   required = TRUE)
parser$add_argument("--scope",    required = TRUE)
parser$add_argument("--contrast", required = TRUE)
parser$add_argument("--tier",     required = TRUE, choices = c("M0", "M1", "M2"))
parser$add_argument("--out-dir",  required = TRUE)
parser$add_argument("--label-column", default = NULL,
                    help = "Override metadata vocabulary.label_column (e.g. 'L2_assignment' for Track C swap)")
args <- parser$parse_args()

dir.create(args$out_dir, showWarnings = FALSE, recursive = TRUE)

message("=== 01_run_contrast.R ===")
message("scope:    ", args$scope)
message("contrast: ", args$contrast)
message("tier:     ", args$tier)

# Resolve YAML
meta <- yaml::read_yaml(args$metadata)
contrast_spec <- meta$contrasts[[args$contrast]]
stopifnot("contrast not found in YAML"  = !is.null(contrast_spec))
stopifnot("scope not in metadata$milos" = !is.null(meta$milos[[args$scope]]))
tier_spec <- contrast_spec$tiers[[args$tier]]
stopifnot("tier not found for contrast" = !is.null(tier_spec))

milo_path_raw <- meta$milos[[args$scope]]$path
if (substr(milo_path_raw, 1, 1) == "/") {
  milo_path <- milo_path_raw
} else {
  milo_path <- normalizePath(file.path(dirname(args$metadata), milo_path_raw),
                             mustWork = FALSE)
}
fdr_w        <- meta$testNhoods$fdr_weighting %||% "graph-overlap"
reduced_dim  <- meta$testNhoods$reducedDim    %||% "scVI_100"
formula_str  <- tier_spec$formula
label_column <- args$label_column %||% (meta$vocabulary$label_column %||% "l1p5_label")
label_set    <- meta$vocabulary$label_set %||% "L1.5"
message("milo:        ", milo_path)
message("formula:     ", formula_str)
message("reduced_dim: ", reduced_dim)
message("label_set:   ", label_set, " (column: ", label_column, ")")

# --- load bridge + covariates ---
bridge <- as.data.frame(fread(args$bridge))
covars <- as.data.frame(fread(args$covars))
message("bridge: ", nrow(bridge), " cells | covars: ", nrow(covars), " libraries")

stopifnot("bridge missing cell_id column"        = "cell_id" %in% colnames(bridge),
          "bridge missing label_column"          = label_column %in% colnames(bridge))

# Compute log_n_cells (locked technical covariate)
covars$log_n_cells <- log10(covars$n_cells)

# If the contrast declares an external class-assignment YAML, load and join
# into covars by library_id BEFORE build_factor_col. This pattern keeps the
# binarization (e.g. median split for epi_class) auditable outside the runner.
if (!is.null(contrast_spec$class_assignments_yaml)) {
  cay_raw <- contrast_spec$class_assignments_yaml
  cay_path <- if (substr(cay_raw, 1, 1) == "/") cay_raw else
              normalizePath(file.path(dirname(args$metadata), cay_raw), mustWork = TRUE)
  message("loading class assignments YAML: ", cay_path)
  cay <- yaml::read_yaml(cay_path)
  stopifnot("class_assignments_yaml missing 'assignments' field" = !is.null(cay$assignments))
  col_name <- contrast_spec$factor_col
  ass_vec <- unlist(cay$assignments, use.names = TRUE)
  covars[[col_name]] <- unname(ass_vec[covars$library_id])
  n_missing <- sum(is.na(covars[[col_name]]))
  message("  joined ", sum(!is.na(covars[[col_name]])),
          "/", nrow(covars), " library assignments into covars$", col_name,
          " (", n_missing, " missing)")
  if (!is.null(cay$metadata)) {
    message("  source method: ", cay$metadata$method %||% "unspecified",
            "  cut: ", cay$metadata$cut %||% "n/a")
  }
}

# Build factor column from rule
covars$contrast_factor <- build_factor_col(covars, contrast_spec$factor_col)
# Also expose under its YAML name so libraries_filter can reference it
covars[[contrast_spec$factor_col]] <- covars$contrast_factor

# Normalize YAML levels: only lowercase if they're YAML booleans (parsed by R as
# logical → "TRUE"/"FALSE" via as.character). Otherwise preserve case so
# libraries_filter expressions still match (Apr 27 fix in joint pipeline).
levels_chr <- as.character(contrast_spec$levels)
if (all(levels_chr %in% c("TRUE", "FALSE"))) {
  levels_chr <- tolower(levels_chr)
}
factor_vals <- as.character(covars$contrast_factor)

# Restrict to the levels listed (drops NA / extraneous values)
covars <- covars[factor_vals %in% levels_chr, , drop = FALSE]
covars$contrast_factor <- factor(as.character(covars$contrast_factor),
                                 levels = levels_chr)
covars[[contrast_spec$factor_col]] <- covars$contrast_factor

# Apply libraries_filter if any
covars <- apply_libraries_filter(covars, contrast_spec$libraries_filter)
message("after contrast filter: ", nrow(covars), " libraries")
print(table(covars$contrast_factor, useNA = "ifany"))

covars$patient_id <- factor(covars$patient_id)
rownames(covars)  <- covars$library_id

# Restrict bridge to scope (compartment subset for non-flex_full scopes).
# Track D uses flex_* prefix instead of joint_*.
if (args$scope != "flex_full") {
  comp_name <- sub("^flex_", "", args$scope)
  comp_match <- list(Epithelial = "Epithelial", Immune = "Immune",
                     Stromal = "Stromal")[[comp_name]]
  if (is.null(comp_match)) stop("Cannot infer compartment from scope: ", args$scope)
  if ("compartment" %in% colnames(bridge)) {
    bridge <- bridge[bridge$compartment == comp_match, ]
    message("scope=", args$scope, " -> ", nrow(bridge), " cells in compartment ", comp_match)
  } else {
    message("scope=", args$scope, " (bridge has no compartment column; relying on Milo cell intersection)")
  }
}

# --- load Milo ---
message("loading milo...")
milo <- readRDS(milo_path)
message("  cells: ", ncol(milo), " | nhoods: ", ncol(nhoods(milo)))

# Spatial_HBCA_DA Milos were built with sample=patient_id (n=4); DA needs library_id
# replication. Detect and recount per library_id derived from cell_id.
# Cell-id format: {library_id}_{barcode-1} (e.g. Pat1_P1_AAAC...AACGGGAA-1).
if (ncol(nhoodCounts(milo)) < 10) {
  message("recountCells: nhoodCounts has ", ncol(nhoodCounts(milo)),
          " columns (patient-level). Recounting per library_id...")
  cd <- as.data.frame(colData(milo))
  cell_ids <- if (!is.null(cd$cell_id)) cd$cell_id else rownames(cd)
  cd$library_id <- sub("_[ACGTN]+-[0-9]+$", "", cell_ids)
  n_libs <- length(unique(cd$library_id))
  message("  derived ", n_libs, " unique library_ids from cell_ids")
  colData(milo)$library_id <- cd$library_id
  milo <- countCells(milo, meta.data = cd, sample = "library_id")
  message("  nhoodCounts now has ", ncol(nhoodCounts(milo)), " columns")
}

# --- align to nhoodCounts ---
common_libs <- intersect(colnames(nhoodCounts(milo)), rownames(covars))
message("libraries: nhoodCounts=", ncol(nhoodCounts(milo)),
        " covars=", nrow(covars), " intersect=", length(common_libs))
stopifnot("No overlapping libraries with nhoodCounts" = length(common_libs) >= 4)

milo_sub <- milo
nhoodCounts(milo_sub) <- nhoodCounts(milo)[, common_libs, drop = FALSE]
design_df <- covars[common_libs, , drop = FALSE]

# --- build design matrix ---
term_remap <- contrast_spec$factor_col
formula_str_resolved <- gsub(paste0("\\b", term_remap, "\\b"), "contrast_factor", formula_str)
message("design formula (after factor remap): ", formula_str_resolved)
design <- model.matrix(as.formula(formula_str_resolved), data = design_df)
message("design matrix: ", nrow(design), " x ", ncol(design))

# --- run testNhoods ---
message("running testNhoods (fdr.weighting=", fdr_w, ", reduced.dim=", reduced_dim, ")...")
da_res <- testNhoods(
  milo_sub,
  design        = design,
  design.df     = design_df,
  reduced.dim   = reduced_dim,
  fdr.weighting = fdr_w
)

# --- annotate with majority L1.5 (label_set parameterized via YAML / --label-column) ---
cell_labels <- setNames(bridge[[label_column]], bridge$cell_id)
ann <- annotate_nhoods(milo, cell_labels)
da_res$nhood_label_l1p5      <- ann$label
da_res$nhood_label_l1p5_frac <- ann$frac

# --- write outputs ---
csv_path <- file.path(args$out_dir, "da_results.csv")
write.csv(da_res, csv_path, row.names = FALSE)
sig01  <- sum(da_res$SpatialFDR < 0.1,  na.rm = TRUE)
sig005 <- sum(da_res$SpatialFDR < 0.05, na.rm = TRUE)
message("significant nhoods: FDR<0.1 = ", sig01, " | FDR<0.05 = ", sig005, " / ", nrow(da_res))

run_log <- list(
  scope            = args$scope,
  contrast         = args$contrast,
  tier             = args$tier,
  formula_yaml     = formula_str,
  formula_resolved = formula_str_resolved,
  contrast_levels  = as.character(contrast_spec$levels),
  positive_level   = as.character(contrast_spec$levels)[2],
  n_libraries      = length(common_libs),
  n_per_group      = as.list(table(design_df$contrast_factor)),
  n_nhoods         = nrow(da_res),
  n_sig_FDR_0.1    = sig01,
  n_sig_FDR_0.05   = sig005,
  fdr_weighting    = fdr_w,
  reduced_dim      = reduced_dim,
  label_set        = label_set,
  label_column     = label_column,
  design_dim       = paste0(nrow(design), " x ", ncol(design)),
  date             = as.character(Sys.time()),
  milo             = milo_path
)
yaml::write_yaml(run_log, file.path(args$out_dir, "run_log.yaml"))
message("wrote ", csv_path)
message("=== done ===")
