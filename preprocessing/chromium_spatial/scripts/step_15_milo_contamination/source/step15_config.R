# =============================================================================
# Step 15 — Configuration
# =============================================================================
# PURPOSE:
# Single source of truth for all paths, column definitions, and constants.
# NO FALLBACKS — fail early with clear errors if required inputs missing.
#
# NOTE (canonical port): paths here describe the dev-tree location of the
# step 14 outputs. For canonical runs the wrapper resolves inputs via the
# CFG_PREPROCESSING_STEP_14 / CFG_PREPROCESSING_STEP_15 env vars before
# calling into these helpers; see scripts/step_15_milo_contamination/wrapper.R.
# =============================================================================

# =============================================================================
# INPUT PATHS (from Step 14)
# =============================================================================

#' Get Step 14 base path
#' @param project_root Project root directory
#' @return Path to Step 14 metadata directory
get_step14_base <- function(project_root) {
    canonical <- Sys.getenv("CFG_PREPROCESSING_STEP_14", "")
    if (nzchar(canonical)) return(file.path(canonical, "metadata"))
    file.path(project_root, "project/01_Preprocessing/outputs/14_CompartmentClassification/metadata")
}

#' Get all Step 14 input paths
get_step14_paths <- function(project_root) {
    base <- get_step14_base(project_root)
    list(
        cell_compartments = file.path(base, "cell_compartments.csv"),
        sample_stats      = file.path(base, "sample_compartment_stats.csv"),
        epidermal_yaml    = file.path(base, "epidermal_samples.yaml")
    )
}

# =============================================================================
# OUTPUT PATHS (Step 15)
# =============================================================================

#' Get Step 15 output paths
get_step15_paths <- function(project_root) {
    canonical <- Sys.getenv("CFG_PREPROCESSING_STEP_15", "")
    base <- if (nzchar(canonical)) canonical else
        file.path(project_root, "project/01_Preprocessing/outputs/15_MiloContamination")
    list(
        base         = base,
        milo_object  = file.path(base, "intermediate/milo_object.rds"),
        metadata_dir = file.path(base, "metadata"),
        viz_dir      = file.path(base, "visualizations"),
        da_results   = file.path(base, "metadata/da_testing_results.csv"),
        cell_flags   = file.path(base, "metadata/contamination_specific_cells.csv"),
        summary_yaml = file.path(base, "metadata/contamination_detection_summary.yaml")
    )
}

# =============================================================================
# COLUMN DEFINITIONS
# =============================================================================
# Explicit column names from Step 14 outputs. NO SEARCHING.

COLS <- list(
    cell_id     = "cell_id",
    sample_id   = "sample_id",
    patient_id  = "patient_id",

    position    = "position",     # Milo colData
    position_id = "position_id",  # Step 14 output column (mapped to position)

    cluster     = "leiden_scvi_1.0",
    celltype    = "cluster_annotation",
    compartment = "compartment",

    has_epidermal = "has_epidermal"
)

# =============================================================================
# CONSTANTS
# =============================================================================

# Contamination threshold (from Step 14)
EPIDERMAL_THRESHOLD_PCT <- 5.0

# P1 (areola) — no clean counterpart in the sampling design
P1_POSITIONS <- c("P1")

# =============================================================================
# VALIDATION
# =============================================================================

validate_step14_inputs <- function(project_root) {
    step14_paths <- get_step14_paths(project_root)
    step15_paths <- get_step15_paths(project_root)

    required_files <- c(
        step14_paths$cell_compartments,
        step14_paths$sample_stats,
        step14_paths$epidermal_yaml,
        step15_paths$milo_object
    )

    missing <- required_files[!file.exists(required_files)]

    if (length(missing) > 0) {
        stop("REQUIRED inputs not found:\n  ",
             paste(missing, collapse = "\n  "),
             "\n\nRun Step 14 (Compartment Classification) and Stage 1 (Milo Object) first.")
    }

    cat("=== INPUT VALIDATION PASSED ===\n")
    cat("All required Step 14 inputs found.\n")
    invisible(TRUE)
}

validate_milo_coldata <- function(milo_obj) {
    cell_meta <- as.data.frame(SummarizedExperiment::colData(milo_obj))

    required_cols <- c(COLS$sample_id, COLS$patient_id, COLS$position)
    missing <- setdiff(required_cols, colnames(cell_meta))

    if (length(missing) > 0) {
        stop("Required columns missing from Milo colData: ",
             paste(missing, collapse = ", "),
             "\n\nAvailable columns: ", paste(colnames(cell_meta), collapse = ", "))
    }

    cat("Milo colData validation passed.\n")
    cat("  Required columns found: ", paste(required_cols, collapse = ", "), "\n")
    invisible(TRUE)
}

# =============================================================================
# DATA LOADING HELPERS
# =============================================================================

load_contaminated_samples <- function(yaml_path) {
    if (!file.exists(yaml_path)) stop("Epidermal samples YAML not found: ", yaml_path)

    step14_output <- yaml::read_yaml(yaml_path)
    contaminated_ids <- sapply(step14_output$epidermal_samples, function(x) x$sample_id)

    cat("Loaded ", length(contaminated_ids), " contaminated samples from Step 14\n", sep = "")
    cat("  Threshold: ", step14_output$threshold_pct, "%\n", sep = "")

    contaminated_ids
}

load_sample_stats <- function(csv_path) {
    if (!file.exists(csv_path)) stop("Sample stats CSV not found: ", csv_path)

    sample_stats <- read.csv(csv_path, stringsAsFactors = FALSE)

    # Map position_id -> position (Step 14 uses position_id, Milo uses position)
    if (COLS$position_id %in% colnames(sample_stats) &&
        !COLS$position %in% colnames(sample_stats)) {
        sample_stats[[COLS$position]] <- sample_stats[[COLS$position_id]]
        cat("Mapped position_id -> position\n")
    }

    cat("Loaded sample stats: ", nrow(sample_stats), " samples\n", sep = "")
    sample_stats
}

load_cell_annotations <- function(csv_path) {
    if (!file.exists(csv_path)) stop("Cell compartments CSV not found: ", csv_path)

    annotations <- read.csv(csv_path, stringsAsFactors = FALSE)

    cat("Loaded cell annotations: ", nrow(annotations), " cells\n", sep = "")
    cat("  Columns: ", paste(colnames(annotations), collapse = ", "), "\n", sep = "")

    annotations
}

get_contaminated_positions <- function(sample_stats) {
    contaminated_df <- sample_stats[sample_stats[[COLS$has_epidermal]] == TRUE, ]

    positions_by_patient <- split(
        contaminated_df[[COLS$position]],
        contaminated_df[[COLS$patient_id]]
    )
    positions_by_patient <- lapply(positions_by_patient, unique)

    cat("Contaminated positions by patient:\n")
    for (pat in names(positions_by_patient)) {
        cat("  ", pat, ": ", paste(positions_by_patient[[pat]], collapse = ", "), "\n", sep = "")
    }

    positions_by_patient
}
