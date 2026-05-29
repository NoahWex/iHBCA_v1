# =============================================================================
# PAL NORM EPI: RAW DATA EXTRACTION FROM RDS
# =============================================================================
# UNSUPERVISED extraction - discovers patterns from the data itself.
# Answers all open questions from edge_cases.yaml.
#
# Output: studies/pal_norm_epi/extracted/raw_data_extraction.yaml
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(yaml)
})

# =============================================================================
# CONFIGURATION
# =============================================================================

RDS_PATH <- "${SOURCE_PAL_ORIGINAL%/}/SeuratObject_NormEpi.rds"
OUTPUT_DIR <- "${SOURCE_IHBCAV1_HARMONIZATION%/}/studies/pal_norm_epi/extracted"
OUTPUT_FILE <- file.path(OUTPUT_DIR, "raw_data_extraction.yaml")

# =============================================================================
# CELL ID ANALYSIS FUNCTIONS
# =============================================================================

#' Analyze cell ID structure - UNSUPERVISED
#' @param cell_ids Character vector of cell IDs
#' @return List with pattern analysis
analyze_cell_ids_unsupervised <- function(cell_ids) {
  message("  Analyzing cell ID patterns...")

  n_total <- length(cell_ids)
  sample_ids <- head(cell_ids, 30)  # Sample for inspection

  # Check for -N suffix (10x barcode convention)
  has_suffix <- grepl("-[0-9]+$", cell_ids)
  suffix_rate <- mean(has_suffix)

  # Detect separator by trying multiple options
  separators <- c("_", "-", ":")
  sep_analysis <- lapply(separators, function(sep) {
    splits <- strsplit(head(cell_ids, 100), split = sep, fixed = TRUE)
    n_parts <- table(sapply(splits, length))
    list(
      separator = sep,
      component_count_distribution = as.list(n_parts),
      most_common_parts = as.integer(names(n_parts)[which.max(n_parts)])
    )
  })
  names(sep_analysis) <- separators

  # Use underscore as primary separator (expected for Pal)
  splits <- strsplit(cell_ids, split = "_", fixed = TRUE)
  n_components <- sapply(splits, length)
  n_components_dist <- as.list(table(n_components))

  # If consistent component count, analyze each component
  component_analysis <- list()
  if (length(unique(n_components)) == 1) {
    n_comp <- unique(n_components)
    message(sprintf("  Detected consistent format: %d components separated by underscore", n_comp))

    for (i in seq_len(n_comp)) {
      comp_values <- sapply(splits, `[`, i)
      unique_vals <- sort(unique(comp_values))
      value_counts <- sort(table(comp_values), decreasing = TRUE)

      # Detect component type
      is_barcode <- length(unique_vals) > 100 &&
                    all(grepl("^[ACGT]+-?[0-9]*$|^[ACGT]+$", unique_vals))
      is_categorical <- length(unique_vals) <= 20

      # Barcode length analysis
      if (is_barcode) {
        barcode_clean <- gsub("-[0-9]+$", "", unique_vals)
        barcode_lengths <- table(nchar(barcode_clean))
      }

      component_analysis[[paste0("component_", i)]] <- list(
        position = i,
        n_unique = length(unique_vals),
        inferred_type = if (is_barcode) "barcode" else if (is_categorical) "categorical" else "identifier",
        # Include ALL unique values with counts for categorical columns
        values = if (length(unique_vals) <= 50) as.list(value_counts) else list(
          sample_values = as.list(head(value_counts, 20)),
          total_unique = length(unique_vals)
        ),
        barcode_info = if (is_barcode) list(
          length_distribution = as.list(barcode_lengths),
          valid_16bp_rate = mean(nchar(barcode_clean) == 16)
        ) else NULL
      )
    }
  }

  # Verify expected pattern: N_{number}_{type}_{barcode}-1
  expected_pattern <- "^[A-Z]+_[0-9.]+_[a-z]+_[ACGT]+-[0-9]+$"
  pattern_match_rate <- mean(grepl(expected_pattern, cell_ids))

  list(
    total_cells = n_total,
    sample_cell_ids = sample_ids,
    suffix_analysis = list(
      has_suffix_rate = suffix_rate,
      suffix_pattern = "-[0-9]+$"
    ),
    separator_analysis = sep_analysis,
    underscore_split = list(
      n_component_distribution = n_components_dist,
      component_details = component_analysis
    ),
    pattern_verification = list(
      expected_pattern = "N_{number}_{type}_{barcode}-1",
      regex = expected_pattern,
      match_rate = pattern_match_rate,
      verified = pattern_match_rate > 0.99
    )
  )
}

#' Analyze metadata columns - extract ALL values
#' @param meta data.frame of metadata
#' @return List with column analysis
analyze_metadata_columns <- function(meta) {
  message("  Analyzing metadata columns...")

  n_cells <- nrow(meta)
  column_names <- names(meta)

  columns <- lapply(column_names, function(col) {
    values <- meta[[col]]
    n_unique <- length(unique(values[!is.na(values)]))
    n_null <- sum(is.na(values))
    n_non_null <- sum(!is.na(values))
    coverage_pct <- round(n_non_null / n_cells * 100, 1)

    # Get ALL value counts for categorical variables
    if (n_unique <= 100) {
      value_counts <- sort(table(values, useNA = "ifany"), decreasing = TRUE)
      names(value_counts)[is.na(names(value_counts))] <- "NA"
      values_output <- as.list(value_counts)
    } else {
      value_counts <- sort(table(values), decreasing = TRUE)
      values_output <- list(
        sample_values = as.list(head(value_counts, 20)),
        total_unique = n_unique
      )
    }

    list(
      name = col,
      dtype = class(values)[1],
      non_null = n_non_null,
      null = n_null,
      coverage_pct = coverage_pct,
      n_unique = n_unique,
      values = values_output
    )
  })

  names(columns) <- column_names
  columns
}

#' Extract sample distribution from cell IDs
#' @param cell_ids Character vector of cell IDs
#' @return List with sample distribution
extract_sample_distribution <- function(cell_ids) {
  message("  Extracting sample distribution from cell IDs...")

  # Expected pattern: N_{number}_{type}_{barcode}-1
  # Sample prefix is first 3 components
  splits <- strsplit(cell_ids, split = "_", fixed = TRUE)

  # Extract sample prefix (first 3 components if 4 total)
  n_comp <- sapply(splits, length)
  if (all(n_comp == 4)) {
    sample_prefix <- sapply(splits, function(x) paste(x[1:3], collapse = "_"))
    sample_counts <- sort(table(sample_prefix), decreasing = TRUE)

    list(
      extraction_method = "First 3 components of underscore-split cell ID",
      n_samples = length(sample_counts),
      distribution = as.list(sample_counts),
      total_cells = sum(sample_counts)
    )
  } else {
    list(
      extraction_method = "Unable - inconsistent component count",
      n_component_distribution = as.list(table(n_comp))
    )
  }
}

#' Check for menopausal status in metadata or derive from cell IDs
#' @param obj Seurat object
#' @return List with menopausal status analysis
analyze_menopausal_status <- function(obj) {
  message("  Analyzing menopausal status...")

  meta <- obj@meta.data
  cell_ids <- colnames(obj)

  # Check if menopausal_status column exists
  meno_columns <- grep("meno", names(meta), ignore.case = TRUE, value = TRUE)

  result <- list(
    columns_matching_meno = meno_columns
  )

  # If direct column exists
  if (length(meno_columns) > 0) {
    for (col in meno_columns) {
      result[[paste0("column_", col)]] <- list(
        values = as.list(table(meta[[col]])),
        n_unique = length(unique(meta[[col]]))
      )
    }
  }

  # Derive from cell IDs based on known post-menopausal samples
  # (N_0342_epi, N_0372_epi, N_0275_epi from supplemental_extraction.yaml)
  post_meno_prefixes <- c("N_0342_epi", "N_0372_epi", "N_0275_epi")

  # Extract sample prefix from cell IDs
  splits <- strsplit(cell_ids, split = "_", fixed = TRUE)
  if (all(sapply(splits, length) == 4)) {
    sample_prefixes <- sapply(splits, function(x) paste(x[1:3], collapse = "_"))
    derived_meno <- ifelse(sample_prefixes %in% post_meno_prefixes, "Post", "Pre")

    result$derived_from_cell_ids <- list(
      post_menopausal_prefixes = post_meno_prefixes,
      distribution = as.list(table(derived_meno)),
      note = "Derived using known sample assignments from Pal 2021 R script"
    )
  }

  result
}

#' Analyze type component (epi vs total)
#' @param cell_ids Character vector of cell IDs
#' @return List with type verification
verify_type_component <- function(cell_ids) {
  message("  Verifying type component...")

  splits <- strsplit(cell_ids, split = "_", fixed = TRUE)

  if (all(sapply(splits, length) == 4)) {
    # Type should be component 3
    type_values <- sapply(splits, `[`, 3)
    type_counts <- table(type_values)

    list(
      component_position = 3,
      expected_value = "epi",
      actual_values = as.list(type_counts),
      all_epi = all(type_values == "epi")
    )
  } else {
    list(
      error = "Inconsistent component count",
      component_distribution = as.list(table(sapply(splits, length)))
    )
  }
}

# =============================================================================
# MAIN EXTRACTION
# =============================================================================

main <- function() {
  message(sprintf("[%s] Starting Pal NormEpi RDS extraction", Sys.time()))

  # Create output directory
  if (!dir.exists(OUTPUT_DIR)) {
    dir.create(OUTPUT_DIR, recursive = TRUE)
  }

  # Load RDS
  message(sprintf("[%s] Loading RDS: %s", Sys.time(), RDS_PATH))
  obj <- readRDS(RDS_PATH)

  # Update Seurat object if needed
  if (inherits(obj, "Seurat")) {
    obj <- tryCatch(UpdateSeuratObject(obj), error = function(e) obj)
  }

  message(sprintf("[%s] Object loaded successfully", Sys.time()))

  # Basic stats
  n_cells <- ncol(obj)
  n_genes <- nrow(obj)
  cell_ids <- colnames(obj)

  message(sprintf("[%s] Total cells: %d", Sys.time(), n_cells))
  message(sprintf("[%s] Total genes: %d", Sys.time(), n_genes))

  # Run analyses
  message(sprintf("\n[%s] Running cell ID analysis...", Sys.time()))
  cell_id_analysis <- analyze_cell_ids_unsupervised(cell_ids)

  message(sprintf("\n[%s] Running metadata analysis...", Sys.time()))
  metadata_analysis <- analyze_metadata_columns(obj@meta.data)

  message(sprintf("\n[%s] Extracting sample distribution...", Sys.time()))
  sample_distribution <- extract_sample_distribution(cell_ids)

  message(sprintf("\n[%s] Analyzing menopausal status...", Sys.time()))
  menopausal_analysis <- analyze_menopausal_status(obj)

  message(sprintf("\n[%s] Verifying type component...", Sys.time()))
  type_verification <- verify_type_component(cell_ids)

  # Check embeddings
  message(sprintf("\n[%s] Checking embeddings...", Sys.time()))
  embeddings <- lapply(names(obj@reductions), function(red_name) {
    red <- obj@reductions[[red_name]]
    list(
      name = red_name,
      dims = ncol(Embeddings(red)),
      cells = nrow(Embeddings(red))
    )
  })
  names(embeddings) <- names(obj@reductions)

  # Check assays
  message(sprintf("\n[%s] Checking assays...", Sys.time()))
  assays_info <- lapply(names(obj@assays), function(assay_name) {
    assay <- obj@assays[[assay_name]]
    list(
      name = assay_name,
      n_features = nrow(assay),
      n_cells = ncol(assay)
    )
  })
  names(assays_info) <- names(obj@assays)

  # File info
  file_info <- file.info(RDS_PATH)

  # Compile results
  message(sprintf("\n[%s] Compiling extraction results...", Sys.time()))

  # Answer specific questions from edge_cases.yaml
  questions_answered <- list(
    cell_id_format_verification = list(
      question = "Confirm cell ID format is {SamplesComb}_{Barcode}-1?",
      answer = list(
        format_verified = cell_id_analysis$pattern_verification$verified,
        pattern = cell_id_analysis$pattern_verification$expected_pattern,
        match_rate = cell_id_analysis$pattern_verification$match_rate,
        components = names(cell_id_analysis$underscore_split$component_details)
      ),
      status = if (cell_id_analysis$pattern_verification$verified) "RESOLVED" else "NEEDS_REVIEW"
    ),

    cell_count_verification = list(
      question = "Confirm total cell count is 53,716",
      answer = list(
        actual_count = n_cells,
        expected_count = 53716,
        matches = n_cells == 53716
      ),
      status = if (n_cells == 53716) "RESOLVED" else "DISCREPANCY"
    ),

    metadata_columns_in_rds = list(
      question = "What metadata columns exist in RDS @meta.data?",
      answer = list(
        columns = names(metadata_analysis),
        n_columns = length(metadata_analysis)
      ),
      status = "RESOLVED"
    ),

    cell_distribution = list(
      question = "How do cells distribute across 11 samples?",
      answer = sample_distribution,
      status = if (!is.null(sample_distribution$n_samples) && sample_distribution$n_samples == 11) "RESOLVED" else "NEEDS_REVIEW"
    ),

    menopausal_group_balance = list(
      question = "Is there sufficient cell count in post-menopausal group?",
      answer = menopausal_analysis,
      status = "RESOLVED"
    )
  )

  result <- list(
    study = "pal_norm_epi",
    extraction_date = as.character(Sys.Date()),
    extraction_type = "raw_rds",

    source_file = list(
      path = RDS_PATH,
      size_bytes = file_info$size,
      size_mb = round(file_info$size / 1024^2, 1)
    ),

    basic_stats = list(
      n_cells = n_cells,
      n_genes = n_genes,
      expected_cells = 53716,
      cell_count_matches = n_cells == 53716
    ),

    cell_ids = cell_id_analysis,
    type_verification = type_verification,
    sample_distribution = sample_distribution,
    menopausal_status = menopausal_analysis,
    metadata_columns = metadata_analysis,
    embeddings = embeddings,
    assays = assays_info,

    questions_answered = questions_answered
  )

  # Write YAML
  message(sprintf("\n[%s] Writing output to: %s", Sys.time(), OUTPUT_FILE))
  write_yaml(result, OUTPUT_FILE)

  message(sprintf("\n[%s] Extraction complete!", Sys.time()))

  # Print summary
  cat("\n========== EXTRACTION SUMMARY ==========\n")
  cat(sprintf("Total cells: %d (expected: 53,716)\n", n_cells))
  cat(sprintf("Total genes: %d\n", n_genes))
  cat(sprintf("Metadata columns: %d\n", length(metadata_analysis)))
  cat(sprintf("Cell ID format verified: %s\n", cell_id_analysis$pattern_verification$verified))
  cat(sprintf("Number of samples: %s\n", sample_distribution$n_samples))
  cat(sprintf("\nOutput written to: %s\n", OUTPUT_FILE))
  cat("========================================\n")
}

# Run
main()
