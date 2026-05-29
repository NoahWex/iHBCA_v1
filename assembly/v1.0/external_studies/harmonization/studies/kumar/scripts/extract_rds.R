# =============================================================================
# KUMAR TRACK 2: UNSUPERVISED RDS EXTRACTION
# =============================================================================
# Extracts ALL data from the ORIGINAL Kumar RDS file.
# Philosophy: Extract everything. Interpret nothing. Let data speak.
#
# Source: ${SOURCE_KUMAR_ORIGINAL%/}/kumar.rds
# Output: studies/kumar/extracted/raw_data_extraction.yaml
#
# Open Questions to Answer:
# 1. missing_samples_explanation - What explains the c01-c65 gap?
# 2. cell_id_unsupervised_analysis - What is the actual cell ID format?
# 3. metadata_column_mapping - How do supplemental columns map to RDS?
# 4. ethnicity_as_condition - Is ethnicity present?
# 5. bmi_as_condition - Is BMI present?
# 6. cell_distribution - How do 714,331 cells distribute across samples?
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(yaml)
})

# =============================================================================
# CONFIGURATION
# =============================================================================

# Source: ORIGINAL Kumar RDS (not iHBCA converted)
KUMAR_RDS_PATH <- "${SOURCE_COMPONENT_STUDIES}/kumar.rds"

# Output
BASE_DIR <- "${SOURCE_IHBCAV1_HARMONIZATION%/}"
OUTPUT_DIR <- file.path(BASE_DIR, "studies", "kumar", "extracted")
OUTPUT_FILE <- file.path(OUTPUT_DIR, "raw_data_extraction.yaml")

# =============================================================================
# UNSUPERVISED CELL ID ANALYSIS
# =============================================================================

#' Perform UNSUPERVISED analysis of cell ID patterns
#' Discovers format without assumptions
analyze_cell_ids_unsupervised <- function(cell_ids) {
  n_total <- length(cell_ids)
  message(sprintf("[%s] Analyzing %d cell IDs...", Sys.time(), n_total))

  # =========================================================================
  # 1. RAW SAMPLES
  # =========================================================================
  samples <- list(
    first_20 = head(cell_ids, 20),
    last_20 = tail(cell_ids, 20),
    random_20 = cell_ids[sample(seq_along(cell_ids), min(20, n_total))]
  )

  # =========================================================================
  # 2. SEPARATOR DETECTION
  # =========================================================================
  # Count occurrences of potential separators in first 1000 cells
  test_sample <- cell_ids[1:min(1000, n_total)]

  separators <- c("_", "-", ":", ".")
  sep_analysis <- lapply(separators, function(sep) {
    counts <- sapply(test_sample, function(id) {
      length(gregexpr(sep, id, fixed = TRUE)[[1]])
    })
    counts[counts < 0] <- 0  # gregexpr returns -1 if no match
    list(
      separator = sep,
      mean_count = mean(counts),
      max_count = max(counts),
      min_count = min(counts),
      consistent = length(unique(counts)) == 1
    )
  })
  names(sep_analysis) <- separators

  # Determine primary separator (most consistent with count > 0)
  primary_sep <- "_"  # default
  for (sep in separators) {
    if (sep_analysis[[sep]]$consistent && sep_analysis[[sep]]$mean_count > 0) {
      if (sep_analysis[[sep]]$mean_count > sep_analysis[[primary_sep]]$mean_count) {
        primary_sep <- sep
      }
    }
  }

  # =========================================================================
  # 3. COMPONENT ANALYSIS (using primary separator)
  # =========================================================================
  splits <- strsplit(cell_ids, primary_sep, fixed = TRUE)
  n_components <- sapply(splits, length)

  # Count distinct patterns
  component_count_dist <- as.list(table(n_components))

  # Get most common pattern count
  most_common_n <- as.integer(names(which.max(table(n_components))))

  # Analyze each component position
  component_details <- list()
  for (i in seq_len(max(n_components))) {
    # Extract component i from all cells that have it
    comp_values <- sapply(splits[n_components >= i], `[`, i)

    if (length(comp_values) == 0) next

    # Unique values
    unique_vals <- unique(comp_values)
    n_unique <- length(unique_vals)

    # Character analysis
    char_lengths <- nchar(comp_values)

    # Detect if barcode (16bp ACGT only)
    is_acgt_only <- all(grepl("^[ACGT]+$", comp_values))
    barcode_rate_16bp <- mean(nchar(comp_values) == 16 & grepl("^[ACGT]+$", comp_values))

    # Get value distribution (full for small, sample for large)
    if (n_unique <= 200) {
      value_counts <- as.list(table(comp_values))
    } else {
      # Sample top 50 and note truncation
      full_table <- sort(table(comp_values), decreasing = TRUE)
      value_counts <- as.list(full_table[1:min(50, length(full_table))])
      value_counts[["__NOTE__"]] <- sprintf("TRUNCATED: showing top 50 of %d unique values", n_unique)
    }

    component_details[[paste0("component_", i)]] <- list(
      position = i,
      n_unique = n_unique,
      n_cells_with_component = length(comp_values),
      char_length_min = min(char_lengths),
      char_length_max = max(char_lengths),
      char_length_mode = as.integer(names(which.max(table(char_lengths)))),
      is_acgt_only = is_acgt_only,
      barcode_rate_16bp = round(barcode_rate_16bp, 4),
      likely_type = if (barcode_rate_16bp > 0.99) {
        "barcode"
      } else if (n_unique < 20 && !is_acgt_only) {
        "categorical"
      } else if (n_unique < 500 && !is_acgt_only) {
        "sample_or_batch_id"
      } else {
        "identifier"
      },
      values = value_counts
    )
  }

  # =========================================================================
  # 4. SUFFIX DETECTION (10x -1 style)
  # =========================================================================
  suffix_pattern <- "-[0-9]+$"
  has_suffix <- grepl(suffix_pattern, cell_ids)
  suffix_rate <- mean(has_suffix)

  if (suffix_rate > 0.01) {
    suffixes <- gsub("^.*(-[0-9]+)$", "\\1", cell_ids[has_suffix])
    suffix_values <- as.list(table(suffixes))
  } else {
    suffix_values <- list()
  }

  # =========================================================================
  # 5. BARCODE EXTRACTION (last 16bp ACGT component)
  # =========================================================================
  # Try to extract barcode: typically last component, 16bp, ACGT only
  potential_barcodes <- sapply(splits, function(x) x[length(x)])
  potential_barcodes_stripped <- gsub("-[0-9]+$", "", potential_barcodes)

  valid_10x <- grepl("^[ACGT]{16}$", potential_barcodes_stripped)
  valid_10x_rate <- mean(valid_10x)

  barcode_analysis <- list(
    extraction_method = "last_component_stripped",
    valid_10x_format_rate = round(valid_10x_rate, 4),
    unique_barcodes = length(unique(potential_barcodes_stripped[valid_10x])),
    total_cells_with_valid_barcode = sum(valid_10x)
  )

  # =========================================================================
  # 6. PATTERN DISCOVERY (unique conformations)
  # =========================================================================
  # Get examples of each unique pattern
  unique_patterns <- unique(n_components)
  pattern_examples <- lapply(unique_patterns, function(n) {
    examples <- cell_ids[n_components == n]
    list(
      n_components = n,
      count = sum(n_components == n),
      examples = head(examples, 5),
      pattern_description = paste(
        sapply(seq_len(n), function(i) sprintf("{comp%d}", i)),
        collapse = primary_sep
      )
    )
  })
  names(pattern_examples) <- paste0("pattern_", unique_patterns, "_components")

  # =========================================================================
  # RETURN
  # =========================================================================
  list(
    total_cells = n_total,
    samples = samples,
    separator_analysis = sep_analysis,
    primary_separator = primary_sep,
    component_count_distribution = component_count_dist,
    most_common_component_count = most_common_n,
    component_details = component_details,
    suffix_analysis = list(
      pattern = suffix_pattern,
      has_suffix_rate = round(suffix_rate, 4),
      suffix_values = suffix_values
    ),
    barcode_analysis = barcode_analysis,
    pattern_examples = pattern_examples
  )
}

# =============================================================================
# METADATA ANALYSIS
# =============================================================================

#' Analyze all metadata columns
#' Pure discovery - extract everything
analyze_metadata_full <- function(meta) {
  n_cells <- nrow(meta)
  col_names <- names(meta)
  message(sprintf("[%s] Analyzing %d metadata columns...", Sys.time(), length(col_names)))

  columns <- list()

  for (col in col_names) {
    values <- meta[[col]]
    n_null <- sum(is.na(values))
    n_non_null <- sum(!is.na(values))

    # Type detection
    dtype <- class(values)[1]

    # For factors, get levels
    if (is.factor(values)) {
      levels_list <- levels(values)
      values <- as.character(values)  # Convert for table()
    } else {
      levels_list <- NULL
    }

    # Unique values
    unique_vals <- unique(values[!is.na(values)])
    n_unique <- length(unique_vals)

    # Value counts (ALL for small, summary for large)
    if (n_unique <= 200) {
      value_table <- table(values, useNA = "ifany")
      # Fix NA key
      names_vec <- names(value_table)
      names_vec[is.na(names_vec)] <- "NA"
      names(value_table) <- names_vec
      value_counts <- as.list(value_table)
    } else if (dtype %in% c("numeric", "integer", "double")) {
      # Numeric summary
      valid_vals <- values[!is.na(values)]
      value_counts <- list(
        "__TYPE__" = "numeric_summary",
        min = min(valid_vals),
        max = max(valid_vals),
        mean = mean(valid_vals),
        median = median(valid_vals),
        sd = sd(valid_vals),
        n_unique = n_unique
      )
    } else {
      # Too many unique values - show top 50
      full_table <- sort(table(values), decreasing = TRUE)
      value_counts <- as.list(full_table[1:min(50, length(full_table))])
      value_counts[["__NOTE__"]] <- sprintf("TRUNCATED: showing top 50 of %d unique values", n_unique)
    }

    columns[[col]] <- list(
      name = col,
      dtype = dtype,
      n_non_null = n_non_null,
      n_null = n_null,
      coverage_pct = round(n_non_null / n_cells * 100, 2),
      n_unique = n_unique,
      factor_levels = levels_list,
      values = value_counts
    )
  }

  list(
    n_columns = length(col_names),
    column_names = col_names,
    columns = columns
  )
}

# =============================================================================
# SAMPLE/DONOR DISTRIBUTION
# =============================================================================

#' Analyze cell distribution across samples/donors
analyze_cell_distribution <- function(meta, cell_ids) {
  message(sprintf("[%s] Analyzing cell distribution...", Sys.time()))

  result <- list()

  # Look for sample/donor columns
  sample_cols <- c("Sample_ID", "sample", "Sample", "orig.ident", "sample_id")
  donor_cols <- c("Patient_ID", "donor", "Donor", "patient", "patient_id", "donor_id")

  # Find sample column
  sample_col <- intersect(sample_cols, names(meta))[1]
  if (!is.na(sample_col)) {
    sample_table <- sort(table(meta[[sample_col]]), decreasing = TRUE)
    result$sample_distribution <- list(
      column_used = sample_col,
      n_samples = length(sample_table),
      sample_counts = as.list(sample_table)
    )

    # Specifically check for hbca_c## format
    sample_ids <- names(sample_table)
    hbca_pattern <- grepl("^hbca_c[0-9]+$", sample_ids)
    if (any(hbca_pattern)) {
      hbca_ids <- sample_ids[hbca_pattern]
      # Extract numeric part
      hbca_nums <- as.integer(gsub("hbca_c", "", hbca_ids))
      result$hbca_sample_analysis <- list(
        n_hbca_samples = sum(hbca_pattern),
        min_id = min(hbca_nums),
        max_id = max(hbca_nums),
        expected_range = "c01-c167 based on supplemental",
        missing_in_c01_c65 = setdiff(1:65, hbca_nums),
        present_in_c01_c65 = intersect(1:65, hbca_nums),
        present_in_c66_c167 = intersect(66:167, hbca_nums)
      )
    }
  }

  # Find donor column
  donor_col <- intersect(donor_cols, names(meta))[1]
  if (!is.na(donor_col)) {
    donor_table <- sort(table(meta[[donor_col]]), decreasing = TRUE)
    result$donor_distribution <- list(
      column_used = donor_col,
      n_donors = length(donor_table),
      donor_counts = as.list(donor_table)
    )
  }

  result
}

# =============================================================================
# OPEN QUESTIONS ANALYSIS
# =============================================================================

#' Generate analysis specifically for open questions
answer_open_questions <- function(meta, cell_id_analysis, distribution) {
  questions <- list()

  # Q1: missing_samples_explanation
  if (!is.null(distribution$hbca_sample_analysis)) {
    questions$missing_samples_explanation <- list(
      status = "ANALYZED",
      findings = distribution$hbca_sample_analysis,
      conclusion = if (length(distribution$hbca_sample_analysis$present_in_c01_c65) > 0) {
        sprintf("Found %d samples in c01-c65 range",
                length(distribution$hbca_sample_analysis$present_in_c01_c65))
      } else {
        "No samples found in c01-c65 range - confirms extraction gap"
      }
    )
  }

  # Q2: cell_id_unsupervised_analysis
  questions$cell_id_unsupervised_analysis <- list(
    status = "COMPLETE",
    primary_separator = cell_id_analysis$primary_separator,
    most_common_pattern = sprintf("%d components", cell_id_analysis$most_common_component_count),
    barcode_detection = cell_id_analysis$barcode_analysis,
    pattern_summary = names(cell_id_analysis$pattern_examples)
  )

  # Q3: metadata_column_mapping
  col_names <- names(meta)
  expected_cols <- c("Sample_ID", "Patient_ID", "Institution", "Tissue_Source",
                     "Age", "Ethnicity", "Parity", "Menopause", "BMI", "BMI_Group")
  found_mapping <- intersect(col_names, expected_cols)
  missing_from_expected <- setdiff(expected_cols, col_names)
  extra_in_rds <- setdiff(col_names, expected_cols)

  questions$metadata_column_mapping <- list(
    status = "COMPLETE",
    total_columns_in_rds = length(col_names),
    all_column_names = col_names,
    expected_from_supplemental = expected_cols,
    found_matches = found_mapping,
    missing_from_expected = missing_from_expected,
    additional_columns = extra_in_rds
  )

  # Q4: ethnicity_as_condition
  ethnicity_cols <- grep("ethnic", names(meta), ignore.case = TRUE, value = TRUE)
  if (length(ethnicity_cols) > 0) {
    eth_col <- ethnicity_cols[1]
    eth_values <- table(meta[[eth_col]])
    questions$ethnicity_as_condition <- list(
      status = "FOUND",
      column_name = eth_col,
      unique_values = names(eth_values),
      value_counts = as.list(eth_values),
      testable = length(eth_values) >= 2
    )
  } else {
    questions$ethnicity_as_condition <- list(
      status = "NOT_FOUND",
      searched_columns = names(meta),
      note = "No column containing 'ethnic' found"
    )
  }

  # Q5: bmi_as_condition
  bmi_cols <- grep("bmi", names(meta), ignore.case = TRUE, value = TRUE)
  if (length(bmi_cols) > 0) {
    bmi_results <- list(status = "FOUND", columns = list())
    for (col in bmi_cols) {
      vals <- meta[[col]]
      bmi_results$columns[[col]] <- list(
        dtype = class(vals)[1],
        n_unique = length(unique(vals[!is.na(vals)])),
        n_non_null = sum(!is.na(vals)),
        values = if (length(unique(vals)) <= 20) as.list(table(vals)) else "numeric"
      )
    }
    questions$bmi_as_condition <- bmi_results
  } else {
    questions$bmi_as_condition <- list(
      status = "NOT_FOUND",
      searched_columns = names(meta),
      note = "No column containing 'bmi' found"
    )
  }

  # Q6: cell_distribution (summary)
  questions$cell_distribution <- list(
    status = "COMPLETE",
    total_cells = nrow(meta),
    sample_count = if (!is.null(distribution$sample_distribution))
      distribution$sample_distribution$n_samples else NA,
    donor_count = if (!is.null(distribution$donor_distribution))
      distribution$donor_distribution$n_donors else NA,
    see = "cell_distribution section for full details"
  )

  questions
}

# =============================================================================
# MAIN
# =============================================================================

main <- function() {
  message("===========================================")
  message("KUMAR TRACK 2: UNSUPERVISED RDS EXTRACTION")
  message(sprintf("Start time: %s", Sys.time()))
  message("===========================================")

  # Check file exists
  if (!file.exists(KUMAR_RDS_PATH)) {
    stop(sprintf("Kumar RDS file not found: %s", KUMAR_RDS_PATH))
  }

  file_info <- file.info(KUMAR_RDS_PATH)
  message(sprintf("Source file: %s", KUMAR_RDS_PATH))
  message(sprintf("File size: %.2f GB", file_info$size / 1e9))

  # Load RDS
  message(sprintf("[%s] Loading Seurat object...", Sys.time()))
  obj <- readRDS(KUMAR_RDS_PATH)

  # Update if needed
  if (inherits(obj, "Seurat")) {
    message(sprintf("[%s] Updating Seurat object...", Sys.time()))
    obj <- tryCatch(UpdateSeuratObject(obj), error = function(e) {
      message("UpdateSeuratObject failed, continuing with original")
      obj
    })
  }

  # Extract cell IDs
  message(sprintf("[%s] Extracting cell IDs...", Sys.time()))
  cell_ids <- colnames(obj)

  # Extract metadata
  message(sprintf("[%s] Extracting metadata...", Sys.time()))
  meta <- obj@meta.data

  # Run analyses
  cell_id_analysis <- analyze_cell_ids_unsupervised(cell_ids)
  metadata_analysis <- analyze_metadata_full(meta)
  distribution_analysis <- analyze_cell_distribution(meta, cell_ids)

  # Answer open questions
  message(sprintf("[%s] Analyzing open questions...", Sys.time()))
  open_questions <- answer_open_questions(meta, cell_id_analysis, distribution_analysis)

  # Check embeddings
  message(sprintf("[%s] Checking embeddings...", Sys.time()))
  embeddings <- list()
  if (length(obj@reductions) > 0) {
    for (red_name in names(obj@reductions)) {
      red <- obj@reductions[[red_name]]
      embeddings[[red_name]] <- list(
        name = red_name,
        dims = ncol(Embeddings(red)),
        n_cells = nrow(Embeddings(red))
      )
    }
  }

  # Assemble result
  result <- list(
    study = "kumar",
    extraction_date = as.character(Sys.Date()),
    extraction_time = as.character(Sys.time()),
    source_file = list(
      path = KUMAR_RDS_PATH,
      size_bytes = file_info$size,
      size_gb = round(file_info$size / 1e9, 2),
      type = "seurat_rds",
      note = "ORIGINAL Kumar RDS (not iHBCA converted)"
    ),
    summary = list(
      total_cells = length(cell_ids),
      n_metadata_columns = ncol(meta),
      n_embeddings = length(embeddings)
    ),
    open_questions_answered = open_questions,
    cell_ids = cell_id_analysis,
    metadata = metadata_analysis,
    cell_distribution = distribution_analysis,
    embeddings = embeddings
  )

  # Write output
  message(sprintf("[%s] Writing output...", Sys.time()))
  if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
  write_yaml(result, OUTPUT_FILE)

  message(sprintf("Output written to: %s", OUTPUT_FILE))
  message("")
  message("===========================================")
  message("EXTRACTION COMPLETE")
  message(sprintf("End time: %s", Sys.time()))
  message("===========================================")

  # Print quick summary
  message("")
  message("QUICK SUMMARY:")
  message(sprintf("  Total cells: %d", result$summary$total_cells))
  message(sprintf("  Metadata columns: %d", result$summary$n_metadata_columns))
  message(sprintf("  Cell ID pattern: %d components with '%s' separator",
                  cell_id_analysis$most_common_component_count,
                  cell_id_analysis$primary_separator))

  if (!is.null(distribution_analysis$sample_distribution)) {
    message(sprintf("  Samples: %d", distribution_analysis$sample_distribution$n_samples))
  }
  if (!is.null(distribution_analysis$donor_distribution)) {
    message(sprintf("  Donors: %d", distribution_analysis$donor_distribution$n_donors))
  }
}

# Run
main()
