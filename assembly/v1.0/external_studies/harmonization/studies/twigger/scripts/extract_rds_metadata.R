#!/usr/bin/env Rscript
# ==============================================================================
# TWIGGER RDS EXTRACTION - UNSUPERVISED ANALYSIS
# ==============================================================================
# Purpose: Extract and analyze Twigger RDS file to resolve open questions
# Priority: HIGH - Critical for sample->cell mapping
#
# Open Questions to Resolve:
# 1. sample_to_cell_mapping: How do sample names (NMC1-7, LMC1-9) map to cells?
# 2. cell_count_discrepancy: Why do only 93,395/110,744 cells match after normalization?
# 3. batch_code_translation: Mapping between metadata Batch (1,2,3) and cell ID batches (HMC, RB1-5)
# 4. cell_id_unsupervised_analysis: Actual cell ID format
# 5. rds_metadata_columns: All @meta.data columns
# 6. lactation_status_in_rds: Is lactation_status available per-cell?
# ==============================================================================

library(Seurat)
library(yaml)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

RDS_PATH <- "${SOURCE_COMPONENT_STUDIES}/twigger_l1_fastmnn.rds"
OUTPUT_PATH <- "${SOURCE_IHBCAV1_HARMONIZATION%/}/studies/twigger/extracted/raw_data_extraction.yaml"

cat("==============================================================================\n")
cat("TWIGGER RDS EXTRACTION - UNSUPERVISED ANALYSIS\n")
cat("==============================================================================\n")
cat("Input:  ", RDS_PATH, "\n")
cat("Output: ", OUTPUT_PATH, "\n")
cat("==============================================================================\n\n")

# ==============================================================================
# LOAD RDS FILE
# ==============================================================================

cat("Loading RDS file...\n")
tryCatch({
  obj <- readRDS(RDS_PATH)
  cat("SUCCESS: RDS file loaded\n\n")
}, error = function(e) {
  cat("ERROR: Failed to load RDS file\n")
  cat("Error message:", e$message, "\n")
  quit(status = 1)
})

# Determine object type
obj_class <- class(obj)[1]
cat("Object class:", obj_class, "\n")

# Handle potential Seurat version issues
if (inherits(obj, "Seurat")) {
  tryCatch({
    obj <- UpdateSeuratObject(obj)
    cat("Seurat object updated successfully\n")
  }, error = function(e) {
    cat("Warning: UpdateSeuratObject failed:", e$message, "\n")
    cat("Proceeding with original object\n")
  })
}

# ==============================================================================
# EXTRACT CELL IDS
# ==============================================================================

cat("\n==============================================================================\n")
cat("EXTRACTING CELL IDS\n")
cat("==============================================================================\n")

cell_ids <- colnames(obj)
total_cells <- length(cell_ids)
cat("Total cells:", total_cells, "\n")
cat("Expected cells: 110,744\n")
cat("Match:", ifelse(total_cells == 110744, "YES", "NO"), "\n\n")

# ==============================================================================
# UNSUPERVISED CELL ID PATTERN ANALYSIS
# ==============================================================================

cat("==============================================================================\n")
cat("UNSUPERVISED CELL ID PATTERN ANALYSIS\n")
cat("==============================================================================\n")

# Store results
pattern_analysis <- list()

# Sample cell IDs for display
cat("\nSample cell IDs (first 20):\n")
head_ids <- head(cell_ids, 20)
for (i in seq_along(head_ids)) {
  cat(sprintf("  %2d: %s\n", i, head_ids[i]))
}

# Analyze by potential separators
analyze_by_separator <- function(ids, sep, sep_name) {
  split_results <- strsplit(ids, sep, fixed = TRUE)
  n_components <- sapply(split_results, length)

  list(
    separator = sep_name,
    component_counts = table(n_components),
    examples = lapply(1:min(5, length(split_results)), function(i) {
      list(cell_id = ids[i], components = split_results[[i]])
    })
  )
}

# Test underscore separator
cat("\n--- Separator: underscore (_) ---\n")
underscore_analysis <- analyze_by_separator(cell_ids, "_", "underscore")
print(underscore_analysis$component_counts)

# Test hyphen separator
cat("\n--- Separator: hyphen (-) ---\n")
hyphen_analysis <- analyze_by_separator(cell_ids, "-", "hyphen")
print(hyphen_analysis$component_counts)

# ==============================================================================
# BARCODE IDENTIFICATION (16bp ACGT)
# ==============================================================================

cat("\n==============================================================================\n")
cat("BARCODE IDENTIFICATION\n")
cat("==============================================================================\n")

# 10x barcodes are 16bp ACGT sequences
barcode_pattern <- "^[ACGT]{16}$"
barcode_with_suffix <- "^[ACGT]{16}-[0-9]+$"

# Check if cell IDs contain barcodes
# Try extracting from different positions

identify_barcode_position <- function(ids) {
  # Split by underscore first
  split_ids <- strsplit(ids, "_", fixed = TRUE)
  n_parts <- sapply(split_ids, length)

  results <- list()

  for (n in unique(n_parts)) {
    subset_idx <- which(n_parts == n)
    if (length(subset_idx) == 0) next

    subset_split <- split_ids[subset_idx]

    for (pos in 1:n) {
      parts <- sapply(subset_split, function(x) if (length(x) >= pos) x[pos] else NA)
      parts <- parts[!is.na(parts)]

      # Check for barcode pattern (strip suffix first)
      stripped <- sub("-[0-9]+$", "", parts)
      is_barcode <- grepl(barcode_pattern, stripped)

      results[[paste0("parts_", n, "_pos_", pos)]] <- list(
        n_parts = n,
        position = pos,
        n_cells = length(parts),
        n_barcodes = sum(is_barcode),
        pct_barcodes = round(100 * sum(is_barcode) / length(parts), 1),
        sample_values = head(unique(parts), 10)
      )
    }
  }

  results
}

barcode_positions <- identify_barcode_position(cell_ids)

cat("\nBarcode position analysis (by underscore split):\n")
for (name in names(barcode_positions)) {
  res <- barcode_positions[[name]]
  if (res$pct_barcodes > 50) {
    cat(sprintf("  %s: %d cells, %d barcodes (%.1f%%)\n",
                name, res$n_cells, res$n_barcodes, res$pct_barcodes))
  }
}

# ==============================================================================
# BATCH CODE IDENTIFICATION - CRITICAL
# ==============================================================================

cat("\n==============================================================================\n")
cat("BATCH CODE IDENTIFICATION (CRITICAL)\n")
cat("==============================================================================\n")

# Known batch codes from edge_cases.yaml: HMC, RB1, RB2, RB3, RB4, RB5
known_batches <- c("HMC", "RB1", "RB2", "RB3", "RB4", "RB5")

# Check for batch codes in cell IDs
batch_analysis <- list()

for (batch in known_batches) {
  # Count cells containing this batch code
  contains_batch <- grepl(batch, cell_ids, fixed = TRUE)

  # Check position of batch code
  # Standard format: BATCH_BARCODE-1 -> batch at start
  # Reversed format: BARCODE-BATCH -> batch at end after hyphen

  starts_with <- grepl(paste0("^", batch), cell_ids)
  ends_with <- grepl(paste0("-", batch, "$"), cell_ids)

  batch_analysis[[batch]] <- list(
    total_containing = sum(contains_batch),
    at_start = sum(starts_with),
    at_end_after_hyphen = sum(ends_with),
    other_position = sum(contains_batch) - sum(starts_with) - sum(ends_with)
  )
}

cat("\nBatch code distribution in cell IDs:\n")
cat(sprintf("%-8s %8s %10s %15s %8s\n",
            "Batch", "Total", "At Start", "End (reversed)", "Other"))
cat(paste(rep("-", 55), collapse = ""), "\n")

for (batch in known_batches) {
  ba <- batch_analysis[[batch]]
  cat(sprintf("%-8s %8d %10d %15d %8d\n",
              batch, ba$total_containing, ba$at_start,
              ba$at_end_after_hyphen, ba$other_position))
}

# Total by format
standard_format_cells <- sum(sapply(batch_analysis, function(x) x$at_start))
reversed_format_cells <- sum(sapply(batch_analysis, function(x) x$at_end_after_hyphen))

cat("\n--- Format Summary ---\n")
cat("Standard format (BATCH_BARCODE-1):", standard_format_cells, "cells\n")
cat("Reversed format (BARCODE-BATCH):", reversed_format_cells, "cells\n")
cat("Total with known batch:", standard_format_cells + reversed_format_cells, "cells\n")
cat("Cells without known batch:", total_cells - standard_format_cells - reversed_format_cells, "cells\n")

# ==============================================================================
# EXTRACT UNIQUE PREFIXES/PATTERNS
# ==============================================================================

cat("\n==============================================================================\n")
cat("UNIQUE PREFIX ANALYSIS\n")
cat("==============================================================================\n")

# Extract first component (before underscore)
first_component <- sapply(strsplit(cell_ids, "_", fixed = TRUE), `[`, 1)
first_unique <- sort(unique(first_component))

cat("Unique first components (before underscore):", length(first_unique), "\n")
cat("\nDistribution:\n")
first_table <- sort(table(first_component), decreasing = TRUE)
for (i in seq_len(min(20, length(first_table)))) {
  cat(sprintf("  %-20s: %7d cells\n", names(first_table)[i], first_table[i]))
}
if (length(first_table) > 20) {
  cat(sprintf("  ... and %d more unique values\n", length(first_table) - 20))
}

# ==============================================================================
# METADATA EXTRACTION - ALL COLUMNS
# ==============================================================================

cat("\n==============================================================================\n")
cat("METADATA EXTRACTION (@meta.data)\n")
cat("==============================================================================\n")

metadata <- obj@meta.data
meta_cols <- colnames(metadata)

cat("Number of metadata columns:", length(meta_cols), "\n")
cat("\nColumn names:\n")
for (col in meta_cols) {
  cat("  -", col, "\n")
}

# ==============================================================================
# METADATA COLUMN ANALYSIS - EXHAUSTIVE
# ==============================================================================

cat("\n==============================================================================\n")
cat("METADATA COLUMN ANALYSIS - ALL UNIQUE VALUES\n")
cat("==============================================================================\n")

metadata_summary <- list()

for (col in meta_cols) {
  values <- metadata[[col]]
  n_unique <- length(unique(values))
  n_na <- sum(is.na(values))

  cat("\n--- Column:", col, "---\n")
  cat("  Type:", class(values)[1], "\n")
  cat("  Unique values:", n_unique, "\n")
  cat("  NA count:", n_na, "\n")

  if (n_unique <= 50) {
    # Show all values with counts
    val_table <- sort(table(values, useNA = "ifany"), decreasing = TRUE)
    cat("  Value distribution:\n")
    for (i in seq_along(val_table)) {
      val_name <- names(val_table)[i]
      if (is.na(val_name)) val_name <- "<NA>"
      cat(sprintf("    %-30s: %7d (%.1f%%)\n",
                  val_name, val_table[i], 100 * val_table[i] / total_cells))
    }
  } else {
    cat("  (Too many unique values to display all)\n")
    cat("  Sample values:", paste(head(unique(as.character(values)), 10), collapse = ", "), "...\n")
    if (is.numeric(values)) {
      cat("  Range:", min(values, na.rm = TRUE), "-", max(values, na.rm = TRUE), "\n")
      cat("  Mean:", round(mean(values, na.rm = TRUE), 2), "\n")
    }
  }

  # Store for YAML output
  if (n_unique <= 50) {
    val_table <- table(values, useNA = "ifany")
    metadata_summary[[col]] <- list(
      type = class(values)[1],
      n_unique = n_unique,
      n_na = n_na,
      values = as.list(val_table)
    )
  } else {
    metadata_summary[[col]] <- list(
      type = class(values)[1],
      n_unique = n_unique,
      n_na = n_na,
      sample_values = head(unique(as.character(values)), 20)
    )
  }
}

# ==============================================================================
# CRITICAL COLUMN SEARCH
# ==============================================================================

cat("\n==============================================================================\n")
cat("CRITICAL COLUMN SEARCH\n")
cat("==============================================================================\n")

# Look for sample/donor columns
sample_patterns <- c("sample", "donor", "orig.ident", "patient", "subject", "individual")
lactation_patterns <- c("lactat", "milk", "breast", "status", "condition")

cat("\n--- Sample/Donor related columns ---\n")
for (col in meta_cols) {
  for (pat in sample_patterns) {
    if (grepl(pat, col, ignore.case = TRUE)) {
      cat("  FOUND:", col, "\n")
      cat("    Unique values:", length(unique(metadata[[col]])), "\n")
      if (length(unique(metadata[[col]])) <= 20) {
        cat("    Values:", paste(unique(metadata[[col]]), collapse = ", "), "\n")
      }
    }
  }
}

cat("\n--- Lactation related columns ---\n")
for (col in meta_cols) {
  for (pat in lactation_patterns) {
    if (grepl(pat, col, ignore.case = TRUE)) {
      cat("  FOUND:", col, "\n")
      cat("    Unique values:", length(unique(metadata[[col]])), "\n")
      if (length(unique(metadata[[col]])) <= 20) {
        cat("    Values:", paste(unique(metadata[[col]]), collapse = ", "), "\n")
      }
    }
  }
}

# ==============================================================================
# CHECK FOR REDUCTIONS/EMBEDDINGS
# ==============================================================================

cat("\n==============================================================================\n")
cat("AVAILABLE REDUCTIONS/EMBEDDINGS\n")
cat("==============================================================================\n")

if (inherits(obj, "Seurat")) {
  reductions <- names(obj@reductions)
  cat("Available reductions:", length(reductions), "\n")
  for (red in reductions) {
    dims <- ncol(obj@reductions[[red]]@cell.embeddings)
    cat(sprintf("  - %s: %d dimensions\n", red, dims))
  }
}

# ==============================================================================
# COMPILE OUTPUT YAML
# ==============================================================================

cat("\n==============================================================================\n")
cat("COMPILING OUTPUT YAML\n")
cat("==============================================================================\n")

# Prepare batch analysis for YAML (convert to simple lists)
batch_yaml <- lapply(batch_analysis, function(x) {
  list(
    total_containing = as.integer(x$total_containing),
    at_start_standard = as.integer(x$at_start),
    at_end_reversed = as.integer(x$at_end_after_hyphen)
  )
})

# Prepare first component distribution for YAML
first_dist <- as.list(head(first_table, 50))

output <- list(
  extraction_info = list(
    source_file = RDS_PATH,
    extraction_date = as.character(Sys.time()),
    object_class = obj_class
  ),

  cell_counts = list(
    total_cells = as.integer(total_cells),
    expected_cells = 110744L,
    count_matches = total_cells == 110744
  ),

  cell_id_patterns = list(
    sample_ids = as.list(head(cell_ids, 20)),
    separator_analysis = list(
      underscore = list(
        component_counts = as.list(underscore_analysis$component_counts)
      ),
      hyphen = list(
        component_counts = as.list(hyphen_analysis$component_counts)
      )
    ),
    first_component_unique_count = length(first_unique),
    first_component_distribution = first_dist
  ),

  batch_code_analysis = list(
    known_batches = known_batches,
    per_batch = batch_yaml,
    summary = list(
      standard_format_cells = as.integer(standard_format_cells),
      reversed_format_cells = as.integer(reversed_format_cells),
      total_with_known_batch = as.integer(standard_format_cells + reversed_format_cells),
      cells_without_known_batch = as.integer(total_cells - standard_format_cells - reversed_format_cells)
    )
  ),

  metadata_columns = list(
    column_names = meta_cols,
    column_count = length(meta_cols),
    column_details = metadata_summary
  ),

  critical_findings = list(
    sample_donor_columns = meta_cols[sapply(meta_cols, function(col) {
      any(sapply(sample_patterns, function(pat) grepl(pat, col, ignore.case = TRUE)))
    })],
    lactation_columns = meta_cols[sapply(meta_cols, function(col) {
      any(sapply(lactation_patterns, function(pat) grepl(pat, col, ignore.case = TRUE)))
    })]
  )
)

# Add reductions if Seurat
if (inherits(obj, "Seurat")) {
  output$reductions <- lapply(names(obj@reductions), function(red) {
    list(name = red, dimensions = ncol(obj@reductions[[red]]@cell.embeddings))
  })
  names(output$reductions) <- names(obj@reductions)
}

# ==============================================================================
# WRITE YAML OUTPUT
# ==============================================================================

cat("Writing output to:", OUTPUT_PATH, "\n")

# Ensure output directory exists
output_dir <- dirname(OUTPUT_PATH)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Write YAML
yaml_output <- yaml::as.yaml(output, indent = 2, indent.mapping.sequence = TRUE)

# Add header
header <- "# =============================================================================
# TWIGGER RDS EXTRACTION - RAW DATA
# =============================================================================
# Generated by: extract_rds_metadata.R
# Source: twigger_l1_fastmnn.rds
# Purpose: Answer open questions from edge_cases.yaml
# =============================================================================

"

writeLines(paste0(header, yaml_output), OUTPUT_PATH)

cat("SUCCESS: Output written to", OUTPUT_PATH, "\n")

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n==============================================================================\n")
cat("EXTRACTION SUMMARY\n")
cat("==============================================================================\n")
cat("Total cells:", total_cells, "\n")
cat("Cell count matches expected:", total_cells == 110744, "\n")
cat("Metadata columns:", length(meta_cols), "\n")
cat("Standard format cells:", standard_format_cells, "\n")
cat("Reversed format cells:", reversed_format_cells, "\n")
cat("Sample/donor columns found:", length(output$critical_findings$sample_donor_columns), "\n")
cat("Lactation columns found:", length(output$critical_findings$lactation_columns), "\n")
cat("\nExtraction complete.\n")
