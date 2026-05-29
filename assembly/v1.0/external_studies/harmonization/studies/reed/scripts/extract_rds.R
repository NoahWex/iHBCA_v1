#!/usr/bin/env Rscript
# =============================================================================
# REED RDS EXTRACTION
# =============================================================================
# Unsupervised extraction from Reed RDS - Reed IS the iHBCA reference,
# so this defines the canonical format for all other studies.
#
# Questions to answer:
# 1. What is the exact cell ID format? (defines canonical)
# 2. Do @meta.data columns match Supplementary Table 1?
# 3. Total cell count (expect 803,283)
# 4. Cell distribution across 58 donors
# 5. Is sample_type preserved?
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(yaml)
})

# Define %+% operator early
`%+%` <- function(a, b) paste0(a, b)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

RDS_PATH <- "${SOURCE_COMPONENT_STUDIES}/reed.rds"
OUTPUT_DIR <- "${SOURCE_IHBCAV1_HARMONIZATION%/}/studies/reed/extracted"

# Ensure output directory exists
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# -----------------------------------------------------------------------------
# Load RDS
# -----------------------------------------------------------------------------

message("=" %+% strrep("=", 79))
message("Loading Reed RDS file...")
message("Source: ", RDS_PATH)
message("=" %+% strrep("=", 79))

load_start <- Sys.time()
srt <- readRDS(RDS_PATH)
load_end <- Sys.time()

message("Loaded in ", round(difftime(load_end, load_start, units = "mins"), 2), " minutes")

# Check if UpdateSeuratObject is needed
if (!inherits(srt, "Seurat") || is.null(srt@meta.data)) {
  message("Attempting to update Seurat object...")
  srt <- UpdateSeuratObject(srt)
}

# -----------------------------------------------------------------------------
# Basic Object Info
# -----------------------------------------------------------------------------

results <- list(
  extraction_timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  source_file = RDS_PATH,
  seurat_version = as.character(packageVersion("Seurat")),
  object_class = class(srt)[1]
)

# Total cell count
n_cells <- ncol(srt)
n_features <- nrow(srt)
results$total_cells <- n_cells
results$total_features <- n_features
results$expected_cells <- 803283
results$cell_count_matches <- n_cells == 803283

message("\n--- BASIC COUNTS ---")
message("Total cells: ", format(n_cells, big.mark = ","))
message("Total features: ", format(n_features, big.mark = ","))
message("Expected cells: 803,283")
message("Match: ", results$cell_count_matches)

# -----------------------------------------------------------------------------
# Cell ID Analysis (CANONICAL FORMAT DEFINITION)
# -----------------------------------------------------------------------------

message("\n--- CELL ID ANALYSIS (Defines Canonical Format) ---")

cell_ids <- colnames(srt)
results$cell_id_analysis <- list()

# Sample cell IDs
set.seed(42)
sample_ids <- cell_ids[sample(length(cell_ids), min(20, length(cell_ids)))]
results$cell_id_analysis$sample_ids <- sample_ids

message("Sample cell IDs:")
for (id in sample_ids[1:10]) {
  message("  ", id)
}

# Check for -N suffix pattern (10x barcode format)
has_suffix <- grepl("-[0-9]+$", cell_ids)
results$cell_id_analysis$has_suffix_fraction <- mean(has_suffix)
message("\nFraction with -N suffix: ", round(mean(has_suffix) * 100, 2), "%")

# Identify components by splitting on common delimiters
# Try underscore first
sample_split_underscore <- strsplit(sample_ids[1], "_")[[1]]
message("\nFirst ID split by underscore: ", paste(sample_split_underscore, collapse = " | "))

# Count underscore-separated components
n_components <- sapply(cell_ids, function(x) length(strsplit(x, "_")[[1]]))
results$cell_id_analysis$component_distribution <- table(n_components)
message("\nComponent count distribution (split by _):")
print(table(n_components))

# If mostly consistent component count, analyze pattern
modal_components <- as.numeric(names(which.max(table(n_components))))
message("\nModal component count: ", modal_components)

# For IDs with modal components, analyze each position
ids_with_modal <- cell_ids[n_components == modal_components]
if (length(ids_with_modal) > 0) {
  # Split all into components
  splits <- strsplit(ids_with_modal, "_")

  # Analyze each position
  results$cell_id_analysis$component_analysis <- list()
  for (i in 1:modal_components) {
    pos_values <- sapply(splits, `[`, i)
    unique_count <- length(unique(pos_values))

    # Check if this looks like barcode (16 chars, ACGT only)
    is_barcode_like <- all(nchar(gsub("-[0-9]+$", "", pos_values)) == 16) &&
                       all(grepl("^[ACGT]+(-[0-9]+)?$", pos_values))

    results$cell_id_analysis$component_analysis[[paste0("position_", i)]] <- list(
      unique_values = unique_count,
      is_barcode_like = is_barcode_like,
      sample_values = unique(pos_values)[1:min(5, unique_count)]
    )

    message(sprintf("\nPosition %d: %d unique values", i, unique_count))
    message("  Sample: ", paste(unique(pos_values)[1:5], collapse = ", "))
    if (is_barcode_like) message("  -> Looks like 10x barcode")
  }
}

# Document the canonical format
message("\n--- CANONICAL CELL ID FORMAT ---")
# Infer format from analysis
if (modal_components >= 2) {
  # Typical pattern: donor_barcode or donor_type_barcode
  splits <- strsplit(sample_ids, "_")
  message("Inferred format structure:")
  for (i in 1:modal_components) {
    pos_values <- sapply(splits, `[`, i)
    # Guess what each position represents
    if (all(nchar(gsub("-[0-9]+$", "", pos_values)) == 16) &&
        all(grepl("^[ACGT]+(-[0-9]+)?$", pos_values))) {
      message(sprintf("  Position %d: BARCODE", i))
    } else if (all(pos_values %in% c("S", "O", "T", "P"))) {
      message(sprintf("  Position %d: SAMPLE_TYPE (S=Supernatant, O=Organoid, T=Tissue, P=PDX?)", i))
    } else {
      unique_ct <- length(unique(pos_values))
      message(sprintf("  Position %d: %d unique values (likely donor or condition)", i, unique_ct))
    }
  }
}

# -----------------------------------------------------------------------------
# Metadata Column Analysis
# -----------------------------------------------------------------------------

message("\n--- METADATA COLUMNS ---")

meta <- srt@meta.data
meta_cols <- colnames(meta)
results$metadata_columns <- meta_cols

message("Total columns: ", length(meta_cols))
message("\nColumn names:")
for (col in meta_cols) {
  message("  - ", col)
}

# Detailed analysis of each column
results$metadata_details <- list()
message("\n--- METADATA COLUMN DETAILS ---")

for (col in meta_cols) {
  values <- meta[[col]]
  col_info <- list(
    class = class(values)[1],
    n_missing = sum(is.na(values)),
    pct_missing = round(mean(is.na(values)) * 100, 2)
  )

  if (is.numeric(values)) {
    col_info$type <- "numeric"
    col_info$min <- min(values, na.rm = TRUE)
    col_info$max <- max(values, na.rm = TRUE)
    col_info$mean <- round(mean(values, na.rm = TRUE), 4)
  } else {
    # Categorical
    value_table <- table(values, useNA = "ifany")
    col_info$type <- "categorical"
    col_info$n_unique <- length(unique(values))
    col_info$value_counts <- as.list(value_table)

    # If few unique values, list them all
    if (col_info$n_unique <= 50) {
      col_info$all_values <- names(value_table)
    }
  }

  results$metadata_details[[col]] <- col_info

  # Print summary
  message("\n", col, " (", col_info$class, "):")
  message("  Missing: ", col_info$n_missing, " (", col_info$pct_missing, "%)")

  if (col_info$type == "numeric") {
    message("  Range: [", col_info$min, ", ", col_info$max, "]")
  } else {
    message("  Unique values: ", col_info$n_unique)
    if (col_info$n_unique <= 20) {
      message("  Values: ", paste(names(value_table), collapse = ", "))
      # Show counts for small categories
      for (v in names(value_table)) {
        message("    ", v, ": ", value_table[v])
      }
    }
  }
}

# -----------------------------------------------------------------------------
# Donor Distribution
# -----------------------------------------------------------------------------

message("\n--- DONOR DISTRIBUTION ---")

# Find donor column (try common names)
donor_col <- NULL
for (candidate in c("donor", "Donor", "donor_id", "sample", "orig.ident")) {
  if (candidate %in% meta_cols) {
    donor_col <- candidate
    break
  }
}

if (!is.null(donor_col)) {
  donor_table <- table(meta[[donor_col]])
  results$donor_distribution <- list(
    column_name = donor_col,
    n_donors = length(donor_table),
    expected_donors = 58,
    matches_expected = length(donor_table) == 58,
    donor_counts = as.list(donor_table)
  )

  message("Donor column: ", donor_col)
  message("Number of donors: ", length(donor_table))
  message("Expected: 58")
  message("Match: ", results$donor_distribution$matches_expected)
  message("\nCells per donor:")
  for (d in names(sort(donor_table))) {
    message(sprintf("  %s: %s", d, format(donor_table[d], big.mark = ",")))
  }
} else {
  message("WARNING: Could not identify donor column")
  results$donor_distribution <- list(error = "Could not identify donor column")
}

# -----------------------------------------------------------------------------
# Sample Type Check
# -----------------------------------------------------------------------------

message("\n--- SAMPLE TYPE CHECK ---")

# Check for sample_type or similar column
sample_type_col <- NULL
for (candidate in c("sample_type", "Sample_type", "sampleType", "type")) {
  if (candidate %in% meta_cols) {
    sample_type_col <- candidate
    break
  }
}

if (!is.null(sample_type_col)) {
  sample_type_table <- table(meta[[sample_type_col]])
  results$sample_type <- list(
    column_name = sample_type_col,
    preserved = TRUE,
    value_counts = as.list(sample_type_table)
  )

  message("Sample type column: ", sample_type_col)
  message("Values:")
  for (st in names(sample_type_table)) {
    message(sprintf("  %s: %s", st, format(sample_type_table[st], big.mark = ",")))
  }
} else {
  message("No explicit sample_type column found")
  message("Checking if sample type is encoded in cell IDs or other columns...")
  results$sample_type <- list(
    column_name = NA,
    preserved = FALSE,
    note = "No explicit sample_type column; may be encoded elsewhere"
  )
}

# -----------------------------------------------------------------------------
# Assay Information
# -----------------------------------------------------------------------------

message("\n--- ASSAY INFORMATION ---")

results$assays <- list(
  available = Assays(srt),
  default = DefaultAssay(srt)
)

message("Available assays: ", paste(Assays(srt), collapse = ", "))
message("Default assay: ", DefaultAssay(srt))

# Check for reductions
if (length(Reductions(srt)) > 0) {
  results$reductions <- Reductions(srt)
  message("Available reductions: ", paste(Reductions(srt), collapse = ", "))
} else {
  results$reductions <- NULL
  message("No dimensional reductions found")
}

# -----------------------------------------------------------------------------
# Summary of Canonical Definitions
# -----------------------------------------------------------------------------

message("\n" %+% strrep("=", 80))
message("CANONICAL DEFINITIONS FROM REED (iHBCA Reference)")
message(strrep("=", 80))

results$canonical_definitions <- list(
  cell_id = list(
    format = "To be documented from analysis above",
    example = sample_ids[1],
    note = "Reed cell IDs define the canonical format for iHBCA"
  ),
  metadata_columns = meta_cols,
  n_donors = if (!is.null(donor_col)) length(unique(meta[[donor_col]])) else NA,
  n_cells = n_cells,
  n_features = n_features
)

# -----------------------------------------------------------------------------
# Write Output YAML
# -----------------------------------------------------------------------------

output_file <- file.path(OUTPUT_DIR, "raw_data_extraction.yaml")
message("\nWriting results to: ", output_file)

# Convert results to YAML-friendly format
# (yaml::write_yaml has issues with some R types)
yaml_output <- yaml::as.yaml(results, indent = 2, indent.mapping.sequence = TRUE)
writeLines(yaml_output, output_file)

message("\nExtraction complete!")
message("Output: ", output_file)
