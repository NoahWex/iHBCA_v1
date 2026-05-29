#!/usr/bin/env Rscript
# =============================================================================
# PAL NORM B1: UNSUPERVISED RDS EXTRACTION
# =============================================================================
# Purpose: Extract ALL cell IDs, metadata, and structure from the NormB1Total
#          Seurat RDS file using UNSUPERVISED pattern analysis.
#
# Open Questions to Answer:
#   1. cell_id_format_b1: Confirm BRCA1 cell ID format (B1_{donor}_{barcode}-1?)
#   2. cell_count_verification: Verify total (59,766), BRCA1 (23,240), Normal (36,526)
#   3. donor_0123_exclusion_verification: Count N_0123_* cells
#   4. metadata_columns_in_rds: What metadata columns exist?
#   5. cell_distribution_by_donor: How do cells distribute across donors?
#   6. ihbca_donor_mapping_verification: Verify donor ID mappings
#
# Output: studies/pal_norm_b1/extracted/raw_data_extraction.yaml
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(yaml)
})

# Configuration
INPUT_FILE <- "${SOURCE_PAL_ORIGINAL%/}/SeuratObject_NormB1Total.rds"
OUTPUT_FILE <- "${SOURCE_IHBCAV1_HARMONIZATION%/}/studies/pal_norm_b1/extracted/raw_data_extraction.yaml"

cat("=" , rep("=", 60), "\n", sep="")
cat("PAL NORM B1: UNSUPERVISED RDS EXTRACTION\n")
cat("=" , rep("=", 60), "\n", sep="")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Input:", INPUT_FILE, "\n")
cat("Output:", OUTPUT_FILE, "\n\n")

# Load RDS
cat("Loading Seurat object...\n")
srt <- tryCatch({
  readRDS(INPUT_FILE)
}, error = function(e) {
  stop(sprintf("Failed to load RDS: %s", e$message))
})

# Update Seurat object if needed
if (inherits(srt, "Seurat")) {
  srt <- tryCatch(UpdateSeuratObject(srt), error = function(e) srt)
}
cat("Loaded successfully.\n\n")

# Basic info
cat("=== Basic Info ===\n")
total_cells <- ncol(srt)
cat("Total cells:", total_cells, "\n")
cat("Total features:", nrow(srt), "\n\n")

# Get cell IDs
cell_ids <- colnames(srt)

# Unsupervised cell ID pattern analysis
cat("=== Cell ID Pattern Analysis ===\n")

# Extract prefix (first component before _)
prefixes <- sapply(strsplit(cell_ids, "_"), `[`, 1)
prefix_counts <- table(prefixes)
cat("Unique prefixes:\n")
print(prefix_counts)
cat("\n")

# Separate by prefix type
b1_cells <- cell_ids[prefixes == "B1"]
n_cells <- cell_ids[prefixes == "N"]

cat("B1_* cells:", length(b1_cells), "\n")
cat("N_* cells:", length(n_cells), "\n\n")

# Check for N_0123_* specifically (collision issue)
n_0123_cells <- n_cells[grepl("^N_0123_", n_cells)]
cat("N_0123_* cells (collision concern):", length(n_0123_cells), "\n\n")

# Analyze cell ID structure per prefix
cat("=== Cell ID Structure Analysis ===\n")

analyze_cell_ids <- function(ids, prefix_name) {
  if (length(ids) == 0) return(NULL)

  cat("\n--- ", prefix_name, " cells (n=", length(ids), ") ---\n", sep="")

  # Show examples
  cat("Examples:\n")
  print(head(ids, 5))

  # Count components
  parts <- strsplit(ids, "_")
  n_parts <- sapply(parts, length)
  cat("\nComponent counts:\n")
  print(table(n_parts))

  # Extract components by position
  max_parts <- max(n_parts)
  components <- list()
  for (i in 1:max_parts) {
    components[[paste0("part", i)]] <- sapply(parts, function(x) if(length(x) >= i) x[i] else NA)
  }

  # Show unique values per component
  cat("\nUnique values per component:\n")
  for (i in 1:max_parts) {
    vals <- na.omit(components[[paste0("part", i)]])
    unique_vals <- unique(vals)
    cat("  Part", i, ":", length(unique_vals), "unique values\n")
    if (length(unique_vals) <= 20) {
      cat("    Values:", paste(head(sort(unique_vals), 20), collapse=", "), "\n")
    }
  }

  # Get per-donor counts (using second component as donor)
  donors <- sapply(parts, function(x) if(length(x) >= 2) x[2] else NA)
  donor_counts <- as.list(table(donors))

  return(list(
    count = length(ids),
    examples = head(ids, 10),
    n_components = as.list(table(n_parts)),
    unique_per_component = lapply(components, function(x) length(unique(na.omit(x)))),
    donor_counts = donor_counts
  ))
}

b1_analysis <- analyze_cell_ids(b1_cells, "B1")
n_analysis <- analyze_cell_ids(n_cells, "N")

# Metadata columns
cat("\n=== Metadata Columns ===\n")
meta <- srt@meta.data
meta_cols <- colnames(meta)
cat("Number of columns:", length(meta_cols), "\n")
cat("Column names:\n")
print(meta_cols)

# For each metadata column, show unique values and counts
cat("\n=== Metadata Column Details ===\n")
meta_details <- list()
for (col in meta_cols) {
  cat("\n--- Column:", col, "---\n")
  vals <- meta[[col]]

  if (is.numeric(vals)) {
    cat("Type: numeric\n")
    cat("Range:", min(vals, na.rm=TRUE), "-", max(vals, na.rm=TRUE), "\n")
    cat("Mean:", mean(vals, na.rm=TRUE), "\n")
    cat("NAs:", sum(is.na(vals)), "\n")
    meta_details[[col]] <- list(
      type = "numeric",
      min = min(vals, na.rm=TRUE),
      max = max(vals, na.rm=TRUE),
      mean = mean(vals, na.rm=TRUE),
      na_count = sum(is.na(vals))
    )
  } else {
    cat("Type: categorical\n")
    val_counts <- table(vals, useNA = "ifany")
    cat("Unique values:", length(val_counts), "\n")
    if (length(val_counts) <= 30) {
      cat("Value counts:\n")
      print(sort(val_counts, decreasing = TRUE))
      meta_details[[col]] <- list(
        type = "categorical",
        n_unique = length(val_counts),
        value_counts = as.list(val_counts)
      )
    } else {
      cat("(Too many unique values to display - showing top 20)\n")
      print(head(sort(val_counts, decreasing = TRUE), 20))
      meta_details[[col]] <- list(
        type = "categorical",
        n_unique = length(val_counts),
        top_values = as.list(head(sort(val_counts, decreasing = TRUE), 20))
      )
    }
  }
}

# Check for risk status or BRCA indicator columns
cat("\n=== Risk Status / BRCA Indicator Columns ===\n")
risk_keywords <- c("risk", "brca", "status", "condition", "group", "type", "category", "mutation")
risk_cols <- meta_cols[sapply(tolower(meta_cols), function(cn) any(sapply(risk_keywords, function(kw) grepl(kw, cn))))]

risk_indicators <- list()
if (length(risk_cols) > 0) {
  cat("Found potential risk/BRCA columns:\n")
  for (rc in risk_cols) {
    cat("\n--- Column:", rc, "---\n")
    val_counts <- table(meta[[rc]], useNA = "ifany")
    print(val_counts)
    risk_indicators[[rc]] <- as.list(val_counts)
  }
} else {
  cat("No columns with risk/BRCA keywords found in metadata.\n")
  cat("Note: Risk status may need to be inferred from cell ID prefix (B1 vs N)\n")
}

# Cell distribution by donor (if donor column exists)
cat("\n=== Cell Distribution by Donor ===\n")
donor_distribution <- NULL

# Try common donor column names
donor_cols <- c("donor", "orig.ident", "sample", "patient", "donor_id", "Sample")
found_donor_col <- NULL
for (dc in donor_cols) {
  if (dc %in% meta_cols) {
    found_donor_col <- dc
    break
  }
}

if (!is.null(found_donor_col)) {
  cat("Using column:", found_donor_col, "\n")
  donor_distribution <- as.list(table(meta[[found_donor_col]]))
  print(table(meta[[found_donor_col]]))
} else {
  cat("No standard donor column found. Extracting from cell IDs...\n")

  # Extract donor from cell IDs based on pattern
  # B1 cells: B1_{donor}_{barcode}
  # N cells: N_{donor}_{type}_{barcode}

  extract_donor <- function(id) {
    parts <- strsplit(id, "_")[[1]]
    if (parts[1] == "B1") {
      return(paste0("B1_", parts[2]))
    } else if (parts[1] == "N") {
      return(paste0("N_", parts[2]))
    } else {
      return(parts[1])
    }
  }

  donors <- sapply(cell_ids, extract_donor)
  donor_distribution <- as.list(table(donors))
  cat("Donor distribution (extracted from cell IDs):\n")
  print(sort(table(donors), decreasing = TRUE))
}

# Create output YAML
cat("\n=== Writing Output ===\n")

output <- list(
  extraction_info = list(
    source_file = INPUT_FILE,
    extraction_date = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    seurat_version = as.character(packageVersion("Seurat"))
  ),

  summary = list(
    total_cells = total_cells,
    total_features = nrow(srt),
    b1_cells = length(b1_cells),
    n_cells = length(n_cells),
    n_0123_cells = length(n_0123_cells)
  ),

  cell_id_patterns = list(
    prefixes = as.list(prefix_counts),
    b1_analysis = b1_analysis,
    n_analysis = n_analysis,
    n_0123_examples = head(n_0123_cells, 10)
  ),

  metadata = list(
    columns = meta_cols,
    column_details = meta_details,
    risk_indicator_columns = if (length(risk_indicators) > 0) risk_indicators else list(note = "No risk/BRCA columns found; infer from cell ID prefix")
  ),

  donor_distribution = donor_distribution,

  # iHBCA donor mapping verification
  ihbca_mapping = list(
    note = "Expected B1 donor mappings (from edge_cases.yaml)",
    expected_mappings = list(
      "B1_0894" = "KCF0894",
      "B1_0033" = "MH0033",
      "B1_0023" = "MH0023",
      "B1_0090" = "MH0090"
    ),
    observed_b1_donors = if (!is.null(b1_analysis)) names(b1_analysis$donor_counts) else NULL
  ),

  answers_to_questions = list(
    cell_id_format_b1 = list(
      question = "Confirm BRCA1 cell ID format is B1_{donor}_{barcode}-1?",
      answer = "See cell_id_patterns.b1_analysis for BRCA1 format discovery",
      b1_examples = if (!is.null(b1_analysis)) b1_analysis$examples else NULL
    ),
    cell_count_verification = list(
      question = "Confirm total (59,766), BRCA1 (23,240), Normal (36,526)",
      observed_total = total_cells,
      expected_total = 59766,
      match_total = total_cells == 59766,
      observed_brca1 = length(b1_cells),
      expected_brca1 = 23240,
      match_brca1 = length(b1_cells) == 23240,
      observed_normal = length(n_cells),
      expected_normal = 36526,
      match_normal = length(n_cells) == 36526
    ),
    donor_0123_exclusion_verification = list(
      question = "How many cells are N_0123_*?",
      n_0123_count = length(n_0123_cells),
      examples = head(n_0123_cells, 5),
      note = "These cells would collide with B1_0023 if both map to MH0023"
    ),
    metadata_columns_in_rds = list(
      question = "What metadata columns exist?",
      columns = meta_cols,
      n_columns = length(meta_cols)
    ),
    cell_distribution_by_donor = list(
      question = "How do cells distribute across BRCA1 vs Normal donors?",
      distribution = donor_distribution,
      b1_donor_counts = if (!is.null(b1_analysis)) b1_analysis$donor_counts else NULL,
      n_donor_counts = if (!is.null(n_analysis)) n_analysis$donor_counts else NULL
    ),
    ihbca_donor_mapping_verification = list(
      question = "Can we verify donor ID mappings?",
      note = "See ihbca_mapping section for expected vs observed donors"
    ),
    risk_status_columns = list(
      question = "Is risk_status or BRCA indicator column present?",
      columns_found = if (length(risk_indicators) > 0) names(risk_indicators) else "None",
      note = if (length(risk_indicators) == 0) "Risk must be inferred from cell ID prefix (B1=BRCA1, N=Normal)" else "See metadata.risk_indicator_columns"
    )
  )
)

# Write YAML
write_yaml(output, OUTPUT_FILE)
cat("Output written to:", OUTPUT_FILE, "\n")

cat("\n", rep("=", 60), "\n", sep="")
cat("EXTRACTION COMPLETE\n")
cat(rep("=", 60), "\n", sep="")
cat("Total cells:", total_cells, "\n")
cat("B1 (BRCA1) cells:", length(b1_cells), "\n")
cat("N (Normal) cells:", length(n_cells), "\n")
cat("N_0123 cells (collision):", length(n_0123_cells), "\n")
cat("Output:", OUTPUT_FILE, "\n")
