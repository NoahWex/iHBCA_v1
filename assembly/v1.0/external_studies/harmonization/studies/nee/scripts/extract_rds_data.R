# =============================================================================
# NEE RDS DATA EXTRACTION
# =============================================================================
# UNSUPERVISED extraction from Nee RDS file
# Focus: Verify cell ID patterns, metadata columns, and donor distribution
#
# Key questions to answer (from edge_cases.yaml):
# 1. cell_id_format_verification: Does Ctrl-*/BRCA1-* pattern hold for ALL cells?
# 2. cell_id_full_format: Complete format beyond prefix?
# 3. patient_id_to_cell_mapping: How do Table 1 Patient IDs map to cell IDs?
# 4. metadata_column_mapping: How do supplemental columns map to @meta.data?
# 5. parity_format_parsing: Is parity in GPMA format?
# 6. age_as_testable_condition: Confirm age distribution (18-71 expected)
# 7. cell_distribution: How do 230,100 cells distribute across 22 donors?
# =============================================================================

library(Seurat)
library(yaml)
library(dplyr)
library(stringr)

# =============================================================================
# CONFIGURATION
# =============================================================================

INPUT_RDS <- "${PROJECT_SPATIAL_HBCA}/project/ReferenceDatasets/iHBCA/component_studies/nee.rds"
OUTPUT_YAML <- "${SOURCE_IHBCAV1_HARMONIZATION%/}/studies/nee/extracted/raw_data_extraction.yaml"

cat("=============================================================================\n")
cat("NEE RDS EXTRACTION - UNSUPERVISED ANALYSIS\n")
cat("=============================================================================\n")
cat("Input:", INPUT_RDS, "\n")
cat("Output:", OUTPUT_YAML, "\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=============================================================================\n\n")

# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading RDS file...\n")
obj <- readRDS(INPUT_RDS)

# Update Seurat object if needed
if (inherits(obj, "Seurat")) {
  obj <- UpdateSeuratObject(obj)
}

cat("Object type:", class(obj)[1], "\n")
cat("Total cells:", ncol(obj), "\n\n")

# =============================================================================
# CELL ID ANALYSIS
# =============================================================================

cat("=============================================================================\n")
cat("CELL ID PATTERN ANALYSIS\n")
cat("=============================================================================\n")

cell_ids <- colnames(obj)
n_cells <- length(cell_ids)

cat("Total cell IDs:", n_cells, "\n\n")

# ----- PREFIX ANALYSIS (CRITICAL FOR BRCA FIX) -----
cat("--- PREFIX ANALYSIS (Ctrl- vs BRCA1-) ---\n")

ctrl_mask <- grepl("^Ctrl-", cell_ids)
brca1_mask <- grepl("^BRCA1-", cell_ids)
other_mask <- !ctrl_mask & !brca1_mask

ctrl_count <- sum(ctrl_mask)
brca1_count <- sum(brca1_mask)
other_count <- sum(other_mask)

cat("Cells with Ctrl- prefix:", ctrl_count, "(", round(100*ctrl_count/n_cells, 2), "%)\n")
cat("Cells with BRCA1- prefix:", brca1_count, "(", round(100*brca1_count/n_cells, 2), "%)\n")
cat("Cells with OTHER prefix:", other_count, "(", round(100*other_count/n_cells, 2), "%)\n")

# If there are "other" cells, show examples
other_examples <- character(0)
if (other_count > 0) {
  other_examples <- head(cell_ids[other_mask], 20)
  cat("\nWARNING: Found cells without Ctrl- or BRCA1- prefix!\n")
  cat("Examples of OTHER patterns:\n")
  print(other_examples)
}

# ----- FULL FORMAT ANALYSIS -----
cat("\n--- FULL CELL ID FORMAT ANALYSIS ---\n")

# Sample cell IDs for pattern analysis
sample_ctrl <- head(cell_ids[ctrl_mask], 10)
sample_brca1 <- head(cell_ids[brca1_mask], 10)

cat("Sample Ctrl- cell IDs:\n")
print(sample_ctrl)

cat("\nSample BRCA1- cell IDs:\n")
print(sample_brca1)

# Parse cell ID components (assuming format: {prefix}-{donor}_{barcode})
parse_cell_id <- function(cell_id) {
  # Try pattern: PREFIX-DONOR_BARCODE
  match <- str_match(cell_id, "^(Ctrl|BRCA1)-([^_]+)_(.+)$")
  if (!is.na(match[1])) {
    return(list(
      prefix = match[2],
      donor = match[3],
      barcode = match[4],
      parsed = TRUE
    ))
  }
  return(list(
    prefix = NA,
    donor = NA,
    barcode = NA,
    parsed = FALSE
  ))
}

# Parse all cell IDs
cat("\n--- PARSING ALL CELL IDs ---\n")
parsed <- lapply(cell_ids, parse_cell_id)
prefixes <- sapply(parsed, `[[`, "prefix")
donors <- sapply(parsed, `[[`, "donor")
barcodes <- sapply(parsed, `[[`, "barcode")
parsed_success <- sapply(parsed, `[[`, "parsed")

cat("Successfully parsed:", sum(parsed_success), "/", n_cells, "\n")
cat("Failed to parse:", sum(!parsed_success), "\n")

if (sum(!parsed_success) > 0) {
  cat("\nExamples of unparseable cell IDs:\n")
  print(head(cell_ids[!parsed_success], 10))
}

# Donor analysis
cat("\n--- DONOR ANALYSIS FROM CELL IDs ---\n")
unique_donors <- unique(donors[!is.na(donors)])
cat("Number of unique donors:", length(unique_donors), "\n")
cat("Unique donors:\n")
print(sort(unique_donors))

# Cells per donor
donor_counts <- table(donors)
cat("\nCells per donor:\n")
print(sort(donor_counts, decreasing = TRUE))

# Donor by prefix
donor_prefix_table <- table(prefixes, donors)
cat("\nDonor distribution by prefix:\n")

ctrl_donors <- unique(donors[prefixes == "Ctrl" & !is.na(prefixes)])
brca1_donors <- unique(donors[prefixes == "BRCA1" & !is.na(prefixes)])
cat("\nCtrl donors (", length(ctrl_donors), "):\n")
print(sort(ctrl_donors))
cat("\nBRCA1 donors (", length(brca1_donors), "):\n")
print(sort(brca1_donors))

# Barcode analysis
cat("\n--- BARCODE ANALYSIS ---\n")
barcode_lengths <- nchar(barcodes[!is.na(barcodes)])
cat("Barcode lengths:\n")
print(table(barcode_lengths))

# Check if barcodes are ACGT only
acgt_only <- grepl("^[ACGT]+$", barcodes[!is.na(barcodes)])
cat("\nBarcodes with ACGT only:", sum(acgt_only), "/", sum(!is.na(barcodes)), "\n")

# Check for -1 suffix
has_suffix <- grepl("-[0-9]+$", barcodes[!is.na(barcodes)])
cat("Barcodes with -N suffix:", sum(has_suffix), "/", sum(!is.na(barcodes)), "\n")

# =============================================================================
# METADATA ANALYSIS
# =============================================================================

cat("\n=============================================================================\n")
cat("METADATA ANALYSIS\n")
cat("=============================================================================\n")

meta <- obj@meta.data
cat("Number of metadata columns:", ncol(meta), "\n")
cat("Metadata column names:\n")
print(colnames(meta))

# Analyze each metadata column
metadata_summary <- list()

for (col in colnames(meta)) {
  cat("\n--- Column:", col, "---\n")
  values <- meta[[col]]

  col_info <- list(
    name = col,
    class = class(values)[1],
    n_missing = sum(is.na(values)),
    n_unique = length(unique(values[!is.na(values)]))
  )

  if (is.numeric(values)) {
    col_info$type <- "numeric"
    col_info$min <- min(values, na.rm = TRUE)
    col_info$max <- max(values, na.rm = TRUE)
    col_info$mean <- mean(values, na.rm = TRUE)
    col_info$median <- median(values, na.rm = TRUE)

    cat("  Type: numeric\n")
    cat("  Range:", col_info$min, "-", col_info$max, "\n")
    cat("  Mean:", round(col_info$mean, 2), "\n")
    cat("  Missing:", col_info$n_missing, "\n")

    # Show distribution for likely categorical numerics
    if (col_info$n_unique <= 30) {
      cat("  Value counts:\n")
      print(table(values, useNA = "ifany"))
    }

  } else {
    col_info$type <- "categorical"
    unique_vals <- unique(values)
    col_info$unique_values <- as.character(head(sort(unique_vals), 50))

    cat("  Type: categorical\n")
    cat("  Unique values:", length(unique_vals), "\n")
    cat("  Missing:", col_info$n_missing, "\n")

    if (length(unique_vals) <= 50) {
      cat("  Value counts:\n")
      print(table(values, useNA = "ifany"))
    } else {
      cat("  Sample values:\n")
      print(head(unique_vals, 20))
    }
  }

  metadata_summary[[col]] <- col_info
}

# =============================================================================
# SPECIFIC QUESTION CHECKS
# =============================================================================

cat("\n=============================================================================\n")
cat("SPECIFIC QUESTION CHECKS\n")
cat("=============================================================================\n")

# Check for risk_status or BRCA column (to verify mislabeling issue)
cat("\n--- BRCA/RISK STATUS COLUMNS ---\n")
brca_cols <- grep("brca|risk|carrier|mutation", colnames(meta), ignore.case = TRUE, value = TRUE)
if (length(brca_cols) > 0) {
  cat("Found potential BRCA-related columns:\n")
  for (col in brca_cols) {
    cat("\n  Column:", col, "\n")
    print(table(meta[[col]], useNA = "ifany"))
  }
} else {
  cat("No columns with 'brca', 'risk', 'carrier', or 'mutation' found\n")
}

# Check for age column
cat("\n--- AGE COLUMN ---\n")
age_cols <- grep("age", colnames(meta), ignore.case = TRUE, value = TRUE)
if (length(age_cols) > 0) {
  cat("Found potential age columns:\n")
  for (col in age_cols) {
    cat("\n  Column:", col, "\n")
    values <- meta[[col]]
    if (is.numeric(values)) {
      cat("  Range:", min(values, na.rm = TRUE), "-", max(values, na.rm = TRUE), "\n")
      cat("  Cells <50:", sum(values < 50, na.rm = TRUE), "\n")
      cat("  Cells >=50:", sum(values >= 50, na.rm = TRUE), "\n")
      cat("  Missing:", sum(is.na(values)), "\n")
    } else {
      print(table(values, useNA = "ifany"))
    }
  }
} else {
  cat("No columns with 'age' found\n")
}

# Check for parity column
cat("\n--- PARITY COLUMN ---\n")
parity_cols <- grep("parity|gravid|para|pregnancy", colnames(meta), ignore.case = TRUE, value = TRUE)
if (length(parity_cols) > 0) {
  cat("Found potential parity columns:\n")
  for (col in parity_cols) {
    cat("\n  Column:", col, "\n")
    print(head(table(meta[[col]], useNA = "ifany"), 30))
  }
} else {
  cat("No columns with 'parity', 'gravid', 'para', or 'pregnancy' found\n")
}

# Check for donor/patient column
cat("\n--- DONOR/PATIENT COLUMN ---\n")
donor_cols <- grep("donor|patient|sample|individual|subject", colnames(meta), ignore.case = TRUE, value = TRUE)
if (length(donor_cols) > 0) {
  cat("Found potential donor columns:\n")
  for (col in donor_cols) {
    cat("\n  Column:", col, "\n")
    values <- meta[[col]]
    unique_vals <- unique(values)
    cat("  Unique values:", length(unique_vals), "\n")
    if (length(unique_vals) <= 50) {
      print(table(values, useNA = "ifany"))
    } else {
      print(head(unique_vals, 30))
    }
  }
} else {
  cat("No columns with 'donor', 'patient', 'sample', 'individual', or 'subject' found\n")
}

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

cat("\n=============================================================================\n")
cat("SUMMARY STATISTICS\n")
cat("=============================================================================\n")

summary_stats <- list(
  total_cells = n_cells,
  expected_cells = 230100,
  cells_match_expected = n_cells == 230100,

  cell_id_pattern = list(
    ctrl_prefix_count = ctrl_count,
    brca1_prefix_count = brca1_count,
    other_prefix_count = other_count,
    pattern_coverage = round(100 * (ctrl_count + brca1_count) / n_cells, 2),
    all_cells_have_prefix = other_count == 0
  ),

  donors = list(
    total_unique = length(unique_donors),
    ctrl_donors = length(ctrl_donors),
    brca1_donors = length(brca1_donors),
    ctrl_donor_ids = sort(ctrl_donors),
    brca1_donor_ids = sort(brca1_donors)
  ),

  parsing = list(
    success_count = sum(parsed_success),
    failure_count = sum(!parsed_success),
    success_rate = round(100 * sum(parsed_success) / n_cells, 2)
  ),

  metadata_columns = colnames(meta)
)

cat("Total cells:", summary_stats$total_cells, "\n")
cat("Expected cells:", summary_stats$expected_cells, "\n")
cat("Match expected:", summary_stats$cells_match_expected, "\n\n")

cat("Cell ID prefix coverage:\n")
cat("  Ctrl-:", summary_stats$cell_id_pattern$ctrl_prefix_count, "\n")
cat("  BRCA1-:", summary_stats$cell_id_pattern$brca1_prefix_count, "\n")
cat("  Other:", summary_stats$cell_id_pattern$other_prefix_count, "\n")
cat("  Pattern coverage:", summary_stats$cell_id_pattern$pattern_coverage, "%\n")
cat("  All cells have Ctrl/BRCA1 prefix:", summary_stats$cell_id_pattern$all_cells_have_prefix, "\n\n")

cat("Donors:\n")
cat("  Total unique:", summary_stats$donors$total_unique, "\n")
cat("  Ctrl donors:", summary_stats$donors$ctrl_donors, "\n")
cat("  BRCA1 donors:", summary_stats$donors$brca1_donors, "\n")

# =============================================================================
# BUILD OUTPUT YAML
# =============================================================================

cat("\n=============================================================================\n")
cat("BUILDING OUTPUT YAML\n")
cat("=============================================================================\n")

# Build comprehensive output structure
output <- list(
  extraction_info = list(
    source_file = INPUT_RDS,
    extraction_date = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    seurat_version = as.character(packageVersion("Seurat")),
    r_version = R.version.string
  ),

  cell_count = list(
    total = n_cells,
    expected = 230100,
    matches_expected = n_cells == 230100
  ),

  cell_id_analysis = list(
    prefix_pattern = list(
      ctrl_count = ctrl_count,
      brca1_count = brca1_count,
      other_count = other_count,
      pattern_coverage_pct = round(100 * (ctrl_count + brca1_count) / n_cells, 4),
      all_have_valid_prefix = other_count == 0,
      conclusion = ifelse(other_count == 0,
                          "CONFIRMED: All cells have Ctrl- or BRCA1- prefix",
                          paste0("WARNING: ", other_count, " cells lack expected prefix"))
    ),

    format = list(
      inferred_pattern = "{prefix}-{donor_id}_{barcode}",
      prefix_values = c("Ctrl", "BRCA1"),
      parse_success_count = sum(parsed_success),
      parse_failure_count = sum(!parsed_success),
      unparseable_examples = if(sum(!parsed_success) > 0) head(cell_ids[!parsed_success], 10) else character(0)
    ),

    barcode = list(
      lengths_observed = as.list(table(barcode_lengths)),
      acgt_only_count = sum(acgt_only),
      has_suffix_count = sum(has_suffix),
      suffix_pattern = ifelse(sum(has_suffix) > 0, "-N suffix present", "No suffix")
    ),

    examples = list(
      ctrl_examples = sample_ctrl,
      brca1_examples = sample_brca1,
      other_examples = other_examples
    )
  ),

  donor_analysis = list(
    total_unique_donors = length(unique_donors),
    ctrl_donor_count = length(ctrl_donors),
    brca1_donor_count = length(brca1_donors),
    expected_donor_count = 22,
    matches_expected = length(unique_donors) == 22,

    ctrl_donor_ids = sort(ctrl_donors),
    brca1_donor_ids = sort(brca1_donors),

    cells_per_donor = as.list(sort(donor_counts, decreasing = TRUE))
  ),

  metadata = list(
    column_count = ncol(meta),
    column_names = colnames(meta),
    column_details = metadata_summary
  ),

  question_answers = list(
    cell_id_format_verification = list(
      question = "Does cell ID prefix pattern (Ctrl-*/BRCA1-*) hold for ALL cells?",
      answer = ifelse(other_count == 0, "YES", "NO"),
      ctrl_cells = ctrl_count,
      brca1_cells = brca1_count,
      other_cells = other_count,
      coverage_pct = round(100 * (ctrl_count + brca1_count) / n_cells, 4)
    ),

    cell_id_full_format = list(
      question = "What is the complete cell ID format?",
      inferred_format = "{prefix}-{donor_id}_{barcode}",
      components = c("prefix (Ctrl or BRCA1)", "donor_id", "barcode"),
      barcode_format = paste0("ACGT only: ", sum(acgt_only) == sum(!is.na(barcodes)),
                              ", with -N suffix: ", sum(has_suffix) > 0)
    ),

    cell_distribution = list(
      question = "Cell distribution across donors?",
      total_cells = n_cells,
      total_donors = length(unique_donors),
      min_cells_per_donor = min(donor_counts),
      max_cells_per_donor = max(donor_counts),
      mean_cells_per_donor = round(mean(donor_counts), 1)
    )
  ),

  brca_fix_validation = list(
    purpose = "Validate that cell ID prefix can be used to fix BRCA mislabeling",

    prefix_coverage = list(
      all_cells_covered = other_count == 0,
      coverage_pct = round(100 * (ctrl_count + brca1_count) / n_cells, 4)
    ),

    expected_after_fix = list(
      ar_donors = length(ctrl_donors),
      hr_donors = length(brca1_donors),
      ar_cells = ctrl_count,
      hr_cells = brca1_count
    ),

    can_apply_fix = other_count == 0 && length(ctrl_donors) > 0 && length(brca1_donors) > 0,

    recommendation = ifelse(
      other_count == 0 && length(ctrl_donors) > 0 && length(brca1_donors) > 0,
      "PROCEED: Cell ID prefix reliably encodes BRCA status for all cells",
      "INVESTIGATE: Some cells may not have valid prefix"
    )
  )
)

# Write YAML
cat("Writing output to:", OUTPUT_YAML, "\n")

# Ensure output directory exists
dir.create(dirname(OUTPUT_YAML), recursive = TRUE, showWarnings = FALSE)

# Write YAML with proper formatting
write_yaml(output, OUTPUT_YAML)

cat("\n=============================================================================\n")
cat("EXTRACTION COMPLETE\n")
cat("=============================================================================\n")
cat("Output written to:", OUTPUT_YAML, "\n")
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
