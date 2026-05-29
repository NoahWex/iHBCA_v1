# =============================================================================
# MURROW RDS EXTRACTION - TRACK 2 WAVE 1
# =============================================================================
# Unsupervised extraction from murrow.rds to answer open questions
# Source: ${SOURCE_MURROW_ORIGINAL%/}/murrow.rds
# =============================================================================

library(Seurat)
library(yaml)

# =============================================================================
# Configuration
# =============================================================================

source_path <- "${SOURCE_COMPONENT_STUDIES}/murrow_l1_original.rds"
output_path <- "${SOURCE_IHBCAV1_HARMONIZATION%/}/studies/murrow/extracted/raw_data_extraction.yaml"

message("=== MURROW RDS EXTRACTION ===")
message("Source: ", source_path)
message("Output: ", output_path)

# =============================================================================
# Load Object
# =============================================================================

message("\n[1/5] Loading RDS file...")
t0 <- Sys.time()
obj <- readRDS(source_path)
message("  Loaded in ", round(difftime(Sys.time(), t0, units = "secs"), 1), " seconds")

# Check if needs updating
if (inherits(obj, "Seurat")) {
  tryCatch({
    obj <- UpdateSeuratObject(obj)
    message("  Updated Seurat object")
  }, error = function(e) {
    message("  Note: ", e$message)
  })
}

# =============================================================================
# Object Summary
# =============================================================================

message("\n[2/5] Extracting object summary...")

file_info <- file.info(source_path)

object_summary <- list(
  class = class(obj)[1],
  seurat_version = if (inherits(obj, "Seurat")) as.character(obj@version) else NA,
  n_cells = ncol(obj),
  n_features = nrow(obj),
  assays = if (inherits(obj, "Seurat")) Assays(obj) else NA,
  reductions = if (inherits(obj, "Seurat")) Reductions(obj) else NA
)

message("  Class: ", object_summary$class)
message("  Cells: ", format(object_summary$n_cells, big.mark = ","))
message("  Features: ", format(object_summary$n_features, big.mark = ","))

# =============================================================================
# Cell ID Extraction
# =============================================================================

message("\n[3/5] Extracting cell IDs...")

cell_ids <- colnames(obj)
n_cells <- length(cell_ids)

# Pattern analysis - split by potential separators
analyze_cell_ids <- function(ids) {
  result <- list()

  # Sample first 1000 for pattern analysis
  sample_ids <- if (length(ids) > 1000) sample(ids, 1000) else ids

  # Try different separators
  separators <- c("_", "-", ":")

  for (sep in separators) {
    splits <- strsplit(sample_ids, sep, fixed = TRUE)
    n_components <- sapply(splits, length)
    result[[paste0("sep_", sep)]] <- table(n_components)
  }

  # Identify barcode component (16bp ACGT)
  barcode_pattern <- "^[ACGT]{16}$"
  barcode_with_suffix <- "^[ACGT]{16}-[0-9]+$"

  # Check each component position for barcode
  underscore_splits <- strsplit(sample_ids, "_", fixed = TRUE)
  max_components <- max(sapply(underscore_splits, length))

  barcode_positions <- list()
  for (i in seq_len(max_components)) {
    components <- sapply(underscore_splits, function(x) if (length(x) >= i) x[i] else NA)
    components <- components[!is.na(components)]

    # Check for barcode pattern
    matches_bare <- sum(grepl(barcode_pattern, components))
    matches_suffix <- sum(grepl(barcode_with_suffix, components))

    if (matches_bare > 0 || matches_suffix > 0) {
      barcode_positions[[paste0("position_", i)]] <- list(
        bare_barcode = matches_bare,
        barcode_with_suffix = matches_suffix,
        total_checked = length(components)
      )
    }
  }

  result$barcode_positions <- barcode_positions

  return(result)
}

pattern_analysis <- analyze_cell_ids(cell_ids)

# Extract unique prefixes (before barcode)
underscore_splits <- strsplit(cell_ids, "_", fixed = TRUE)
n_components <- sapply(underscore_splits, length)
component_counts <- table(n_components)

# Extract prefixes for most common pattern
most_common_n <- as.integer(names(which.max(component_counts)))

# Get prefix (all components except last which is likely barcode)
prefixes <- sapply(underscore_splits, function(x) {
  if (length(x) >= 2) {
    paste(x[1:(length(x)-1)], collapse = "_")
  } else {
    x[1]
  }
})
prefix_counts <- sort(table(prefixes), decreasing = TRUE)

# Extract barcode component (last underscore-separated component)
barcodes <- sapply(underscore_splits, function(x) x[length(x)])

# Analyze barcode suffixes
barcode_suffixes <- sub("^[ACGT]+", "", barcodes)
suffix_counts <- table(barcode_suffixes)

# Barcode lengths (excluding suffix)
barcode_cores <- sub("-[0-9]+$", "", barcodes)
barcode_lengths <- nchar(barcode_cores)
barcode_length_dist <- table(barcode_lengths)

cell_id_data <- list(
  total_count = n_cells,

  patterns = list(
    list(
      pattern = paste0(most_common_n, " components separated by underscore"),
      example = cell_ids[1],
      count = as.integer(component_counts[as.character(most_common_n)]),
      fraction = round(100 * as.integer(component_counts[as.character(most_common_n)]) / n_cells, 2)
    )
  ),

  components = list(
    n_components_distribution = as.list(component_counts),

    prefix = list(
      unique_count = length(prefix_counts),
      values = as.list(prefix_counts)
    ),

    barcode = list(
      length_distribution = as.list(barcode_length_dist),
      suffix_pattern = as.list(suffix_counts)
    )
  ),

  pattern_analysis = pattern_analysis,

  examples = list(
    first_10 = head(cell_ids, 10),
    last_10 = tail(cell_ids, 10)
  )
)

message("  Unique prefixes: ", length(prefix_counts))
message("  Most common prefix: ", names(prefix_counts)[1], " (", prefix_counts[1], " cells)")

# =============================================================================
# Metadata Extraction
# =============================================================================

message("\n[4/5] Extracting metadata...")

meta <- obj@meta.data

metadata_result <- list(
  column_count = ncol(meta),
  row_count = nrow(meta),
  columns = list()
)

message("  Columns: ", ncol(meta))

for (col in colnames(meta)) {
  values <- meta[[col]]

  # Get dtype
  dtype <- class(values)[1]

  # Handle factors
  if (is.factor(values)) {
    values_for_table <- as.character(values)
  } else {
    values_for_table <- values
  }

  # Count non-null
  non_null <- sum(!is.na(values))
  null_count <- sum(is.na(values))
  coverage_pct <- round(100 * non_null / length(values), 2)

  # Get unique values with counts
  if (is.numeric(values) && length(unique(values)) > 50) {
    # For continuous numeric, provide summary instead of all values
    value_counts <- list(
      min = min(values, na.rm = TRUE),
      max = max(values, na.rm = TRUE),
      mean = round(mean(values, na.rm = TRUE), 4),
      median = median(values, na.rm = TRUE),
      n_unique = length(unique(values[!is.na(values)])),
      type = "continuous_numeric"
    )
  } else {
    # For categorical/discrete, list all values
    tbl <- table(values_for_table, useNA = "no")
    value_counts <- as.list(tbl)
  }

  metadata_result$columns[[col]] <- list(
    dtype = dtype,
    non_null = non_null,
    null = null_count,
    coverage_pct = coverage_pct,
    unique = length(unique(values[!is.na(values)])),
    values = value_counts
  )

  message("    - ", col, " (", dtype, "): ", length(unique(values[!is.na(values)])), " unique, ", coverage_pct, "% coverage")
}

# =============================================================================
# Assay and Reduction Info
# =============================================================================

message("\n[5/5] Extracting assay and reduction info...")

assay_info <- list()
if (inherits(obj, "Seurat")) {
  for (assay_name in Assays(obj)) {
    assay_obj <- obj[[assay_name]]
    assay_info[[assay_name]] <- list(
      n_features = nrow(assay_obj),
      n_cells = ncol(assay_obj),
      layers = names(assay_obj),
      feature_examples = list(
        first_10 = head(rownames(assay_obj), 10)
      )
    )
    message("  Assay '", assay_name, "': ", nrow(assay_obj), " features")
  }
}

reduction_info <- list()
if (inherits(obj, "Seurat")) {
  for (red_name in Reductions(obj)) {
    red_obj <- obj[[red_name]]
    reduction_info[[red_name]] <- list(
      dimensions = ncol(Embeddings(red_obj)),
      n_cells = nrow(Embeddings(red_obj))
    )
    message("  Reduction '", red_name, "': ", ncol(Embeddings(red_obj)), " dimensions")
  }
}

# =============================================================================
# Additional Observations
# =============================================================================

observations <- list()

# Check for sample/donor column candidates
sample_candidates <- c("orig.ident", "sample", "donor", "Sample_ID", "patient", "subject")
found_sample_cols <- intersect(sample_candidates, colnames(meta))
if (length(found_sample_cols) > 0) {
  observations <- c(observations, paste0("Potential sample ID columns found: ", paste(found_sample_cols, collapse = ", ")))
}

# Cell distribution by most likely sample column
if ("orig.ident" %in% colnames(meta)) {
  dist_by_orig <- table(meta$orig.ident)
  observations <- c(observations, paste0("Cell distribution by orig.ident: ", length(dist_by_orig), " unique values"))
}

# =============================================================================
# Compile Final Output
# =============================================================================

output <- list(
  study = "murrow",
  extraction_date = format(Sys.Date(), "%Y-%m-%d"),
  extraction_type = "complete",
  source_file = source_path,
  source_size_bytes = file_info$size,

  object_summary = object_summary,
  cell_ids = cell_id_data,
  metadata = metadata_result,
  assays = assay_info,
  reductions = reduction_info,
  observations = observations
)

# =============================================================================
# Write YAML
# =============================================================================

message("\n=== Writing output to YAML ===")

# Ensure output directory exists
output_dir <- dirname(output_path)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Write YAML
write_yaml(output, output_path)

message("Output written to: ", output_path)
message("\n=== EXTRACTION COMPLETE ===")
message("Total cells: ", format(n_cells, big.mark = ","))
message("Metadata columns: ", ncol(meta))
message("Unique prefixes (potential batches): ", length(prefix_counts))
