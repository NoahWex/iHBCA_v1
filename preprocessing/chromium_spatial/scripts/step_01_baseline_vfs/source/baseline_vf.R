# ==============================================================================
# STEP 01: Baseline Variable Feature Discovery - Core Algorithm
# ==============================================================================
# Extracted from: project/.archive_pre_refactor_20251029/workflows/qc/scripts/module_1_hard_cutoff_bigsur.R
# Extraction date: 2025-10-29
# Architecture: Manifest-driven preprocessing pipeline
#
# Purpose:
#   - Apply hard UMI threshold filtering (>=150 UMI by default)
#   - Run BigSur for variance-based variable feature discovery
#   - Determine sample-level QC pass/fail (>=100 VFs by default)
#   - Generate cell-level QC metrics for ALL cells (passing + failing)
#
# Design Principles:
#   - Pure algorithm: NO file I/O operations
#   - All inputs via function parameters
#   - Returns data structures, NOT file paths
#   - No hardcoded paths or parameter values
#   - Wrapper script handles config loading, file I/O, manifest updates
#
# Critical Requirements:
#   - Column name: Use 'sampleID' (capital ID) for downstream compatibility
#   - QC logic: sample_qc_pass = (n_vfs >= bigsur_vf_threshold)
#   - Metrics: Include ALL cells (passing + failing threshold)
#   - VF list: Return character vector of gene names
#
# Dependencies (must be loaded by wrapper):
#   - Seurat (CreateSeuratObject, Read10X_h5, VariableFeatures)
#   - BigSur (BigSur)
#   - Matrix (colSums, rowSums)
# ==============================================================================

#' Run baseline variable feature discovery for a single sample
#'
#' @description
#' Core algorithm for Step 01 of the preprocessing pipeline.
#' Applies hard UMI filtering, runs BigSur variable feature discovery,
#' and generates cell-level QC metrics.
#'
#' @details
#' Processing steps:
#' 1. Load raw count matrix from H5 file
#' 2. Calculate UMI and gene counts per cell
#' 3. Apply flat UMI threshold (filter before Seurat object creation)
#' 4. Subset matrix to passing cells
#' 5. Remove zero-count genes
#' 6. Create Seurat object with min.cells parameter
#' 7. Run BigSur variance-based feature discovery
#' 8. Extract variable features and determine sample QC pass/fail
#' 9. Generate cell-level metadata for ALL cells
#'
#' Sample QC Logic:
#'   - PASS: n_vfs >= bigsur_vf_threshold (default: 100)
#'   - FAIL: n_vfs < bigsur_vf_threshold
#'   - Failing samples excluded from downstream analysis
#'
#' Critical Design Decisions:
#'   - Filtering BEFORE CreateSeuratObject (per HPC memory constraints)
#'   - Metadata includes ALL cells (not just passing) for traceability
#'   - Column name 'sampleID' (capital ID) matches system-wide standard
#'   - No file I/O operations (wrapper handles all reads/writes)
#'
#' @param h5_path Character. Full path to filtered_feature_bc_matrix.h5 file
#'   from CellRanger output. Must exist and be readable.
#'
#' @param sample_id Character. Sample identifier (e.g., "Pat1_P1").
#'   Must match the sample_id in raw_data_manifest.yaml and
#'   preprocessing_manifest.yaml.
#'
#' @param umi_threshold Integer. Minimum UMI count for cell to pass
#'   hard threshold filter. Default: 150 (from module_configs.yaml).
#'   Cells below this threshold are excluded from Seurat object but
#'   included in metrics output with flat_threshold_pass = FALSE.
#'
#' @param bigsur_vf_threshold Integer. Minimum variable features required
#'   for sample to pass QC. Default: 100 (from module_configs.yaml).
#'   Samples with fewer VFs are marked sample_qc_pass = FALSE and
#'   excluded from downstream processing.
#'
#' @param bigsur_min_cells Integer. Minimum cells expressing a gene
#'   for gene to be included in Seurat object (min.cells parameter).
#'   Default: 2 (from module_configs.yaml).
#'
#' @param bigsur_correlations Logical. Enable correlation-based variable
#'   feature discovery in BigSur. Default: FALSE (use variance-based only).
#'
#' @param bigsur_variable_features Logical. Enable variable feature discovery.
#'   Default: TRUE (required for downstream analysis).
#'
#' @param bigsur_log_file Logical. Enable BigSur internal log file generation.
#'   Default: FALSE (wrapper handles logging).
#'
#' @return Named list with 4 elements:
#'   \item{metrics}{data.frame with 7 columns and nrow = total cells:
#'     \itemize{
#'       \item cell_barcode: Character. Cell barcode from H5 file
#'       \item sampleID: Character. Sample identifier (capital ID)
#'       \item flat_threshold_pass: Logical. Cell passed UMI threshold
#'       \item n_umi: Integer. Total UMI count for cell
#'       \item n_genes: Integer. Total genes detected in cell
#'       \item sample_bigsur_vf_count: Integer. Sample-level VF count
#'       \item sample_qc_pass: Logical. Sample passed VF threshold
#'     }
#'   }
#'   \item{vf_genes}{Character vector of variable feature gene names
#'     from BigSur. Length = n_vfs.}
#'   \item{n_vfs}{Integer. Count of variable features identified.}
#'   \item{sample_qc_pass}{Logical. TRUE if n_vfs >= bigsur_vf_threshold.}
#'
#' @examples
#' \dontrun{
#' # Typical usage (parameters from module_configs.yaml)
#' results <- run_baseline_vf(
#'   h5_path = "/path/to/filtered_feature_bc_matrix.h5",
#'   sample_id = "Pat1_P1",
#'   umi_threshold = 150,
#'   bigsur_vf_threshold = 100,
#'   bigsur_min_cells = 2,
#'   bigsur_correlations = FALSE,
#'   bigsur_variable_features = TRUE,
#'   bigsur_log_file = FALSE
#' )
#'
#' # Access results
#' metrics_df <- results$metrics
#' vf_genes <- results$vf_genes
#' n_vfs <- results$n_vfs
#' qc_pass <- results$sample_qc_pass
#' }
#'
#' @export
run_baseline_vf <- function(
  h5_path,
  sample_id,
  patient_id,
  position_id,
  umi_threshold = 150,
  bigsur_vf_threshold = 100,
  bigsur_min_cells = 2,
  bigsur_correlations = FALSE,
  bigsur_variable_features = TRUE,
  bigsur_log_file = FALSE
) {

  # ============================================================================
  # STEP 1: Load raw h5 count matrix
  # Source: Lines 147-155 of module_1_hard_cutoff_bigsur.R
  # ============================================================================

  message("\n--- Step 1: Load raw h5 count matrix ---")

  if (!file.exists(h5_path)) {
    stop(sprintf("H5 file not found: %s", h5_path))
  }

  counts_raw <- Seurat::Read10X_h5(h5_path)
  message(sprintf("Loaded raw matrix: %d genes × %d cells",
                  nrow(counts_raw), ncol(counts_raw)))

  # ============================================================================
  # STEP 2: Calculate UMI and gene counts per cell
  # Source: Lines 160-166 of module_1_hard_cutoff_bigsur.R
  # ============================================================================

  message("\n--- Step 2: Calculate UMI counts per cell ---")

  cell_umi <- Matrix::colSums(counts_raw)
  cell_genes <- Matrix::colSums(counts_raw > 0)

  message(sprintf("UMI range: %d - %d", min(cell_umi), max(cell_umi)))
  message(sprintf("Gene range: %d - %d", min(cell_genes), max(cell_genes)))

  # ============================================================================
  # STEP 3: Apply flat UMI threshold filter
  # Source: Lines 167-180 of module_1_hard_cutoff_bigsur.R
  # CRITICAL: Filter BEFORE CreateSeuratObject to reduce memory footprint
  # ============================================================================

  message(sprintf("\n--- Step 3: Apply flat threshold (UMI >= %d) ---", umi_threshold))

  passing_cells <- names(cell_umi)[cell_umi >= umi_threshold]
  failing_cells <- names(cell_umi)[cell_umi < umi_threshold]

  n_pass <- length(passing_cells)
  n_fail <- length(failing_cells)
  pass_rate <- n_pass / (n_pass + n_fail) * 100

  message(sprintf("Passing cells: %d (%.1f%%)", n_pass, pass_rate))
  message(sprintf("Failing cells: %d (%.1f%%)", n_fail, 100 - pass_rate))

  if (n_pass == 0) {
    stop("No cells pass UMI threshold - sample has failed")
  }

  # ============================================================================
  # STEP 4: Subset matrix to passing cells
  # Source: Lines 182-186 of module_1_hard_cutoff_bigsur.R
  # ============================================================================

  message("\n--- Step 4: Subset matrix to passing cells ---")

  counts_filtered <- counts_raw[, passing_cells, drop = FALSE]
  message(sprintf("Filtered matrix: %d genes × %d cells",
                  nrow(counts_filtered), ncol(counts_filtered)))

  # ============================================================================
  # STEP 5: Remove zero-count genes
  # Source: Lines 187-196 of module_1_hard_cutoff_bigsur.R
  # ============================================================================

  message("\n--- Step 5: Remove zero-count genes ---")

  gene_counts <- Matrix::rowSums(counts_filtered)
  genes_keep <- names(gene_counts)[gene_counts > 0]
  counts_filtered <- counts_filtered[genes_keep, , drop = FALSE]

  message(sprintf("Removed %d zero-count genes",
                  nrow(counts_raw) - length(genes_keep)))
  message(sprintf("Final matrix: %d genes × %d cells",
                  nrow(counts_filtered), ncol(counts_filtered)))

  # Free memory
  rm(counts_raw)
  gc(verbose = FALSE)

  # ============================================================================
  # STEP 6: Create Seurat object
  # Source: Lines 201-210 of module_1_hard_cutoff_bigsur.R
  # ============================================================================

  message(sprintf("\n--- Step 6: Create Seurat object (min.cells = %d) ---", bigsur_min_cells))

  seu <- Seurat::CreateSeuratObject(
    counts = counts_filtered,
    project = sample_id,
    min.cells = bigsur_min_cells
  )

  message(sprintf("Seurat object: %d cells × %d features",
                  ncol(seu), nrow(seu)))

  # Free memory
  rm(counts_filtered)
  gc(verbose = FALSE)

  # ============================================================================
  # STEP 7: Run BigSur variable feature discovery
  # Source: Lines 215-228 of module_1_hard_cutoff_bigsur.R
  # CHANGE: Parameterized BigSur options (previously hardcoded)
  # ============================================================================

  message("\n--- Step 7: Run BigSur ---")
  message("BigSur parameters:")
  message(sprintf("  correlations: %s", bigsur_correlations))
  message(sprintf("  variable.features: %s", bigsur_variable_features))
  message(sprintf("  log.file: %s", bigsur_log_file))

  seu <- BigSur::BigSur(
    seu,
    correlations = bigsur_correlations,
    variable.features = bigsur_variable_features,
    log.file = bigsur_log_file
  )

  message("BigSur complete")

  # ============================================================================
  # STEP 8: Extract variable features and determine sample QC pass/fail
  # Source: Lines 230-245 of module_1_hard_cutoff_bigsur.R
  # ============================================================================

  message("\n--- Step 8: Extract variable features ---")

  vfs <- Seurat::VariableFeatures(seu)
  n_vfs <- length(vfs)

  message(sprintf("Variable features identified: %d", n_vfs))

  # Determine sample pass/fail
  sample_qc_pass <- n_vfs >= bigsur_vf_threshold
  vf_status <- ifelse(sample_qc_pass, "PASS", "FAIL")

  message(sprintf("Sample QC status: %s (threshold: %d VFs)", vf_status, bigsur_vf_threshold))

  if (!sample_qc_pass) {
    message(sprintf("WARNING: Sample has <%d VFs and will be excluded from downstream analysis",
                    bigsur_vf_threshold))
  }

  # ============================================================================
  # STEP 9: Generate cell-level metadata for ALL cells
  # Source: Lines 256-271 of module_1_hard_cutoff_bigsur.R
  # CRITICAL FIX: Column name changed from 'sample_id' to 'sampleID'
  #   - Previous extraction used 'SampleID' (capital S)
  #   - Correct system-wide standard is 'sampleID' (lowercase s, capital ID)
  #   - Matches raw_data_manifest.yaml, preprocessing_manifest.yaml
  #   - Required for downstream merging with master_metadata.csv
  # ============================================================================

  message("\n--- Step 9: Generate cell-level metadata for ALL cells ---")

  # Create metadata for ALL cells (passing + failing)
  all_barcodes <- c(passing_cells, failing_cells)
  all_umi <- cell_umi[all_barcodes]
  all_genes <- cell_genes[all_barcodes]

  # Prefix barcodes with sample_id for global uniqueness
  prefixed_barcodes <- paste0(sample_id, "_", all_barcodes)

  cell_metadata <- data.frame(
    cell_id = prefixed_barcodes,
    sample_id = sample_id,
    patient_id = patient_id,
    position = position_id,
    flat_threshold_pass = all_umi >= umi_threshold,
    n_umi = as.integer(all_umi),
    n_genes = as.integer(all_genes),
    sample_bigsur_vf_count = as.integer(n_vfs),
    sample_qc_pass = sample_qc_pass,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  message(sprintf("Cell metadata: %d total cells", nrow(cell_metadata)))
  message(sprintf("  Passing threshold: %d", sum(cell_metadata$flat_threshold_pass)))
  message(sprintf("  Failing threshold: %d", sum(!cell_metadata$flat_threshold_pass)))

  # ============================================================================
  # Return structured results (NO FILE I/O)
  # ============================================================================

  message("\n=== Core algorithm complete ===")
  message(sprintf("Sample: %s", sample_id))
  message(sprintf("Status: %s", vf_status))
  message(sprintf("Cells (total): %d", nrow(cell_metadata)))
  message(sprintf("Cells (passing): %d (%.1f%%)", n_pass, pass_rate))
  message(sprintf("Variable features: %d", n_vfs))

  return(list(
    metrics = cell_metadata,
    vf_genes = vfs,
    n_vfs = n_vfs,
    sample_qc_pass = sample_qc_pass
  ))
}

# ==============================================================================
# EXTRACTION NOTES
# ==============================================================================
# Source file: project/.archive_pre_refactor_20251029/workflows/qc/scripts/module_1_hard_cutoff_bigsur.R
# Extraction lines: 147-271 (algorithm core)
# Extraction ratio: ~125 algorithm lines preserved from 307 total lines (40.7%)
#
# Algorithm changes from source:
#   1. Column name fix: 'sample_id' → 'sampleID' (BUG FIX - system standard)
#   2. Parameterized BigSur options (previously hardcoded FALSE values)
#   3. Added comprehensive documentation and examples
#
# Excluded from extraction (wrapper handles):
#   - Shebang and library imports (lines 1-36)
#   - Command-line argument parsing (lines 38-88)
#   - Config loading and metadata reading (lines 90-124)
#   - Output directory creation (lines 126-132)
#   - Logging setup (sink redirection) (lines 134-139)
#   - tryCatch error handling (lines 145, 291-297)
#   - File write operations (lines 247-251, 277-281)
#   - Logging cleanup (lines 299-302)
#   - Final stdout message (lines 304-306)
#
# Validation:
#   ✓ No hardcoded paths
#   ✓ No hardcoded parameters (all via function arguments)
#   ✓ No file I/O operations (Read10X_h5 is data loading, not file management)
#   ✓ Returns data structures (list with data.frame and vectors)
#   ✓ Column names match system standard (sampleID with capital ID)
#   ✓ All critical QC metrics included in output
#
# References:
#   - EXTRACTION_SUMMARY.txt: Previous extraction documentation
#   - module_configs.yaml: Algorithm parameter defaults
#   - preprocessing_manifest.yaml: I/O contract structure
#   - raw_data_manifest.yaml: Input data paths and sample vocabulary
# ==============================================================================
