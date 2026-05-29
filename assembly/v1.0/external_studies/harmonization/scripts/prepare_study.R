#!/usr/bin/env Rscript
# =============================================================================
# prepare_study.R
# =============================================================================
# Phase 1 of parallel Milo build pipeline.
# Runs steps 1-6 only: load RDS → rename cells → attach metadata → extract
# embeddings → save. No Milo build — that's handled by build_milo_worker.R.
#
# Sources build_study_objects.R for STUDY_CONFIG and helper functions.
#
# Usage:
#   Rscript prepare_study.R --study gray
#   Rscript prepare_study.R --study gray --skip-joint
#
# Outputs (in outputs/study_objects/{study}/):
#   seurat.rds, metadata.csv, embedding_native.csv, embedding_joint.csv,
#   prepare_summary.yaml
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(argparse)
  library(dplyr)
  library(yaml)
  library(data.table)
})

# Source config and helpers from build_study_objects.R
# (guarded by .is_main_script() so CLI block won't execute)
SCRIPT_DIR <- dirname(sub("--file=", "",
  commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))]))
if (length(SCRIPT_DIR) == 0) SCRIPT_DIR <- "."
source(file.path(SCRIPT_DIR, "build_study_objects.R"))

# =============================================================================
# Main Pipeline — Steps 1-6 only
# =============================================================================

# Milo is expensive. Prep is cheap. That's the whole point.
# -- ancient HPC proverb (ca. 2026)

prepare_study <- function(study, skip_joint = FALSE) {
  if (!study %in% names(STUDY_CONFIG)) {
    stop("Unknown study: ", study, ". Available: ",
         paste(names(STUDY_CONFIG), collapse = ", "))
  }

  config <- STUDY_CONFIG[[study]]
  output_dir <- file.path(DATA_PREP_PATH, "outputs", "study_objects", study)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  message("\n", strrep("=", 60))
  message("Preparing study: ", toupper(study), " (Phase 1 — no Milo)")
  message("  Source: ", basename(config$file))
  message("  Output: ", output_dir)
  message(strrep("=", 60))

  # [1] Load cell mapping
  message("\n[1/6] Loading cell mapping...")
  cell_mapping <- load_cell_mapping(study)

  # [2] Load donor metadata
  message("\n[2/6] Loading donor metadata...")
  donor_metadata <- load_donor_metadata()
  study_key <- tolower(sub("_norm.*", "", study))
  donor_metadata <- donor_metadata[tolower(donor_metadata$study) == study_key, ]
  message("  Filtered to ", nrow(donor_metadata), " donors for ", study,
          " (study_key: ", study_key, ")")

  # [3] Load component RDS
  message("\n[3/6] Loading component RDS...")
  obj <- load_rds_object(config$file)

  # [4] Rename cells to iHBCA format
  message("\n[4/6] Renaming cells to iHBCA format...")
  obj <- rename_cells_to_ihbca(obj, cell_mapping, study)

  # [5] Attach harmonized metadata
  message("\n[5/6] Attaching harmonized metadata...")
  obj <- attach_metadata(obj, cell_mapping, donor_metadata, config$study_prefix)

  # [6] Rename reduction slot to canonical name if needed, then extract
  message("\n[6/6] Extracting embeddings...")
  embedding_key <- config$embedding_key
  available <- Reductions(obj)
  if (!embedding_key %in% available) {
    base_slot <- sub("^[a-z]+_", "", embedding_key)
    if (base_slot != embedding_key && base_slot %in% available) {
      message("  Renaming reduction '", base_slot, "' -> '", embedding_key, "'")
      obj[[embedding_key]] <- obj[[base_slot]]
      obj[[base_slot]] <- NULL
    }
  }
  embedding_native <- extract_embedding(obj, config$embedding_key,
                                        config$embedding_dims)

  embedding_joint <- NULL
  if (!skip_joint) {
    embedding_joint <- load_joint_embedding(colnames(obj))
    if (!is.null(embedding_joint)) {
      common_cells <- intersect(colnames(obj), rownames(embedding_joint))
      if (length(common_cells) < ncol(obj)) {
        message("  Warning: ", ncol(obj) - length(common_cells),
                " cells missing from joint embedding")
      }
      embedding_joint <- embedding_joint[common_cells, , drop = FALSE]
    }
  } else {
    message("  Skipping joint embedding (--skip-joint)")
  }

  cell_metadata <- obj@meta.data
  cell_metadata$cell_id <- rownames(cell_metadata)

  # --- Save outputs ---
  message("\n[SAVE] Saving outputs...")

  meta_file <- file.path(output_dir, "metadata.csv")
  write.csv(cell_metadata, meta_file, row.names = FALSE)
  message("  Saved: ", basename(meta_file))

  native_file <- file.path(output_dir, "embedding_native.csv")
  write.csv(embedding_native, native_file)
  message("  Saved: ", basename(native_file))

  if (!is.null(embedding_joint)) {
    joint_file <- file.path(output_dir, "embedding_joint.csv")
    write.csv(embedding_joint, joint_file)
    message("  Saved: ", basename(joint_file))
  }

  seurat_file <- file.path(output_dir, "seurat.rds")
  saveRDS(obj, seurat_file)
  message("  Saved: ", basename(seurat_file))

  # Summary (no milo_params or n_neighborhoods — that's Phase 2)
  summary <- list(
    study = study,
    phase = "prepare",
    created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    source_file = config$file,
    n_cells = ncol(obj),
    n_genes = nrow(obj),
    metadata_match_rate = round(obj@misc$metadata_match_rate, 4),
    embedding_native = list(key = config$embedding_key,
                            dims = config$embedding_dims),
    embedding_joint = list(
      available = !is.null(embedding_joint),
      dims = if (!is.null(embedding_joint)) ncol(embedding_joint) else NA
    ),
    outputs = list(
      metadata = "metadata.csv",
      embedding_native = "embedding_native.csv",
      embedding_joint = if (!is.null(embedding_joint)) "embedding_joint.csv" else NA,
      seurat = "seurat.rds"
    )
  )

  summary_file <- file.path(output_dir, "prepare_summary.yaml")
  write_yaml(summary, summary_file)
  message("  Saved: ", basename(summary_file))

  message("\n", strrep("=", 60))
  message("SUCCESS: ", study, " prepared (Phase 1)")
  message("  Cells: ", summary$n_cells)
  message("  Metadata match: ", round(100 * summary$metadata_match_rate, 1), "%")
  message("  Joint embedding: ",
          if (summary$embedding_joint$available) "YES" else "NO")
  message("  Output: ", output_dir)
  message("  Next: run build_milo_worker.R --study ", study, " --mode native|joint")
  message(strrep("=", 60), "\n")

  return(invisible(summary))
}

# =============================================================================
# CLI Entry Point
# =============================================================================

.is_main_script_prepare <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  script_arg <- args[grep("--file=", args)]
  if (length(script_arg) == 0) return(FALSE)
  basename(sub("--file=", "", script_arg)) == "prepare_study.R"
}

if (!interactive() && .is_main_script_prepare()) {
  parser <- ArgumentParser(description = "Prepare study objects (Phase 1 — no Milo)")
  parser$add_argument("--study", type = "character", required = TRUE,
                      help = "Study name (gray, kumar, murrow, nee, twigger, pal_norm_*, reed)")
  parser$add_argument("--skip-joint", action = "store_true",
                      help = "Skip joint embedding extraction")

  args <- parser$parse_args()
  prepare_study(args$study, skip_joint = args$skip_joint)
}
