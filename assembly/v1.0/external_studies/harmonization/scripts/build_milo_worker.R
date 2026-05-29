#!/usr/bin/env Rscript
# =============================================================================
# build_milo_worker.R
# =============================================================================
# Phase 2 of parallel Milo build pipeline.
# Unified Milo builder for both native and joint embedding modes.
# Uses BiocParallel for multi-core k-NN graph construction.
#
# Expects Phase 1 outputs in outputs/study_objects/{study}/:
#   seurat.rds, embedding_joint.csv (for joint mode)
#
# Usage:
#   Rscript build_milo_worker.R --study gray --mode native --cpus 8
#   Rscript build_milo_worker.R --study gray --mode joint --cpus 8
#   Rscript build_milo_worker.R --study gray --mode native --dry-run
#
# Outputs:
#   native mode: milo.rds, build_summary.yaml
#   joint mode:  milo_joint.rds, build_joint_summary.yaml
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(miloR)
  library(SingleCellExperiment)
  library(BiocParallel)
  library(argparse)
  library(yaml)
  library(data.table)
})

# Source config from build_study_objects.R for STUDY_CONFIG
SCRIPT_DIR <- dirname(sub("--file=", "",
  commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))]))
if (length(SCRIPT_DIR) == 0) SCRIPT_DIR <- "."
source(file.path(SCRIPT_DIR, "build_study_objects.R"))

# Joint embedding dims (fixed at 100 for scVI)
JOINT_DIMS <- 100L

# =============================================================================
# Milo Build Functions
# =============================================================================

# Why did the k-NN graph go to therapy?
# Because it had too many close neighbors and no personal space.

build_milo_with_bpparam <- function(sce, embedding_key, dims, k, prop, cpus) {
  bp <- if (cpus > 1L) {
    message("  Registering MulticoreParam with ", cpus, " workers")
    MulticoreParam(workers = cpus)
  } else {
    SerialParam()
  }

  message("  Building k-NN graph (k=", k, ", d=", dims,
          ", embedding=", embedding_key, ")...")
  milo <- Milo(sce)
  # BPPARAM must be passed explicitly — global register() doesn't propagate
  milo <- buildGraph(milo, k = k, d = dims, reduced.dim = embedding_key,
                     BPPARAM = bp)

  message("  Computing neighborhoods (prop=", prop, ")...")
  # makeNhoods is serial — no BPPARAM arg
  milo <- makeNhoods(milo, prop = prop, k = k, d = dims,
                     refined = TRUE, reduced_dims = embedding_key)

  n_nhoods <- ncol(nhoods(milo))
  message("  Created ", n_nhoods, " neighborhoods")

  message("  Computing neighborhood distances...")
  # calcNhoodDistance is serial — no BPPARAM arg
  milo <- calcNhoodDistance(milo, d = dims, reduced.dim = embedding_key)

  return(milo)
}

# =============================================================================
# Main Pipeline
# =============================================================================

build_milo_from_inputs <- function(study, mode, cpus = 1L, dry_run = FALSE) {
  if (!study %in% names(STUDY_CONFIG)) {
    stop("Unknown study: ", study)
  }
  if (!mode %in% c("native", "joint")) {
    stop("Mode must be 'native' or 'joint', got: ", mode)
  }

  config <- STUDY_CONFIG[[study]]
  study_dir <- file.path(DATA_PREP_PATH, "outputs", "study_objects", study)

  message("\n", strrep("=", 60))
  message("Building Milo: ", toupper(study), " [", mode, "]",
          if (dry_run) " (DRY RUN)" else "")
  message("  Input: ", study_dir)
  message("  CPUs: ", cpus)
  message(strrep("=", 60))

  # --- Load Seurat object ---
  seurat_file <- file.path(study_dir, "seurat.rds")
  if (!file.exists(seurat_file)) {
    stop("seurat.rds not found: ", seurat_file,
         "\n  Run prepare_study.R --study ", study, " first")
  }

  message("\n[1] Loading Seurat object...")
  seurat_obj <- readRDS(seurat_file)
  message("  Loaded: ", ncol(seurat_obj), " cells x ", nrow(seurat_obj), " genes")

  # --- Mode-specific setup ---
  if (mode == "native") {
    message("\n[2] Validating native embedding...")
    sce <- as.SingleCellExperiment(seurat_obj)
    embedding_key <- config$embedding_key
    dims <- config$embedding_dims

    # Resolve case-insensitive embedding name
    available <- reducedDimNames(sce)
    if (!embedding_key %in% available) {
      matches <- available[tolower(available) == tolower(embedding_key)]
      if (length(matches) == 0) {
        stop("Embedding '", embedding_key, "' not found. Available: ",
             paste(available, collapse = ", "))
      }
      embedding_key <- matches[1]
    }

    actual_dims <- min(dims, ncol(reducedDim(sce, embedding_key)))
    message("  Embedding: ", embedding_key, " (", actual_dims, " dims)")
    message("  Cells: ", ncol(sce))

    if (dry_run) {
      message("\n[DRY RUN] Validation passed. Would build native Milo with:")
      message("  k=", MILO_K, ", d=", actual_dims, ", prop=", MILO_PROP)
      message("  cpus=", cpus)
      message("  Output: milo.rds, build_summary.yaml")
      message("\nDRY RUN complete — exiting 0")
      return(invisible(NULL))
    }

    message("\n[3] Building Milo...")
    milo <- build_milo_with_bpparam(sce, embedding_key, actual_dims,
                                    MILO_K, MILO_PROP, cpus)

    # Save
    milo_file <- file.path(study_dir, "milo.rds")
    saveRDS(milo, milo_file)
    message("  Saved: ", basename(milo_file))

    # Merge with prepare_summary if it exists
    prepare_summary_file <- file.path(study_dir, "prepare_summary.yaml")
    prepare_info <- if (file.exists(prepare_summary_file)) {
      read_yaml(prepare_summary_file)
    } else {
      list()
    }

    summary <- list(
      study = study,
      created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      source_file = config$file,
      n_cells = ncol(seurat_obj),
      n_genes = nrow(seurat_obj),
      n_neighborhoods = ncol(nhoods(milo)),
      metadata_match_rate = if (!is.null(prepare_info$metadata_match_rate)) {
        prepare_info$metadata_match_rate
      } else {
        round(seurat_obj@misc$metadata_match_rate %||% NA_real_, 4)
      },
      embedding_native = list(key = config$embedding_key,
                              dims = config$embedding_dims),
      embedding_joint = if (!is.null(prepare_info$embedding_joint)) {
        prepare_info$embedding_joint
      } else {
        list(available = NA, dims = NA)
      },
      milo_params = list(k = MILO_K, prop = MILO_PROP),
      build_cpus = cpus,
      outputs = list(
        metadata = "metadata.csv",
        embedding_native = "embedding_native.csv",
        embedding_joint = if (!is.null(prepare_info$outputs$embedding_joint)) {
          prepare_info$outputs$embedding_joint
        } else {
          NA
        },
        seurat = "seurat.rds",
        milo = "milo.rds"
      )
    )

    summary_file <- file.path(study_dir, "build_summary.yaml")
    write_yaml(summary, summary_file)
    message("  Saved: ", basename(summary_file))

  } else {
    # --- Joint mode ---
    message("\n[2] Loading joint embedding...")
    joint_file <- file.path(study_dir, "embedding_joint.csv")
    if (!file.exists(joint_file)) {
      stop("embedding_joint.csv not found: ", joint_file,
           "\n  Run prepare_study.R --study ", study, " first (without --skip-joint)")
    }

    emb <- fread(joint_file)
    cell_ids <- emb[[1]]
    joint_mat <- as.matrix(emb[, -1, with = FALSE])
    rownames(joint_mat) <- cell_ids
    message("  Joint embedding: ", nrow(joint_mat), " cells x ",
            ncol(joint_mat), " dims")

    # Cell matching
    common_cells <- intersect(colnames(seurat_obj), rownames(joint_mat))
    n_common <- length(common_cells)
    message("  Cell matching: Seurat=", ncol(seurat_obj),
            ", Joint=", nrow(joint_mat), ", Common=", n_common)

    if (n_common == 0) {
      stop("No common cells between Seurat object and joint embedding")
    }

    if (n_common < ncol(seurat_obj)) {
      message("  Subsetting to ", n_common, " common cells")
      seurat_obj <- seurat_obj[, common_cells]
    }

    joint_mat <- joint_mat[colnames(seurat_obj), , drop = FALSE]

    if (dry_run) {
      message("\n[DRY RUN] Validation passed. Would build joint Milo with:")
      message("  k=", MILO_K, ", d=", JOINT_DIMS, ", prop=", MILO_PROP)
      message("  cpus=", cpus)
      message("  Cells: ", n_common)
      message("  Output: milo_joint.rds, build_joint_summary.yaml")
      message("\nDRY RUN complete — exiting 0")
      return(invisible(NULL))
    }

    message("\n[3] Building Milo on joint embedding...")
    sce <- as.SingleCellExperiment(seurat_obj)
    reducedDim(sce, "scVI_joint") <- joint_mat

    milo <- build_milo_with_bpparam(sce, "scVI_joint", JOINT_DIMS,
                                    MILO_K, MILO_PROP, cpus)

    # Save
    milo_file <- file.path(study_dir, "milo_joint.rds")
    saveRDS(milo, milo_file)
    message("  Saved: ", basename(milo_file))

    summary <- list(
      study = study,
      created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      type = "joint",
      n_cells = ncol(milo),
      n_neighborhoods = ncol(nhoods(milo)),
      embedding = list(name = "scVI_joint", dims = JOINT_DIMS,
                       source = "embedding_joint.csv"),
      milo_params = list(k = MILO_K, prop = MILO_PROP),
      build_cpus = cpus
    )

    summary_file <- file.path(study_dir, "build_joint_summary.yaml")
    write_yaml(summary, summary_file)
    message("  Saved: ", basename(summary_file))
  }

  message("\n", strrep("=", 60))
  message("SUCCESS: ", study, " [", mode, "] Milo built")
  message("  Cells: ", ncol(milo))
  message("  Neighborhoods: ", ncol(nhoods(milo)))
  message(strrep("=", 60), "\n")

  return(invisible(summary))
}

# =============================================================================
# CLI Entry Point
# =============================================================================

.is_main_script_worker <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  script_arg <- args[grep("--file=", args)]
  if (length(script_arg) == 0) return(FALSE)
  basename(sub("--file=", "", script_arg)) == "build_milo_worker.R"
}

if (!interactive() && .is_main_script_worker()) {
  parser <- ArgumentParser(description = "Build Milo object (Phase 2 — parallel)")
  parser$add_argument("--study", type = "character", required = TRUE,
                      help = "Study name")
  parser$add_argument("--mode", type = "character", required = TRUE,
                      choices = c("native", "joint"),
                      help = "Embedding mode: native or joint")
  parser$add_argument("--cpus", type = "integer", default = 1L,
                      help = "Number of CPUs for BiocParallel (default: 1)")
  parser$add_argument("--dry-run", action = "store_true",
                      help = "Validate inputs and print plan without building")

  args <- parser$parse_args()
  build_milo_from_inputs(args$study, args$mode, cpus = args$cpus,
                         dry_run = args$dry_run)
}
