#!/usr/bin/env Rscript
# =============================================================================
# build_joint_milo.R
# =============================================================================
# Build Milo object on joint scVI embedding (100d) for integration assessment.
# This script builds milo_joint.rds from existing seurat.rds + embedding_joint.csv.
#
# Usage:
#   Rscript build_joint_milo.R --study gray
#   Rscript build_joint_milo.R --all
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(miloR)
  library(SingleCellExperiment)
  library(argparse)
  library(yaml)
  library(data.table)
})

# =============================================================================
# Configuration
# =============================================================================

BASE_PATH <- "${SOURCE_IHBCAV1_EXTERNAL_STUDIES%/}"
DATA_PREP_PATH <- file.path(BASE_PATH, "harmonization")

STUDIES <- c("gray", "kumar", "murrow", "nee", "twigger",
             "pal_norm_epi", "pal_norm_total", "pal_norm_b1", "reed")

MILO_K <- 30
MILO_PROP <- 0.1
JOINT_DIMS <- 100

# =============================================================================
# Helper Functions
# =============================================================================

load_joint_embedding <- function(study_dir) {
  emb_file <- file.path(study_dir, "embedding_joint.csv")
  if (!file.exists(emb_file)) {
    stop("Joint embedding not found: ", emb_file)
  }
  message("  Loading joint embedding...")
  emb <- fread(emb_file)
  cell_ids <- emb[[1]]
  emb_mat <- as.matrix(emb[, -1, with = FALSE])
  rownames(emb_mat) <- cell_ids
  message("  Joint embedding: ", nrow(emb_mat), " cells x ", ncol(emb_mat), " dims")
  return(emb_mat)
}

build_joint_milo <- function(seurat_obj, joint_embedding) {
  message("  Building Milo on joint embedding...")
  
  common_cells <- intersect(colnames(seurat_obj), rownames(joint_embedding))
  n_seurat <- ncol(seurat_obj)
  n_joint <- nrow(joint_embedding)
  n_common <- length(common_cells)
  
  message("  Cell matching: Seurat=", n_seurat, ", Joint=", n_joint, ", Common=", n_common)
  
  if (n_common == 0) {
    stop("No common cells between Seurat object and joint embedding")
  }
  
  if (n_common < n_seurat) {
    message("  Subsetting to ", n_common, " common cells")
    seurat_obj <- seurat_obj[, common_cells]
  }
  
  joint_embedding <- joint_embedding[colnames(seurat_obj), , drop = FALSE]
  
  sce <- as.SingleCellExperiment(seurat_obj)
  reducedDim(sce, "scVI_joint") <- joint_embedding
  milo <- Milo(sce)
  
  message("  Building k-NN graph (k=", MILO_K, ", d=", JOINT_DIMS, ")...")
  milo <- buildGraph(milo, k = MILO_K, d = JOINT_DIMS, reduced.dim = "scVI_joint")
  
  message("  Computing neighborhoods (prop=", MILO_PROP, ")...")
  milo <- makeNhoods(milo, prop = MILO_PROP, k = MILO_K, d = JOINT_DIMS,
                     refined = TRUE, reduced_dims = "scVI_joint")
  
  n_nhoods <- ncol(nhoods(milo))
  message("  Created ", n_nhoods, " neighborhoods")
  
  message("  Computing neighborhood distances...")
  milo <- calcNhoodDistance(milo, d = JOINT_DIMS, reduced.dim = "scVI_joint")
  
  return(milo)
}

save_outputs <- function(study, milo, study_dir) {
  milo_file <- file.path(study_dir, "milo_joint.rds")
  saveRDS(milo, milo_file)
  message("  Saved: ", basename(milo_file))
  
  summary <- list(
    study = study,
    created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    type = "joint",
    n_cells = ncol(milo),
    n_neighborhoods = ncol(nhoods(milo)),
    embedding = list(name = "scVI_joint", dims = JOINT_DIMS, source = "embedding_joint.csv"),
    milo_params = list(k = MILO_K, prop = MILO_PROP)
  )
  
  summary_file <- file.path(study_dir, "build_joint_summary.yaml")
  write_yaml(summary, summary_file)
  message("  Saved: ", basename(summary_file))
  return(summary)
}

# =============================================================================
# Main Pipeline
# =============================================================================

build_joint_milo_for_study <- function(study) {
  study_dir <- file.path(DATA_PREP_PATH, "outputs", "study_objects", study)
  seurat_file <- file.path(study_dir, "seurat.rds")
  joint_file <- file.path(study_dir, "embedding_joint.csv")
  
  if (!file.exists(seurat_file)) {
    stop("Seurat object not found: ", seurat_file)
  }
  if (!file.exists(joint_file)) {
    stop("Joint embedding not found: ", joint_file)
  }
  
  message("\n", strrep("=", 60))
  message("Building Joint Milo for: ", toupper(study))
  message("  Input: ", study_dir)
  message(strrep("=", 60))
  
  message("\n[1/3] Loading Seurat object...")
  seurat_obj <- readRDS(seurat_file)
  message("  Loaded: ", ncol(seurat_obj), " cells")
  
  message("\n[2/3] Loading joint embedding...")
  joint_embedding <- load_joint_embedding(study_dir)
  
  message("\n[3/3] Building Milo object...")
  milo <- build_joint_milo(seurat_obj, joint_embedding)
  
  message("\n[SAVE] Saving outputs...")
  summary <- save_outputs(study, milo, study_dir)
  
  message("\n", strrep("=", 60))
  message("SUCCESS: ", study, " joint Milo built")
  message("  Cells: ", summary$n_cells)
  message("  Neighborhoods: ", summary$n_neighborhoods)
  message("  Output: ", file.path(study_dir, "milo_joint.rds"))
  message(strrep("=", 60), "\n")
  
  return(invisible(summary))
}

# =============================================================================
# CLI Entry Point
# =============================================================================

.is_main_script <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  script_arg <- args[grep("--file=", args)]
  if (length(script_arg) == 0) return(FALSE)
  basename(sub("--file=", "", script_arg)) == "build_joint_milo.R"
}

if (!interactive() && .is_main_script()) {
  parser <- ArgumentParser(description = "Build Milo on joint scVI embedding")
  parser$add_argument("--study", type = "character", help = "Study name")
  parser$add_argument("--all", action = "store_true", help = "Build for all studies")
  
  args <- parser$parse_args()
  
  if (args$all) {
    for (study in STUDIES) {
      tryCatch({
        build_joint_milo_for_study(study)
      }, error = function(e) {
        message("\nERROR building ", study, ": ", e$message, "\n")
      })
    }
  } else if (!is.null(args$study)) {
    if (!args$study %in% STUDIES) {
      stop("Unknown study: ", args$study)
    }
    build_joint_milo_for_study(args$study)
  } else {
    stop("Must specify --study or --all")
  }
}
