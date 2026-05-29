#!/usr/bin/env Rscript
# =============================================================================
# build_study_objects.R
# =============================================================================
# Build complete study objects for Phase B analysis.
# Produces: metadata.csv, embedding_native.csv, embedding_joint.csv,
#           seurat.rds, milo.rds
#
# Usage:
#   Rscript build_study_objects.R --study gray
#   Rscript build_study_objects.R --study gray --skip-joint  # Skip joint embedding
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(miloR)
  library(SingleCellExperiment)
  library(argparse)
  library(dplyr)
  library(yaml)
  library(data.table)
})

# =============================================================================
# Configuration
# =============================================================================

BASE_PATH <- "${SOURCE_IHBCAV1_EXTERNAL_STUDIES%/}"
DATA_PREP_PATH <- file.path(BASE_PATH, "harmonization")

STUDY_CONFIG <- list(
  gray = list(
    file = "${SOURCE_COMPONENT_STUDIES}/gray.rds",
    embedding_key = "pca",
    embedding_dims = 50,
    study_prefix = "Gray"
  ),
  kumar = list(
    file = "${SOURCE_COMPONENT_STUDIES}/kumar.rds",
    embedding_key = "cca_pca",
    embedding_dims = 30,
    study_prefix = "Kumar"
  ),
  murrow = list(
    file = "${SOURCE_COMPONENT_STUDIES}/murrow.rds",
    embedding_key = "cca_pca",
    embedding_dims = 30,
    study_prefix = "Murrow"
  ),
  nee = list(
    file = "${PROJECT_SPATIAL_HBCA}/project/ReferenceDatasets/iHBCA/component_studies/nee.rds",
    embedding_key = "cca_pca",
    embedding_dims = 30,
    study_prefix = "Nee"
  ),
  twigger = list(
    file = "${SOURCE_COMPONENT_STUDIES}/twigger_l1_fastmnn.rds",
    embedding_key = "fastMNN",
    embedding_dims = 50,
    study_prefix = "Twigger"
  ),
  pal_norm_epi = list(
    file = "${SOURCE_PAL_ORIGINAL%/}/SeuratObject_NormEpi.rds",
    embedding_key = "cca_pca",
    embedding_dims = 30,
    study_prefix = "Pal"
  ),
  pal_norm_total = list(
    file = "${SOURCE_PAL_ORIGINAL%/}/SeuratObject_NormTotal.rds",
    embedding_key = "cca_pca",
    embedding_dims = 30,
    study_prefix = "Pal"
  ),
  pal_norm_b1 = list(
    file = "${SOURCE_PAL_ORIGINAL%/}/SeuratObject_NormB1Total.rds",
    embedding_key = "cca_pca",
    embedding_dims = 30,
    study_prefix = "Pal"
  ),
  reed = list(
    file = "${SOURCE_COMPONENT_STUDIES}/reed.rds",
    embedding_key = "scvi",
    embedding_dims = 20,
    study_prefix = "Reed"
  )
)

# Joint embedding paths
JOINT_EMBEDDING_PATH <- "${SOURCE_AUTHOR_SHARE}/X_scVI100.csv"
CELL_ANNOTATIONS_PATH <- "${SOURCE_AUTHOR_SHARE}/ihbca_level1.5_annotations.csv"

# iHBCA cell inventory (cell_id -> patient_id mapping)
CELL_INVENTORY_PATH <- file.path(DATA_PREP_PATH, "outputs", "ihbca_reference", "ihbca_cell_inventory.csv")

# Milo parameters
MILO_K <- 30
MILO_PROP <- 0.1

# =============================================================================
# Helper Functions
# =============================================================================

load_cell_mapping <- function(study) {
  mapping_file <- file.path(DATA_PREP_PATH, "studies", study, "outputs",
                            paste0(study, "_cells.csv"))
  if (!file.exists(mapping_file)) {
    stop("Cell mapping not found: ", mapping_file)
  }
  mapping <- fread(mapping_file, stringsAsFactors = FALSE)
  message("  Loaded ", nrow(mapping), " cell mappings from ", basename(mapping_file))
  return(as.data.frame(mapping))
}

load_donor_metadata <- function() {
  meta_file <- file.path(DATA_PREP_PATH, "outputs", "harmonized_metadata",
                         "harmonized_donor_metadata.csv")
  if (!file.exists(meta_file)) {
    stop("Harmonized donor metadata not found: ", meta_file)
  }
  # Use fread for consistency with inventory (both strip whitespace)
  metadata <- fread(meta_file, stringsAsFactors = FALSE)
  message("  Loaded ", nrow(metadata), " donors from harmonized metadata")
  return(as.data.frame(metadata))
}

load_rds_object <- function(path) {
  message("  Loading RDS object from: ", basename(path))
  obj <- readRDS(path)
  if (inherits(obj, "Seurat")) {
    tryCatch({
      obj <- UpdateSeuratObject(obj)
      message("  Updated to Seurat v5")
    }, error = function(e) {
      message("  Already compatible with current Seurat version")
    })
  }
  message("  Cells: ", ncol(obj), ", Genes: ", nrow(obj))
  return(obj)
}

extract_embedding <- function(obj, embedding_key, embedding_dims) {
  available <- Reductions(obj)
  if (embedding_key %in% available) {
    actual_key <- embedding_key
  } else {
    # Canonical name not found — check for base slot to rename
    # e.g. cca_pca -> look for pca; scvi_pca -> look for pca
    base_slot <- sub("^[a-z]+_", "", embedding_key)
    if (base_slot != embedding_key && base_slot %in% available) {
      message("  Renaming reduction '", base_slot, "' -> '", embedding_key, "'")
      obj[[embedding_key]] <- obj[[base_slot]]
      obj[[base_slot]] <- NULL
      actual_key <- embedding_key
    } else {
      matches <- available[tolower(available) == tolower(embedding_key)]
      if (length(matches) == 0) {
        stop("Embedding '", embedding_key, "' not found. Available: ",
             paste(available, collapse = ", "))
      }
      actual_key <- matches[1]
    }
  }
  emb <- Embeddings(obj, reduction = actual_key)
  if (ncol(emb) >= embedding_dims) {
    emb <- emb[, 1:embedding_dims]
  }
  message("  Extracted embedding '", actual_key, "': ", nrow(emb), " cells x ",
          ncol(emb), " dims")
  return(emb)
}

load_joint_embedding <- function(cell_ids) {
  message("  Loading joint embedding (X_scVI100.csv)...")
  if (!file.exists(JOINT_EMBEDDING_PATH)) {
    warning("Joint embedding file not found: ", JOINT_EMBEDDING_PATH)
    return(NULL)
  }
  if (!file.exists(CELL_ANNOTATIONS_PATH)) {
    warning("Cell annotations file not found: ", CELL_ANNOTATIONS_PATH)
    return(NULL)
  }
  message("  Loading cell IDs from annotations file...")
  annotations <- fread(CELL_ANNOTATIONS_PATH, select = "cellID")
  all_cell_ids <- annotations$cellID
  target_indices <- which(all_cell_ids %in% cell_ids)
  message("  Found ", length(target_indices), " / ", length(cell_ids),
          " cells in joint embedding")
  if (length(target_indices) == 0) {
    warning("No cells matched in joint embedding")
    return(NULL)
  }
  message("  Reading joint embedding CSV (2.7GB)...")
  joint <- fread(JOINT_EMBEDDING_PATH)
  joint_subset <- joint[target_indices, ]
  target_cell_ids <- all_cell_ids[target_indices]
  joint_mat <- as.matrix(joint_subset[, -1, with = FALSE])
  rownames(joint_mat) <- target_cell_ids
  message("  Joint embedding subset: ", nrow(joint_mat), " cells x ",
          ncol(joint_mat), " dims")
  return(joint_mat)
}

rename_cells_to_ihbca <- function(obj, mapping, study) {
  original_cells <- colnames(obj)
  n_original <- length(original_cells)
  component_col <- paste0(study, "_cell_id")
  if (!component_col %in% names(mapping)) {
    component_col <- names(mapping)[grep("cell_id", names(mapping))[2]]
  }
  id_lookup <- setNames(mapping$ihbca_cell_id, mapping[[component_col]])
  new_cells <- id_lookup[original_cells]
  n_mapped <- sum(!is.na(new_cells))
  n_unmapped <- sum(is.na(new_cells))
  message("  Cell renaming: ", n_mapped, " mapped, ", n_unmapped, " unmapped")
  if (n_unmapped > 0) {
    message("  Filtering out ", n_unmapped, " unmapped cells")
    keep_idx <- !is.na(new_cells)
    obj <- obj[, keep_idx]
    new_cells <- new_cells[keep_idx]
  }
  obj <- RenameCells(obj, new.names = new_cells)
  return(obj)
}

# Fixed attach_metadata using iHBCA cell inventory lookup
attach_metadata <- function(obj, cell_mapping, donor_metadata, study_prefix) {
  ihbca_cellids <- colnames(obj)
  
  # Load iHBCA cell inventory for cell_id -> patient_id mapping
  if (!file.exists(CELL_INVENTORY_PATH)) {
    stop("Cell inventory not found: ", CELL_INVENTORY_PATH)
  }
  
  message("  Loading cell inventory for donor ID lookup...")
  inventory <- fread(CELL_INVENTORY_PATH, stringsAsFactors = FALSE)
  
  # Create lookup: cell_id -> patient_id
  cell_to_patient <- setNames(inventory$patient_id, inventory$cell_id)
  
  # Look up donor IDs from inventory
  donor_ids <- cell_to_patient[ihbca_cellids]
  
  # Check inventory lookup coverage
  n_found <- sum(!is.na(donor_ids))
  n_missing <- sum(is.na(donor_ids))
  if (n_missing > 0) {
    warning("  ", n_missing, " cells not found in iHBCA inventory")
    missing_ids <- ihbca_cellids[is.na(donor_ids)]
    message("  Sample missing cell IDs: ", paste(head(missing_ids, 3), collapse = ", "))
  }
  
  # Create cell metadata
  cell_meta <- data.frame(
    cell_id = ihbca_cellids,
    ihbca_donor_id = donor_ids,
    stringsAsFactors = FALSE
  )
  
  # Normalize whitespace before join (guards against CSV whitespace artifacts)
  cell_meta$ihbca_donor_id <- trimws(cell_meta$ihbca_donor_id)
  donor_metadata$ihbca_donor_id <- trimws(donor_metadata$ihbca_donor_id)

  # Join with donor metadata
  cell_meta <- left_join(cell_meta, donor_metadata, by = "ihbca_donor_id")
  rownames(cell_meta) <- cell_meta$cell_id
  
  # Check join success
  condition_cols <- c("age_binary", "parity_binary", "risk_status_binary", "menopausal_status_binary")
  check_col <- NULL
  for (col in condition_cols) {
    if (col %in% names(cell_meta) && any(!is.na(cell_meta[[col]]))) {
      check_col <- col
      break
    }
  }
  
  if (is.null(check_col)) {
    n_matched <- sum(!is.na(cell_meta$ihbca_donor_id) &
                     cell_meta$ihbca_donor_id %in% donor_metadata$ihbca_donor_id)
    n_total <- nrow(cell_meta)
    match_rate <- n_matched / n_total
    message("  Metadata join (donor ID match): ", n_matched, "/", n_total, " cells (",
            round(100 * match_rate, 1), "%)")
  } else {
    n_matched <- sum(!is.na(cell_meta[[check_col]]))
    n_total <- nrow(cell_meta)
    match_rate <- n_matched / n_total
    message("  Metadata join (", check_col, "): ", n_matched, "/", n_total, " cells (",
            round(100 * match_rate, 1), "%)")
  }
  
  if (match_rate < 0.9) {
    warning("Low metadata match rate - check donor IDs")
    message("  Sample donor IDs from cells: ",
            paste(head(unique(na.omit(donor_ids)), 5), collapse = ", "))
    message("  Donor IDs in metadata: ",
            paste(head(donor_metadata$ihbca_donor_id, 5), collapse = ", "))
  }
  
  obj <- AddMetaData(obj, metadata = cell_meta)
  obj@misc$metadata_match_rate <- match_rate
  return(obj)
}

build_milo <- function(seurat_obj, embedding_key, embedding_dims) {
  message("  Building Milo object...")
  sce <- as.SingleCellExperiment(seurat_obj)
  milo <- Milo(sce)
  available <- reducedDimNames(milo)
  if (!embedding_key %in% available) {
    matches <- available[tolower(available) == tolower(embedding_key)]
    if (length(matches) > 0) {
      embedding_key <- matches[1]
    } else {
      stop("Embedding not found in SCE: ", embedding_key)
    }
  }
  actual_dims <- min(embedding_dims, ncol(reducedDim(milo, embedding_key)))
  message("  Building k-NN graph (k=", MILO_K, ")...")
  milo <- buildGraph(milo, k = MILO_K, d = actual_dims, reduced.dim = embedding_key)
  message("  Computing neighborhoods (prop=", MILO_PROP, ")...")
  milo <- makeNhoods(milo, prop = MILO_PROP, k = MILO_K, d = actual_dims,
                     refined = TRUE, reduced_dims = embedding_key)
  n_nhoods <- ncol(nhoods(milo))
  message("  Created ", n_nhoods, " neighborhoods")
  message("  Computing neighborhood distances...")
  milo <- calcNhoodDistance(milo, d = actual_dims, reduced.dim = embedding_key)
  return(milo)
}

save_outputs <- function(study, seurat_obj, milo, embedding_native,
                         embedding_joint, cell_metadata, config, output_dir) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
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
  saveRDS(seurat_obj, seurat_file)
  message("  Saved: ", basename(seurat_file))
  
  milo_file <- file.path(output_dir, "milo.rds")
  saveRDS(milo, milo_file)
  message("  Saved: ", basename(milo_file))
  
  summary <- list(
    study = study,
    created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    source_file = config$file,
    n_cells = ncol(seurat_obj),
    n_genes = nrow(seurat_obj),
    n_neighborhoods = ncol(nhoods(milo)),
    metadata_match_rate = round(seurat_obj@misc$metadata_match_rate, 4),
    embedding_native = list(key = config$embedding_key, dims = config$embedding_dims),
    embedding_joint = list(
      available = !is.null(embedding_joint),
      dims = if (!is.null(embedding_joint)) ncol(embedding_joint) else NA
    ),
    milo_params = list(k = MILO_K, prop = MILO_PROP),
    outputs = list(
      metadata = "metadata.csv",
      embedding_native = "embedding_native.csv",
      embedding_joint = if (!is.null(embedding_joint)) "embedding_joint.csv" else NA,
      seurat = "seurat.rds",
      milo = "milo.rds"
    )
  )
  
  summary_file <- file.path(output_dir, "build_summary.yaml")
  write_yaml(summary, summary_file)
  message("  Saved: ", basename(summary_file))
  return(summary)
}

# =============================================================================
# Main Pipeline
# =============================================================================

build_study_objects <- function(study, skip_joint = FALSE, start_step = 1L) {
  if (!study %in% names(STUDY_CONFIG)) {
    stop("Unknown study: ", study, ". Available: ", paste(names(STUDY_CONFIG), collapse = ", "))
  }
  
  config <- STUDY_CONFIG[[study]]
  output_dir <- file.path(DATA_PREP_PATH, "outputs", "study_objects", study)
  
  message("\n", strrep("=", 60))
  message("Building study objects for: ", toupper(study))
  message("  Source: ", basename(config$file))
  message("  Output: ", output_dir)
  if (start_step > 1) message("  Checkpoint: resuming from step ", start_step)
  message(strrep("=", 60))

  if (start_step >= 7) {
    # --- Checkpoint resume: load existing seurat.rds, skip steps 1-6 ---
    seurat_file <- file.path(output_dir, "seurat.rds")
    if (!file.exists(seurat_file)) {
      stop("No seurat.rds found for checkpoint resume at: ", seurat_file)
    }
    message("\n[CHECKPOINT] Loading existing seurat.rds...")
    obj <- readRDS(seurat_file)
    message("  Loaded: ", ncol(obj), " cells x ", nrow(obj), " genes")

    # Reload embeddings from disk for save_outputs
    message("[CHECKPOINT] Loading saved embeddings...")
    native_file <- file.path(output_dir, "embedding_native.csv")
    if (!file.exists(native_file)) {
      stop("No embedding_native.csv found for checkpoint resume at: ", native_file)
    }
    embedding_native <- read.csv(native_file, row.names = 1)

    joint_file <- file.path(output_dir, "embedding_joint.csv")
    embedding_joint <- if (file.exists(joint_file)) {
      read.csv(joint_file, row.names = 1)
    } else {
      NULL
    }

    cell_metadata <- obj@meta.data
    cell_metadata$cell_id <- rownames(cell_metadata)

    message("  Native embedding: ", nrow(embedding_native), " x ", ncol(embedding_native))
    if (!is.null(embedding_joint)) {
      message("  Joint embedding: ", nrow(embedding_joint), " x ", ncol(embedding_joint))
    }
  } else {
    # --- Full pipeline: steps 1-6 ---
    message("\n[1/7] Loading cell mapping...")
    cell_mapping <- load_cell_mapping(study)

    message("\n[2/7] Loading donor metadata...")
    donor_metadata <- load_donor_metadata()
    # Map study name to harmonized metadata study column value
    # Pal substudies (pal_norm_*) all use "pal" in harmonized metadata
    study_key <- tolower(sub("_norm.*", "", study))
    donor_metadata <- donor_metadata[tolower(donor_metadata$study) == study_key, ]
    message("  Filtered to ", nrow(donor_metadata), " donors for ", study, " (study_key: ", study_key, ")")

    message("\n[3/7] Loading component RDS...")
    obj <- load_rds_object(config$file)

    message("\n[4/7] Renaming cells to iHBCA format...")
    obj <- rename_cells_to_ihbca(obj, cell_mapping, study)

    message("\n[5/7] Attaching harmonized metadata...")
    obj <- attach_metadata(obj, cell_mapping, donor_metadata, config$study_prefix)

    message("\n[6/7] Extracting embeddings...")
    embedding_native <- extract_embedding(obj, config$embedding_key, config$embedding_dims)

    embedding_joint <- NULL
    if (!skip_joint) {
      embedding_joint <- load_joint_embedding(colnames(obj))
      if (!is.null(embedding_joint)) {
        common_cells <- intersect(colnames(obj), rownames(embedding_joint))
        if (length(common_cells) < ncol(obj)) {
          message("  Warning: ", ncol(obj) - length(common_cells), " cells missing from joint embedding")
        }
        embedding_joint <- embedding_joint[common_cells, , drop = FALSE]
      }
    } else {
      message("  Skipping joint embedding (--skip-joint)")
    }

    cell_metadata <- obj@meta.data
    cell_metadata$cell_id <- rownames(cell_metadata)
  }

  message("\n[7/7] Building Milo object...")
  milo <- build_milo(obj, config$embedding_key, config$embedding_dims)
  
  message("\n[SAVE] Saving outputs...")
  summary <- save_outputs(study, obj, milo, embedding_native, embedding_joint,
                          cell_metadata, config, output_dir)
  
  message("\n", strrep("=", 60))
  message("SUCCESS: ", study, " objects built")
  message("  Cells: ", summary$n_cells)
  message("  Neighborhoods: ", summary$n_neighborhoods)
  message("  Metadata match: ", round(100 * summary$metadata_match_rate, 1), "%")
  message("  Joint embedding: ", if (summary$embedding_joint$available) "YES" else "NO")
  message("  Output: ", output_dir)
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
  basename(sub("--file=", "", script_arg)) == "build_study_objects.R"
}

if (!interactive() && .is_main_script()) {
  parser <- ArgumentParser(description = "Build study objects for Phase B analysis")
  parser$add_argument("--study", type = "character", required = TRUE,
                      help = "Study name (gray, kumar, murrow, nee, twigger, pal_norm_*, reed)")
  parser$add_argument("--skip-joint", action = "store_true",
                      help = "Skip joint embedding extraction (faster for testing)")
  parser$add_argument("--start-step", type = "integer", default = 1L,
                      help = "Resume from step N (1=full, 7=milo-only)")
  parser$add_argument("--all", action = "store_true", help = "Build all studies")
  
  args <- parser$parse_args()
  
  if (args$all) {
    for (study in names(STUDY_CONFIG)) {
      tryCatch({
        build_study_objects(study, skip_joint = args$skip_joint, start_step = args$start_step)
      }, error = function(e) {
        message("\nERROR building ", study, ": ", e$message, "\n")
      })
    }
  } else {
    build_study_objects(args$study, skip_joint = args$skip_joint, start_step = args$start_step)
  }
}
