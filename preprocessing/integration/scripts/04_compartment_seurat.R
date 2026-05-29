#!/usr/bin/env Rscript
# Per-compartment Seurat integration — comparative scib benchmark only.
#
# Three methods supported via --method:
#   sketch   Seurat v5 sketch-based RPCA integration on a balanced ncells-per-
#            dataset sketch, then ProjectIntegration to recover the full data.
#   rpca     Reciprocal PCA on the full data.
#   harmony  Built-in Seurat v5 HarmonyIntegration on the full data.
#
# Inputs: a 10x HDF5 v3 counts file (genes x cells) + per-compartment metadata
# CSV. The HDF5 path replaces an earlier MTX route — Matrix::readMM hits R's
# int32 scan() length limit on the epithelial compartment (nnz > 2^31), while
# 10x HDF5 v3 stores indices as int64 and is read natively by
# Seurat::Read10X_h5.
#
# Included as a comparative integration baseline in the scib benchmark suite.
# Not part of the canonical pipeline — scANVI at n_latent=50 is the published
# winner.
#
# Outputs: RDS + CSV exports (embedding, UMAP, Leiden across 10 resolutions).
#
# Usage:
#   Rscript 04_compartment_seurat.R \
#     --counts-h5  <{comp}_counts.h5> \
#     --metadata   <{comp}_metadata.csv> \
#     --compartment Immune \
#     --method sketch \
#     --output-dir <out_dir>

suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
})

# Seurat v5 IntegrateLayers parallelizes across batches via the future
# framework. The default 500 MB serialization cap is well below the size of a
# 1M-cell sparse matrix passed as a global to workers, so raise it to 100 GB.
options(future.globals.maxSize = 100 * 1024^3)

# %||% null-coalesce helper (defined early so it's available in args block below)
`%||%` <- function(x, y) if (is.null(x)) y else x

# Base R argument parser (avoids argparse dependency).
# Supports `--key value` and `--key=value` forms.
parse_args <- function(argv) {
    args <- list()
    i <- 1
    while (i <= length(argv)) {
        a <- argv[i]
        if (grepl("^--", a)) {
            key <- sub("^--", "", a)
            if (grepl("=", key)) {
                kv <- strsplit(key, "=", fixed = TRUE)[[1]]
                args[[kv[1]]] <- kv[2]
                i <- i + 1
            } else if (i + 1 <= length(argv) && !grepl("^--", argv[i + 1])) {
                args[[key]] <- argv[i + 1]
                i <- i + 2
            } else {
                args[[key]] <- TRUE  # boolean flag
                i <- i + 1
            }
        } else {
            i <- i + 1
        }
    }
    args
}

raw_args <- parse_args(commandArgs(trailingOnly = TRUE))

# Validate required + apply defaults
required <- c("counts-h5", "metadata", "compartment", "method", "output-dir")
for (r in required) {
    if (is.null(raw_args[[r]])) stop(sprintf("Missing required arg: --%s", r))
}

valid_compartments <- c("Immune", "Epithelial", "Stromal")
if (!raw_args$compartment %in% valid_compartments) {
    stop(sprintf("--compartment must be one of: %s",
                 paste(valid_compartments, collapse = ", ")))
}
if (!raw_args$method %in% c("sketch", "rpca", "harmony")) {
    stop("--method must be 'sketch', 'rpca', or 'harmony'")
}

# Build args list with R-friendly underscore keys + integer defaults
args <- list(
    counts_h5    = raw_args[["counts-h5"]],
    metadata     = raw_args$metadata,
    compartment  = raw_args$compartment,
    method       = raw_args$method,
    output_dir   = raw_args[["output-dir"]],
    n_hvg        = as.integer(raw_args[["n-hvg"]] %||% 4000L),
    sketch_cells = as.integer(raw_args[["sketch-cells"]] %||% 5000L),
    n_dims       = as.integer(raw_args[["n-dims"]] %||% 50L),
    seed         = as.integer(raw_args$seed %||% 42L)
)

set.seed(args$seed)
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

cat("=================================================================\n")
cat(sprintf("Per-Compartment Seurat (%s): %s\n", args$method, args$compartment))
cat("=================================================================\n")

# --- Load 10x HDF5 via Seurat::Read10X_h5 ---
# The HDF5 file is produced by npz_to_10x_h5.py with:
#   matrix/features/id    = Ensembl IDs     (rownames when use.names=FALSE)
#   matrix/features/name  = gene symbols
#   matrix/barcodes       = cell IDs        (colnames)
# use.names=FALSE returns a dgCMatrix with Ensembl IDs as rownames.
cat(sprintf("Loading 10x HDF5: %s\n", args$counts_h5))
t0 <- Sys.time()
mat <- Seurat::Read10X_h5(args$counts_h5, use.names = FALSE)
cat(sprintf("  Loaded in %.1f s\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
cat(sprintf("  Matrix class: %s, dim: %d x %d, nnz: %d\n",
            class(mat)[[1]], nrow(mat), ncol(mat), length(mat@x)))
# Read10X_h5 returns dgCMatrix directly with rownames/colnames set — no
# coercion or manual label assignment needed.

# --- Load metadata ---
# Per-compartment MTX is already compartment-specific (one MTX per compartment).
# The metadata file matches it row-for-row by cell_id.
cat("Loading metadata...\n")
meta <- read.csv(args$metadata, stringsAsFactors = FALSE)
required_cols <- c("dataset", "level0_annotation")
stopifnot(all(required_cols %in% colnames(meta)))
stopifnot(nrow(meta) == ncol(mat))
# Confirm the metadata is for the right compartment
meta_compartments <- unique(meta$level0_annotation)
if (length(meta_compartments) != 1 || meta_compartments[1] != args$compartment) {
    stop(sprintf("Metadata compartment mismatch: expected '%s', got '%s'",
                 args$compartment, paste(meta_compartments, collapse = ",")))
}
cat(sprintf("  Metadata: %d rows, compartment: %s\n",
            nrow(meta), args$compartment))

# --- Build Seurat object ---
cat("Building Seurat object...\n")
seu <- CreateSeuratObject(counts = mat, project = "iHBCA")
seu$dataset <- meta$dataset
rm(mat); gc()
cat(sprintf("  Seurat: %d cells x %d features\n", ncol(seu), nrow(seu)))

# --- Normalize + HVG ---
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, nfeatures = args$n_hvg, verbose = FALSE)

# --- Method-specific integration ---
if (args$method == "sketch") {
    cat("\n--- Seurat v5 sketch-based integration ---\n")
    # Per the Seurat v5 sketch integration vignette: split layers BEFORE
    # SketchData so the sketch samples ncells per dataset (not globally),
    # guaranteeing balanced representation across studies and avoiding
    # k.weight failures on small-study compartments.
    seu[["RNA"]] <- split(seu[["RNA"]], f = seu$dataset)
    seu <- FindVariableFeatures(seu, nfeatures = args$n_hvg, verbose = FALSE)

    seu <- SketchData(
        object = seu,
        ncells = args$sketch_cells,
        method = "LeverageScore",
        sketched.assay = "sketch",
        var.name = "leverage.score",
        seed = args$seed
    )
    DefaultAssay(seu) <- "sketch"

    # The sketch assay inherits the layer split from RNA — no manual split needed
    seu <- FindVariableFeatures(seu, nfeatures = args$n_hvg, verbose = FALSE)
    seu <- ScaleData(seu, verbose = FALSE)
    seu <- RunPCA(seu, npcs = args$n_dims, verbose = FALSE)

    # Integrate sketch
    seu <- IntegrateLayers(
        object = seu,
        method = RPCAIntegration,
        orig.reduction = "pca",
        new.reduction = "integrated.rpca",
        verbose = FALSE
    )

    # Project full data back from sketch integration
    seu <- ProjectIntegration(
        object = seu,
        sketched.assay = "sketch",
        assay = "RNA",
        reduction = "integrated.rpca"
    )
    integrated_reduction <- "integrated.rpca.full"

    # Re-join RNA layers for downstream consistency. FindNeighbors and
    # FindClusters can behave inconsistently with split layers even when
    # reading from a reduction; rejoining is the safe invariant.
    seu <- JoinLayers(seu)

} else if (args$method == "rpca") {
    cat("\n--- Seurat v5 reciprocal PCA integration ---\n")
    seu <- ScaleData(seu, verbose = FALSE)
    seu <- RunPCA(seu, npcs = args$n_dims, verbose = FALSE)

    # Split layers by dataset for v5 integration
    seu[["RNA"]] <- split(seu[["RNA"]], f = seu$dataset)

    seu <- IntegrateLayers(
        object = seu,
        method = RPCAIntegration,
        orig.reduction = "pca",
        new.reduction = "integrated.rpca",
        verbose = FALSE
    )

    # Re-join layers
    seu <- JoinLayers(seu)
    integrated_reduction <- "integrated.rpca"

} else if (args$method == "harmony") {
    cat("\n--- Seurat v5 Harmony integration ---\n")
    # Seurat v5 has Harmony as a built-in IntegrateLayers method
    # (Seurat::HarmonyIntegration). No external harmony package needed.
    seu <- ScaleData(seu, verbose = FALSE)
    seu <- RunPCA(seu, npcs = args$n_dims, verbose = FALSE)

    seu[["RNA"]] <- split(seu[["RNA"]], f = seu$dataset)

    seu <- IntegrateLayers(
        object = seu,
        method = HarmonyIntegration,
        orig.reduction = "pca",
        new.reduction = "integrated.harmony",
        verbose = FALSE
    )

    seu <- JoinLayers(seu)
    integrated_reduction <- "integrated.harmony"
}

cat(sprintf("  Integration done. Reduction: %s\n", integrated_reduction))

# Force DefaultAssay to RNA before downstream neighbors/clustering. The sketch
# branch leaves DefaultAssay on "sketch" after integration; the rpca/harmony
# branches stay on RNA. Resetting unconditionally guarantees downstream code
# runs on RNA so FindNeighbors and FindClusters reference the same graph.
DefaultAssay(seu) <- "RNA"
cat(sprintf("  DefaultAssay forced to: %s\n", DefaultAssay(seu)))

# --- UMAP + Leiden ---
cat("\nUMAP + clustering...\n")
seu <- RunUMAP(seu, reduction = integrated_reduction,
               dims = 1:args$n_dims, verbose = FALSE)

# Sanity check: confirm the saved UMAP embedding lives on the expected assay
# (RunUMAP infers assay from the reduction).
cat(sprintf("  UMAP cells: %d, DefaultAssay still: %s\n",
            nrow(Embeddings(seu, "umap")), DefaultAssay(seu)))

seu <- FindNeighbors(seu, reduction = integrated_reduction,
                     dims = 1:args$n_dims, verbose = FALSE)

resolutions <- c(0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0)
for (res in resolutions) {
    seu <- FindClusters(seu, resolution = res, verbose = FALSE,
                        cluster.name = paste0("leiden_", res))
    n_clust <- length(unique(seu[[paste0("leiden_", res)]][, 1]))
    cat(sprintf("    res %.1f: %d clusters\n", res, n_clust))
}

# --- Save outputs ---
comp_short <- c("Immune" = "imm", "Epithelial" = "epi", "Stromal" = "stro")[
    args$compartment
]
out_prefix <- file.path(args$output_dir, paste0(comp_short, "_seurat_", args$method))

# Save the Seurat object as RDS (R-native format; converting to h5ad in R is fragile)
saveRDS(seu, paste0(out_prefix, ".rds"))
cat(sprintf("  RDS saved: %s.rds\n", out_prefix))

# Export integrated embedding + UMAP + leiden as CSVs
embedding_df <- data.frame(
    cell_id = colnames(seu),
    Embeddings(seu, reduction = integrated_reduction)
)
write.csv(embedding_df, paste0(out_prefix, "_embedding.csv"), row.names = FALSE)

umap_df <- data.frame(
    cell_id = colnames(seu),
    UMAP_1 = Embeddings(seu, "umap")[, 1],
    UMAP_2 = Embeddings(seu, "umap")[, 2]
)
write.csv(umap_df, paste0(out_prefix, "_umap.csv"), row.names = FALSE)

leiden_cols <- paste0("leiden_", resolutions)
leiden_df <- data.frame(cell_id = colnames(seu),
                        seu[[leiden_cols]])
write.csv(leiden_df, paste0(out_prefix, "_leiden.csv"), row.names = FALSE)

cat("\nDone.\n")
