#!/usr/bin/env Rscript
# plot_da_feature.R
# Per-cell DA feature plots over UMAP_scVI.
# Pattern from clinical_cohort_summary_20260426/scripts/22_aging_umap_preview.R:
# project per-nhood logFC to cells via nhood membership matrix, restricted to
# significant nhoods, then plot UMAP coloured by per-cell mean logFC.
#
# Output (per contrast):
#   reports/figures/da_umap_<contrast>.pdf

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(SummarizedExperiment); library(miloR)
  library(Matrix); library(dplyr); library(readr); library(yaml); library(ggplot2)
  library(scales)
})
have_repel <- requireNamespace("ggrepel", quietly = TRUE)

`%||%` <- function(a, b) if (is.null(a)) b else a

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "plot_da_feature.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths
inquiry <- cfg$inquiry
fig_dir <- file.path(paths$outputs$reports, "figures")
ensure_dir(fig_dir)

cat(sprintf("=== DA per-cell UMAP plots for '%s' ===\n", inquiry$inquiry))

# Load milo + UMAP + nhood membership ONCE
cat("Loading milo + patches...\n")
source(paths$inputs$patch_milor)
milo <- readRDS(paths$inputs$milo)
nh_mat <- nhoods(milo)
cat(sprintf("milo cells=%d nhoods=%d\n", ncol(milo), ncol(nh_mat)))

umap_key <- intersect(c("UMAP_scVI", "UMAP", "X_umap", "umap"),
                      reducedDimNames(milo))[1]
if (is.na(umap_key)) stop("No UMAP reducedDim. Available: ",
                          paste(reducedDimNames(milo), collapse = ", "))
umap_df <- as.data.frame(reducedDim(milo, umap_key))
colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")
umap_df$cell_id <- colnames(milo)
cat(sprintf("UMAP key: %s\n", umap_key))

# Load L1-grain labels for centroid overlay (publication/.../labels_full.csv).
# L1 = `lineage` column (~20 broad cell types). Centroids = median UMAP per L1.
labels_path <- paths$inputs$labels
l1_centroids <- NULL
if (file.exists(labels_path)) {
  cat(sprintf("Loading L1 labels from %s\n", labels_path))
  lab <- read_csv(labels_path, show_col_types = FALSE,
                  col_select = c("cell_id", "compartment", "lineage", "is_artifact"))
  lab <- lab %>% filter(!isTRUE(as.logical(is_artifact)),
                         !is.na(lineage), lineage != "")
  umap_lab <- umap_df %>% inner_join(lab, by = "cell_id")
  cat(sprintf("UMAP cells joined to labels: %d / %d\n",
              nrow(umap_lab), nrow(umap_df)))
  l1_centroids <- umap_lab %>%
    group_by(compartment, lineage) %>%
    summarise(UMAP_1 = median(UMAP_1, na.rm = TRUE),
              UMAP_2 = median(UMAP_2, na.rm = TRUE),
              n_cells = dplyr::n(),
              .groups = "drop") %>%
    filter(n_cells >= 200)
  cat(sprintf("L1 centroids: %d (>= 200 cells each)\n", nrow(l1_centroids)))
} else {
  cat("WARN: labels file not found — skipping L1 centroid overlay\n")
}

# Per-contrast plot
contrast_dirs <- list.dirs(paths$outputs$stageD, recursive = FALSE)
for (cdir in contrast_dirs) {
  cname <- basename(cdir)
  da_csv <- file.path(cdir, "da_results.csv")
  if (!file.exists(da_csv)) next
  da <- read_csv(da_csv, show_col_types = FALSE)
  rs <- yaml::read_yaml(file.path(cdir, "run_summary.yaml"))
  fdr_t <- rs$spatial_fdr %||% 0.05

  da <- da[order(da$Nhood), ]
  if (!all(da$Nhood == seq_len(nrow(da)))) {
    cat(sprintf("WARN %s: Nhood column not 1..N — projection may be off\n", cname))
    next
  }

  sig_mask <- !is.na(da$SpatialFDR) & da$SpatialFDR < fdr_t
  sig_idx <- which(sig_mask)
  cat(sprintf("\n--- %s: %d sig nhoods (FDR<%.2f) ---\n",
              cname, length(sig_idx), fdr_t))
  if (length(sig_idx) == 0) {
    cat("  no sig nhoods — skipping\n"); next
  }

  nh_sig <- nh_mat[, sig_idx, drop = FALSE]
  logfc  <- da$logFC[sig_idx]
  num <- as.numeric(nh_sig %*% logfc)
  den <- as.numeric(Matrix::rowSums(nh_sig))
  cell_logfc <- ifelse(den > 0, num / den, NA_real_)

  d <- umap_df
  d$logFC <- cell_logfc
  d_have <- d[!is.na(d$logFC), ]
  d_bg   <- d[is.na(d$logFC), ]
  cat(sprintf("  cells with sig membership: %d / %d (%.1f%%)\n",
              nrow(d_have), nrow(d), 100 * nrow(d_have) / nrow(d)))

  lim <- max(abs(quantile(d_have$logFC, c(0.01, 0.99), na.rm = TRUE)))
  d_have$logFC_capped <- pmin(pmax(d_have$logFC, -lim), lim)

  if (nrow(d_have) > 250000) d_have <- d_have[sample.int(nrow(d_have), 250000), ]
  if (nrow(d_bg)   > 200000) d_bg   <- d_bg[sample.int(nrow(d_bg), 200000), ]

  p <- ggplot() +
    geom_point(data = d_bg,
               aes(x = UMAP_1, y = UMAP_2),
               color = "gray85", size = 0.12, alpha = 0.4) +
    geom_point(data = d_have,
               aes(x = UMAP_1, y = UMAP_2, color = logFC_capped),
               size = 0.18, alpha = 0.7) +
    scale_color_gradient2(low = "#2166AC", mid = "gray92", high = "#B2182B",
                          midpoint = 0,
                          limits = c(-lim, lim), oob = scales::squish,
                          name = "logFC")
  if (!is.null(l1_centroids)) {
    if (have_repel) {
      p <- p + ggrepel::geom_text_repel(
        data = l1_centroids,
        aes(x = UMAP_1, y = UMAP_2, label = lineage),
        size = 2.4, fontface = "bold", color = "black",
        bg.color = "white", bg.r = 0.12,
        max.overlaps = Inf, segment.size = 0.15,
        segment.color = "grey50", min.segment.length = 0,
        box.padding = 0.25, point.padding = 0.15
      )
    } else {
      p <- p + geom_text(
        data = l1_centroids,
        aes(x = UMAP_1, y = UMAP_2, label = lineage),
        size = 2.4, fontface = "bold", color = "black"
      )
    }
  }
  p <- p +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme_void(base_size = 8) +
    theme(legend.position = "right")

  out <- file.path(fig_dir, sprintf("da_umap_%s.pdf", cname))
  ggsave(out, p, width = 9, height = 7, dpi = 150, device = "pdf")
  cat(sprintf("  wrote: %s\n", out))
}

cat("\n=== DA UMAP feature plotting done ===\n")
