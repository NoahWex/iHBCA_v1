#!/usr/bin/env Rscript
# plot_nhoodgroup_umap.R
# Per-cell NhoodGroup projection onto UMAP, one panel per contrast.
#
# Method (mirror of plot_da_feature.R logic for categorical labels):
#   1. Load milo + UMAP + nhood membership matrix.
#   2. For each contrast with stageF1 nhood_groups.csv, build per-nhood
#      NhoodGroup vector.
#   3. For each cell, take the FIRST nhood it belongs to whose NhoodGroup is
#      assigned (i.e., the dominant signal-bearing nhood). Cells in no
#      assigned nhood -> NA / gray.
#   4. Plot UMAP with high-contrast palette, no legend.
#
# Outputs:
#   reports/figures/da_umap_nhoodgroup_<contrast>.pdf  per contrast

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(SummarizedExperiment); library(miloR)
  library(Matrix); library(dplyr); library(readr); library(yaml); library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "plot_nhoodgroup_umap.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths
fig_dir <- file.path(paths$outputs$reports, "figures")
ensure_dir(fig_dir)

cat("Loading milo + patches...\n")
source(paths$inputs$patch_milor)
milo <- readRDS(paths$inputs$milo)
nh_mat <- nhoods(milo)
cat(sprintf("milo cells=%d nhoods=%d\n", ncol(milo), ncol(nh_mat)))

umap_key <- intersect(c("UMAP_scVI","UMAP","X_umap","umap"), reducedDimNames(milo))[1]
if (is.na(umap_key)) stop("No UMAP reducedDim. Found: ",
                          paste(reducedDimNames(milo), collapse=","))
umap_df <- as.data.frame(reducedDim(milo, umap_key))
colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")

# Stable contrasting palette: hash group ID -> hue, distribute across the
# wheel so adjacent group numbers don't collapse to similar colors.
make_palette <- function(group_ids) {
  ids <- sort(unique(group_ids[!is.na(group_ids)]))
  n <- length(ids)
  # Use Polychrome if available, else HCL with stride to avoid neighbor clash
  if (requireNamespace("Polychrome", quietly = TRUE)) {
    pal <- Polychrome::createPalette(n, c("#5A8BCB","#E07B39","#9B2226"))
    names(pal) <- as.character(ids)
  } else {
    # Manual hue stride: shuffle to spread similar IDs
    set.seed(7)
    hues <- (seq(0, 1, length.out = n + 1)[-(n+1)] + runif(n, 0, 0.05)) %% 1
    pal <- hcl(h = hues * 360, c = 75, l = 60)
    names(pal) <- as.character(ids[sample.int(n)])
    pal <- pal[as.character(ids)]
  }
  pal
}

# Iterate every Stage F.1 contrast directory
f1_dir <- paths$outputs$stageF1
contrast_dirs <- list.dirs(f1_dir, recursive = FALSE)
if (length(contrast_dirs) == 0) {
  cat("No F.1 outputs at", f1_dir, "\n"); quit(status = 0)
}

VIABLE_MIN_NHOODS <- 10

# Render two variants per contrast: all groups (existing) + limma-viable subset
# (groups with n_nhoods >= VIABLE_MIN_NHOODS, the F.3 pseudobulk viability proxy).
# Variants are emitted by looping over a small spec list; rest of the per-contrast
# logic is shared.
variants <- list(
  list(suffix = "",              filter_fn = function(ng) ng,
       label  = "all groups"),
  list(suffix = "_limma_viable", filter_fn = function(ng) ng %>%
                                              filter(n_nhoods_in_group >= VIABLE_MIN_NHOODS),
       label  = sprintf("groups with n_nhoods >= %d", VIABLE_MIN_NHOODS))
)

for (cdir in contrast_dirs) {
  cname <- basename(cdir)
  ng_csv <- file.path(cdir, "nhood_groups.csv")
  if (!file.exists(ng_csv) || file.size(ng_csv) < 200) {
    cat(sprintf("  skip %s: no nhood_groups.csv\n", cname)); next
  }
  ng_full <- read_csv(ng_csv, show_col_types = FALSE)
 for (v in variants) {
  ng <- v$filter_fn(ng_full) %>%
    select(Nhood, NhoodGroup) %>%
    arrange(Nhood)

  # Build per-nhood group vector aligned to nh_mat columns (Nhood = 1..N).
  nh_groups <- rep(NA_integer_, ncol(nh_mat))
  nh_groups[ng$Nhood] <- as.integer(ng$NhoodGroup)
  n_assigned_nh <- sum(!is.na(nh_groups))
  n_groups <- n_distinct(nh_groups, na.rm = TRUE)
  cat(sprintf("\n--- %s: %d assigned nhoods in %d groups ---\n",
              cname, n_assigned_nh, n_groups))
  if (n_assigned_nh == 0) {
    cat("  no assigned nhoods, skipping\n"); next
  }

  # Per-cell: first assigned-group nhood it belongs to.
  # nh_mat is cells x nhoods. For each cell row, find the index of the first
  # nhood column where the cell is a member AND the nhood has a group.
  assigned_cols <- which(!is.na(nh_groups))
  nh_assigned <- nh_mat[, assigned_cols, drop = FALSE]
  group_per_assigned_col <- nh_groups[assigned_cols]

  # For each cell: first column index in nh_assigned where membership is TRUE
  # Convert to triplets, take min(j) per i.
  trip <- summary(nh_assigned)  # data.frame i, j, x
  trip <- trip[trip$x > 0, c("i","j")]
  # Keep first j per i (smallest column → arbitrary but deterministic)
  trip <- trip[order(trip$i, trip$j), ]
  trip <- trip[!duplicated(trip$i), ]
  cell_group <- rep(NA_integer_, nrow(nh_assigned))
  cell_group[trip$i] <- group_per_assigned_col[trip$j]
  cat(sprintf("  cells with group: %d / %d (%.1f%%)\n",
              sum(!is.na(cell_group)), length(cell_group),
              100 * sum(!is.na(cell_group)) / length(cell_group)))

  d <- umap_df
  d$NhoodGroup <- factor(cell_group)
  d_have <- d[!is.na(d$NhoodGroup), ]
  d_bg   <- d[is.na(d$NhoodGroup), ]
  if (nrow(d_have) > 250000) d_have <- d_have[sample.int(nrow(d_have), 250000), ]
  if (nrow(d_bg)   > 200000) d_bg   <- d_bg[sample.int(nrow(d_bg), 200000), ]

  pal <- make_palette(d_have$NhoodGroup)

  p <- ggplot() +
    geom_point(data = d_bg, aes(x = UMAP_1, y = UMAP_2),
               color = "gray88", size = 0.10, alpha = 0.4) +
    geom_point(data = d_have,
               aes(x = UMAP_1, y = UMAP_2, color = NhoodGroup),
               size = 0.18, alpha = 0.85) +
    scale_color_manual(values = pal, na.value = "gray88", guide = "none") +
    labs(title = sprintf("%s — NhoodGroups projected onto UMAP (%s)",
                          cname, v$label),
         subtitle = sprintf("%d groups across %d nhoods. Gray = cell in no assigned group.",
                            n_groups, n_assigned_nh),
         x = "UMAP 1", y = "UMAP 2") +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(color = "gray40", size = 9))

  out <- file.path(fig_dir,
                    sprintf("da_umap_nhoodgroup_%s%s.pdf", cname, v$suffix))
  ggsave(out, p, width = 9, height = 7, dpi = 150, device = "pdf")
  cat(sprintf("  wrote: %s\n", out))
 }  # end variants loop
}

cat("\n=== nhoodgroup UMAP done ===\n")
