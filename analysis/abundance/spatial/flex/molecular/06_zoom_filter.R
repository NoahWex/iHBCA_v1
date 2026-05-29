#!/usr/bin/env Rscript
# 06_zoom_filter.R — extract fibroblast + basal NhoodGroup zoom from canonical
# L2.0s markers.csv. Splits canonical markers into compartment-specific files
# and copies within-parent pairwise outputs into the zoom dir.
#
# Inputs:
#   --canonical-dir   outputs/markers_canonical/{contrast}/{scope}/{tier}/
#                     (must contain markers.csv, group_summary.csv,
#                      pairwise_fibro_major_split.csv, pairwise_bmyo_myo_split.csv)
#   --out-dir         outputs/markers_zoom/{contrast}/{tier}/
#
# Scope definition (L2.0s parent_label):
#   fibroblast scope: Fibro-major, Fibro-myo, Fibro-prematrix
#   basal scope:      BMYO-myo, BMYO-basal

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
})

parser <- ArgumentParser()
parser$add_argument("--canonical-dir", required = TRUE)
parser$add_argument("--out-dir",       required = TRUE)
args <- parser$parse_args()

dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)

FIBRO_PARENTS <- c("Fibro-major", "Fibro-myo", "Fibro-prematrix")
BASAL_PARENTS <- c("BMYO-myo", "BMYO-basal")

gs <- fread(file.path(args$canonical_dir, "group_summary.csv"))
fibro_groups <- gs[parent_label %in% FIBRO_PARENTS & passed_outlier_filter == TRUE,
                   NhoodGroup]
basal_groups <- gs[parent_label %in% BASAL_PARENTS & passed_outlier_filter == TRUE,
                   NhoodGroup]
cat(sprintf("[zoom] fibroblast NhoodGroups (%d): %s\n",
            length(fibro_groups), paste(fibro_groups, collapse = ", ")))
cat(sprintf("[zoom] basal NhoodGroups (%d): %s\n",
            length(basal_groups), paste(basal_groups, collapse = ", ")))

cat("[zoom] reading canonical markers.csv...\n")
markers <- fread(file.path(args$canonical_dir, "markers.csv"))
cat(sprintf("[zoom] canonical markers: %d rows across %d groups\n",
            nrow(markers), uniqueN(markers$NhoodGroup)))

fibro_markers <- markers[NhoodGroup %in% fibro_groups]
basal_markers <- markers[NhoodGroup %in% basal_groups]

fwrite(fibro_markers, file.path(args$out_dir, "fibroblast_markers.csv"))
fwrite(basal_markers, file.path(args$out_dir, "basal_markers.csv"))
cat(sprintf("[zoom] wrote fibroblast_markers.csv (%d rows)\n", nrow(fibro_markers)))
cat(sprintf("[zoom] wrote basal_markers.csv (%d rows)\n", nrow(basal_markers)))

# Copy pairwise files (within-parent splits already computed)
pairwise_files <- c(
  fibroblast = "pairwise_fibro_major_split.csv",
  basal      = "pairwise_bmyo_myo_split.csv"
)
for (tag in names(pairwise_files)) {
  src <- file.path(args$canonical_dir, pairwise_files[[tag]])
  dst <- file.path(args$out_dir, sprintf("%s_pairwise.csv", tag))
  if (file.exists(src)) {
    file.copy(src, dst, overwrite = TRUE)
    cat(sprintf("[zoom] copied %s -> %s\n", basename(src), basename(dst)))
  } else {
    cat(sprintf("[zoom] WARN: pairwise source missing: %s\n", src))
  }
}

# Manifest
manifest <- list(
  fibro_groups = fibro_groups,
  basal_groups = basal_groups,
  fibro_n_rows = nrow(fibro_markers),
  basal_n_rows = nrow(basal_markers),
  source = args$canonical_dir,
  generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
)
writeLines(
  c(sprintf("source: %s", manifest$source),
    sprintf("generated: %s", manifest$generated),
    sprintf("fibroblast_groups: [%s]", paste(manifest$fibro_groups, collapse = ", ")),
    sprintf("basal_groups: [%s]", paste(manifest$basal_groups, collapse = ", ")),
    sprintf("fibroblast_marker_rows: %d", manifest$fibro_n_rows),
    sprintf("basal_marker_rows: %d", manifest$basal_n_rows)),
  file.path(args$out_dir, "zoom_manifest.yaml")
)
cat(sprintf("\n[zoom] outputs written to %s\n", args$out_dir))
