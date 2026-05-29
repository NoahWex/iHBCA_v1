#!/usr/bin/env Rscript
# Prepare data for the LBridge differential-abundance beeswarm (Fig 4a / Fig 5).
#
# Loads Milo DA results for both the joint Chromium+Xenium cohort and the
# FLEX-only cohort, assigns per-nhood cell types via the LBridge cross-platform
# vocabulary, and serialises the processed data frame for fast iterative
# rendering.  Run once per contrast; the render script reads the cached RDS.
#
# The LBridge vocabulary maps Xenium L1.5 labels (joint rows) and iHBCA L2.0
# / FLEX Track-C labels (FLEX L2S rows) to shared family groups.  Families are
# ordered by joint-cohort median signed logFC within each LBridge compartment.
# FLEX L2S rows whose display label matches the joint row in the same family are
# suppressed to avoid duplication.
#
# Contrasts
#   c1  p1/p2 (NAC-proximal) vs. p3 (peripheral)  →  Fig 4a
#   c3  UOQ vs. non-UOQ                            →  Fig 5
#
# Outputs: <out-dir>/plot_data_<contrast>.rds

suppressPackageStartupMessages({
  library(argparse)
  library(yaml)
  library(data.table)
  library(dplyr)
})

parser <- ArgumentParser()
parser$add_argument("--project-root",   required = TRUE)
parser$add_argument("--contrast",       required = TRUE,
                    choices = c("c1", "c2a", "c2b", "c3", "c3a", "c3b", "c3c"))
parser$add_argument("--purity-min",     type = "double",  default = 0.3)
parser$add_argument("--min-l2s-nhoods", type = "integer", default = 10)
parser$add_argument("--out-dir",        default = NULL)
args <- parser$parse_args()

root    <- normalizePath(args$project_root)
pub     <- file.path(root, "publication")
out_dir <- if (!is.null(args$out_dir)) args$out_dir else
           file.path(root, "drafting_space", "fig4", "outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

paths_yaml <- yaml::read_yaml(file.path(pub, "config", "paths.yaml"))
user_root  <- paths_yaml$roots$hpc
resolve <- function(token) {
  val <- paths_yaml$sources[[token]]
  if (is.null(val)) stop(sprintf("Token not found in paths.yaml: %s", token))
  file.path(user_root, val)
}

contrast          <- args$contrast
joint_da_path     <- resolve(paste0("joint_da_results_", contrast))
# FLEX parity track only carries c1, c3, c4. Within-P3 sub-contrasts
# (c2a, c2b, c3a, c3b, c3c) have no FLEX equivalent. tryCatch so missing
# token degrades to joint-only without breaking prepare.
flex_da_path      <- tryCatch(
  resolve(paste0("flex_da_results_", contrast)),
  error = function(e) NULL
)
has_flex          <- !is.null(flex_da_path) && file.exists(flex_da_path)
counts_path       <- resolve("joint_nhood_counts")
l1p5_path         <- resolve("xenium_joint_l1p5")
lbridge_path      <- resolve("lbridge_vocabulary")
flex_l2s_map_path <- resolve("flex_l2s_lbridge_map")
flex_l2s_path     <- resolve("flex_nhood_l2s")
cat(sprintf("[prepare] has_flex   : %s\n", has_flex))

cat(sprintf("[prepare] contrast : %s\n", contrast))

ARTIFACT_LABELS <- c(
  "BMYO_artifact", "BMYO-NC", "Epithelial_artifact",
  "Stromal_artifact", "Immune_artifact", "Macrophage_artifact",
  "MastCell_artifact"
)
SIG_FDR        <- 0.1
PURITY_MIN     <- args$purity_min
MIN_L2S_NHOODS <- args$min_l2s_nhoods

# ── Lookups ────────────────────────────────────────────────────────────────────

cat("[prepare] Loading compartment map...\n")
l1p5 <- fread(l1p5_path, select = c("l1p5_label", "compartment", "is_artifact"))
compartment_map <- l1p5[is_artifact == FALSE,
  .(compartment = names(which.max(table(compartment)))), by = l1p5_label]

cat("[prepare] Loading LBridge vocabulary...\n")
lbridge_vocab <- fread(lbridge_path)
flex_l2s_map  <- fread(flex_l2s_map_path,
                       select = c("flex_l2s_label", "lbridge_family", "display_l20"))

# Canonical display names for joint (Xenium L1.5) and FLEX L2S rows.
# display_l15 / display_l20 columns in the vocabulary align with fig3d panel labels.
joint_display <- lbridge_vocab[xenium_l15_label != "NA" & !is.na(xenium_l15_label) &
                                display_l15 != "NA" & !is.na(display_l15),
                               setNames(display_l15, xenium_l15_label)]
joint_display["BMYO"] <- "BMYO"

flex_display <- c(
  lbridge_vocab[ihbca_l20_label != "NA" & !is.na(ihbca_l20_label) &
                 display_l20 != "NA" & !is.na(display_l20),
                setNames(display_l20, ihbca_l20_label)],
  setNames(flex_l2s_map$display_l20, flex_l2s_map$flex_l2s_label)
)

l15_lbridge <- lbridge_vocab[
  !is.na(xenium_l15_label) & xenium_l15_label != "NA",
  .(l1p5_label = xenium_l15_label, lbridge_family)
][, .SD[1], by = l1p5_label]
l15_lbridge <- rbind(l15_lbridge,
                     data.table(l1p5_label = "BMYO", lbridge_family = "BMYO"))

l2s_lbridge <- lbridge_vocab[
  !is.na(ihbca_l20_label) & ihbca_l20_label != "NA",
  .(l2s_label = ihbca_l20_label, lbridge_family)
][, .SD[1], by = l2s_label]

# Explicit family assignments for FLEX Track-C labels not in the LBridge L2.0 column.
# These arise because FLEX-only Track-C annotation has higher gene-space resolution
# than Xenium-panel features; they compress to different families at L2S level.
L2S_FAMILY_OVERRIDE <- setNames(flex_l2s_map$lbridge_family,
                                 flex_l2s_map$flex_l2s_label)

# ── Platform classification ────────────────────────────────────────────────────

cat("[prepare] Loading nhood counts (platform classification)...\n")
nc <- fread(counts_path)
xenium_cols <- grep("_xenium_", colnames(nc), value = TRUE)
flex_cols   <- setdiff(colnames(nc), c("Nhood", xenium_cols))
nc[, xenium_cells := rowSums(.SD), .SDcols = xenium_cols]
nc[, flex_cells   := rowSums(.SD), .SDcols = flex_cols]
nc[, platform     := fifelse(xenium_cells >= flex_cells, "Xenium", "FLEX")]
nc_sub <- nc[, .(Nhood, platform)]

# ── Joint layer ────────────────────────────────────────────────────────────────

cat("[prepare] Loading joint DA results...\n")
da_joint <- fread(joint_da_path)

joint_df <- da_joint |>
  as.data.table() |>
  filter(!nhood_label_l1p5 %in% ARTIFACT_LABELS,
         !is.na(nhood_label_l1p5), nhood_label_l1p5 != "NA",
         nhood_label_l1p5_frac >= PURITY_MIN) |>
  merge(nc_sub, by = "Nhood", all.x = TRUE) |>
  filter(!is.na(platform)) |>
  merge(l15_lbridge, by.x = "nhood_label_l1p5", by.y = "l1p5_label", all.x = TRUE) |>
  filter(!is.na(lbridge_family), lbridge_family != "NA") |>
  mutate(type_label = ifelse(nhood_label_l1p5 %in% names(joint_display),
                             joint_display[nhood_label_l1p5],
                             nhood_label_l1p5),
         layer = "joint") |>
  select(Nhood, logFC, SpatialFDR, type_label, lbridge_family, platform, layer)

cat(sprintf("[prepare] Joint: %d nhoods, %d types\n",
    nrow(joint_df), n_distinct(joint_df$type_label)))

# ── FLEX L2S layer ─────────────────────────────────────────────────────────────

# FLEX-only nhoods annotated at L2S resolution by majority-vote from
# Track-C (iHBCA L2.0) labels.  Purity filter retains nhoods with a
# dominant label fraction >= PURITY_MIN.
#
# Skipped entirely when has_flex is FALSE (within-P3 sub-contrasts c2a/c2b/
# c3a/c3b/c3c). flex_df is emitted as an empty 7-column data.table so the
# downstream bind_rows + lbridge_stats pipeline runs unchanged.

if (has_flex) {

cat("[prepare] Loading FLEX DA results...\n")
da_flex <- fread(flex_da_path)

cat("[prepare] Loading FLEX per-nhood L2S annotation...\n")
nhood_l2s <- fread(flex_l2s_path,
                   select = c("nhood_id", "dominant_label", "label_purity"))

flex_df <- da_flex |>
  as.data.table() |>
  filter(!nhood_label_l1p5 %in% ARTIFACT_LABELS,
         nhood_label_l1p5_frac >= PURITY_MIN) |>
  merge(nhood_l2s, by.x = "Nhood", by.y = "nhood_id", all.x = TRUE) |>
  filter(!is.na(dominant_label), dominant_label != "NA",
         !grepl("ARTIFACT", dominant_label, ignore.case = TRUE),
         !is.na(label_purity), label_purity >= PURITY_MIN) |>
  merge(l2s_lbridge, by.x = "dominant_label", by.y = "l2s_label", all.x = TRUE) |>
  mutate(
    lbridge_family = case_when(
      !is.na(lbridge_family) & lbridge_family != "NA" ~ lbridge_family,
      dominant_label %in% names(L2S_FAMILY_OVERRIDE) ~
        L2S_FAMILY_OVERRIDE[dominant_label],
      TRUE ~ l15_lbridge$lbridge_family[
               match(nhood_label_l1p5, l15_lbridge$l1p5_label)]
    )
  ) |>
  filter(!is.na(lbridge_family), lbridge_family != "NA") |>
  mutate(type_label = ifelse(dominant_label %in% names(flex_display),
                             flex_display[dominant_label],
                             dominant_label),
         layer = "flex_l2s", platform = "FLEX") |>
  group_by(type_label) |>
  filter(n() >= MIN_L2S_NHOODS) |>
  ungroup() |>
  select(Nhood, logFC, SpatialFDR, type_label, lbridge_family, platform, layer)

# Suppress FLEX L2S rows whose display label matches the joint row in the same
# family (1:1 identity cases where L2S resolution adds no additional information).
joint_labels_by_family <- joint_df |>
  distinct(type_label, lbridge_family) |>
  rename(joint_type_label = type_label)
flex_df <- flex_df |>
  left_join(joint_labels_by_family, by = "lbridge_family") |>
  filter(is.na(joint_type_label) | type_label != joint_type_label) |>
  select(-joint_type_label)

cat(sprintf("[prepare] FLEX L2S: %d nhoods, %d types\n",
    nrow(flex_df), n_distinct(flex_df$type_label)))

} else {
  cat("[prepare] FLEX L2S: skipped (no flex_da_results_*. token for this contrast)\n")
  flex_df <- data.table(
    Nhood          = integer(0),
    logFC          = double(0),
    SpatialFDR     = double(0),
    type_label     = character(0),
    lbridge_family = character(0),
    platform       = character(0),
    layer          = character(0)
  )
}

# ── LBridge family ordering ─────────────────────────────────────────────────────

# Families are ordered by compartment (Epithelial / Immune / Stromal), then by
# median signed logFC of the joint rows within each family.

plot_df_full <- bind_rows(joint_df, flex_df) |>
  mutate(sig = SpatialFDR < SIG_FDR,
         platform = factor(platform, levels = c("FLEX", "Xenium")))

type_stats <- plot_df_full |>
  group_by(type_label, layer, lbridge_family) |>
  summarise(med_lfc = median(if (any(sig)) logFC[sig] else logFC, na.rm = TRUE),
            .groups = "drop")

family_compartment <- lbridge_vocab[
  , .(lbridge_compartment = lbridge_compartment[1]), by = lbridge_family]

# Fixed compartment-stack order, constant across contrasts. ggplot's y-axis
# convention plots factor levels bottom-to-top, so a HIGHER rank value puts
# the compartment HIGHER on the panel. Target visual stack (top → bottom):
# Epithelial → Stromal → Immune. Map: Immune=1 (bottom), Stromal=2,
# Epithelial=3 (top).
COMPARTMENT_STACK_ORDER <- c("Immune" = 1, "Stromal" = 2, "Epithelial" = 3)

lbridge_stats <- type_stats |>
  filter(layer == "joint") |>
  group_by(lbridge_family) |>
  summarise(lbridge_med = median(med_lfc, na.rm = TRUE), .groups = "drop") |>
  left_join(as.data.frame(family_compartment), by = "lbridge_family") |>
  mutate(lbridge_compartment = coalesce(lbridge_compartment, "Unknown")) |>
  mutate(compartment_rank = unname(
           COMPARTMENT_STACK_ORDER[as.character(lbridge_compartment)])) |>
  # Unknown / unmapped compartments fall to rank 0 so they stack below Immune.
  mutate(compartment_rank = ifelse(is.na(compartment_rank), 0L, compartment_rank)) |>
  arrange(compartment_rank, lbridge_med) |>
  mutate(lbridge_rank = row_number())

# ── Save ────────────────────────────────────────────────────────────────────────

out_path <- file.path(out_dir, sprintf("plot_data_%s.rds", contrast))
saveRDS(list(
  joint_df      = joint_df,
  flex_df       = flex_df,
  lbridge_stats = lbridge_stats,
  SIG_FDR       = SIG_FDR,
  contrast      = contrast
), out_path)
cat(sprintf("[prepare] Saved: %s\n", out_path))
