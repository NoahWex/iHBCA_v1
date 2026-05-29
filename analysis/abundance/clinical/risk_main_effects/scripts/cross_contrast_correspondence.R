#!/usr/bin/env Rscript
# cross_contrast_correspondence.R
# Cross-contrast nhood-membership correspondence over all (L2 + viable
# NhoodGroup) entities. Produces:
#   - correspondence_long.csv  long pairwise table
#   - correspondence_matrix.csv square Jaccard matrix
#   - reports/figures/correspondence_clustered_global.pdf
#   - reports/figures/correspondence_per_parent_L2.pdf
#
# Entities:
#   - (L2) one per cell-type L2 ident, nhood set := all nhoods with that L2
#     label in the milo (read from any F.1 file, since milo is contrast-shared)
#   - (contrast, NhoodGroup_renamed) one per viable group (n_nhoods >= 10)
#
# Distance: 1 - Jaccard(nhood_set_i, nhood_set_j)

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(ggplot2);
  library(Matrix); library(stringr)
})

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "cross_contrast_correspondence.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths
fig_dir <- file.path(paths$outputs$reports, "figures")
ensure_dir(fig_dir)
out_dir <- paths$inquiry_root

ALL_CONTRASTS <- c("parity_in_AR", "parity_x_HR_BRCA1",
                   "parity_x_HR_sporadic", "parity_x_HR_BRCA2")
VIABLE_MIN_NHOODS <- 10

# ------------------------------------------------------------------------
# 1. Load per-contrast nhood->NhoodGroup map; build entity registry
# ------------------------------------------------------------------------
nhood_maps <- list()
for (cname in ALL_CONTRASTS) {
  csv <- file.path(paths$outputs$stageF1, cname, "nhood_groups.csv")
  if (!file.exists(csv)) next
  nhood_maps[[cname]] <- read_csv(csv, show_col_types = FALSE,
                                    col_types = cols(.default = "c")) %>%
    mutate(Nhood = as.integer(Nhood),
            n_nhoods_in_group = suppressWarnings(as.integer(n_nhoods_in_group)))
}
stopifnot(length(nhood_maps) >= 1)

# L2 entities (contrast-independent) from the first contrast's nhood file
ref <- nhood_maps[[1]]
l2_to_nhoods <- ref %>%
  filter(!is.na(L2_joint)) %>%
  group_by(L2_joint, compartment, label) %>%
  summarise(nhood_ids = list(sort(unique(Nhood))),
            n_nhoods = length(unique(Nhood)),
            .groups = "drop") %>%
  mutate(entity_id = paste0("L2::", L2_joint),
         entity_type = "L2",
         entity_contrast = NA_character_,
         entity_parent_L2 = L2_joint)

cat(sprintf("L2 entities: %d (median nhoods: %.0f)\n",
            nrow(l2_to_nhoods), median(l2_to_nhoods$n_nhoods)))

# NhoodGroup entities (per contrast, viable only)
ng_rows <- list()
for (cname in names(nhood_maps)) {
  nm <- nhood_maps[[cname]]
  ng <- nm %>%
    filter(!is.na(NhoodGroup_renamed),
           n_nhoods_in_group >= VIABLE_MIN_NHOODS) %>%
    group_by(NhoodGroup_renamed, parent_L2_joint, parent_compartment, parent_label) %>%
    summarise(nhood_ids = list(sort(unique(Nhood))),
              n_nhoods = length(unique(Nhood)),
              .groups = "drop") %>%
    mutate(entity_id = paste0(cname, "::", NhoodGroup_renamed),
           entity_type = "NhoodGroup",
           entity_contrast = cname,
           entity_parent_L2 = parent_L2_joint)
  ng_rows[[cname]] <- ng
}
ng_entities <- bind_rows(ng_rows)
cat(sprintf("NhoodGroup entities: %d (per contrast: %s)\n",
            nrow(ng_entities),
            paste(sprintf("%s=%d", names(nhood_maps),
                            sapply(ng_rows, nrow)), collapse=", ")))

# Combined entity registry
ent <- bind_rows(
  l2_to_nhoods  %>% select(entity_id, entity_type, entity_contrast,
                             entity_parent_L2, n_nhoods, nhood_ids),
  ng_entities    %>% select(entity_id, entity_type, entity_contrast,
                             entity_parent_L2, n_nhoods, nhood_ids)
)
cat(sprintf("Total entities: %d\n", nrow(ent)))

# ------------------------------------------------------------------------
# 2. Build sparse logical membership matrix (entity x nhood_id)
# ------------------------------------------------------------------------
all_nhoods <- sort(unique(unlist(ent$nhood_ids)))
n_ent <- nrow(ent)
n_nh <- length(all_nhoods)
nhood_idx <- setNames(seq_along(all_nhoods), as.character(all_nhoods))

i_idx <- integer(0); j_idx <- integer(0)
for (k in seq_len(n_ent)) {
  ids <- ent$nhood_ids[[k]]
  jj <- nhood_idx[as.character(ids)]
  i_idx <- c(i_idx, rep.int(k, length(jj)))
  j_idx <- c(j_idx, jj)
}
M <- sparseMatrix(i = i_idx, j = j_idx, x = 1, dims = c(n_ent, n_nh))
rownames(M) <- ent$entity_id

# Intersection counts via tcrossprod
INT <- tcrossprod(M)
size <- rowSums(M)
# Union: |A| + |B| - |A∩B|
SUM <- outer(size, size, FUN = "+")
UNION <- SUM - INT
J <- as.matrix(INT / UNION)
diag(J) <- 1
J[is.na(J)] <- 0
cat(sprintf("Jaccard matrix: %d x %d, density (>0): %.3f\n",
            nrow(J), ncol(J), mean(J[upper.tri(J)] > 0)))

# ------------------------------------------------------------------------
# 3. Long correspondence table (upper triangle only, exclude self)
# ------------------------------------------------------------------------
ut <- which(upper.tri(J), arr.ind = TRUE)
long <- tibble(
  entity_a = ent$entity_id[ut[, 1]],
  entity_b = ent$entity_id[ut[, 2]],
  type_a   = ent$entity_type[ut[, 1]],
  type_b   = ent$entity_type[ut[, 2]],
  contrast_a = ent$entity_contrast[ut[, 1]],
  contrast_b = ent$entity_contrast[ut[, 2]],
  parent_L2_a = ent$entity_parent_L2[ut[, 1]],
  parent_L2_b = ent$entity_parent_L2[ut[, 2]],
  n_a = ent$n_nhoods[ut[, 1]],
  n_b = ent$n_nhoods[ut[, 2]],
  n_inter = INT[ut],
  jaccard = J[ut]
)
long <- long %>% filter(jaccard > 0)
out_long <- file.path(out_dir, "outputs", "correspondence_long.csv")
write_csv(long, out_long)
cat(sprintf("Wrote: %s (%d non-zero pairs)\n", out_long, nrow(long)))

# Square matrix CSV
out_mat <- file.path(out_dir, "outputs", "correspondence_matrix.csv")
write.csv(round(J, 4), out_mat)
cat(sprintf("Wrote: %s\n", out_mat))

# ------------------------------------------------------------------------
# 4. Global clustered heatmap via base hclust + ggplot
# ------------------------------------------------------------------------
D <- as.dist(1 - J)
hc <- hclust(D, method = "average")
ord <- hc$order
ent_ord <- ent[ord, ]

# Long form for ggplot tile
J_df <- as.data.frame(J) %>%
  tibble::rownames_to_column("a") %>%
  pivot_longer(-a, names_to = "b", values_to = "jaccard") %>%
  mutate(a = factor(a, levels = ent_ord$entity_id),
          b = factor(b, levels = ent_ord$entity_id))

# Side annotation: contrast color
contrast_colors <- c(
  parity_in_AR = "#1B9E77",
  parity_x_HR_BRCA1 = "#D95F02",
  parity_x_HR_sporadic = "#7570B3",
  parity_x_HR_BRCA2 = "#E7298A"
)
ent_ord <- ent_ord %>%
  mutate(contrast_label = ifelse(is.na(entity_contrast), "L2_anchor", entity_contrast),
          contrast_label = factor(contrast_label,
                                    levels = c("L2_anchor", ALL_CONTRASTS)))

p_global <- ggplot(J_df, aes(x = b, y = a, fill = jaccard)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#08306B",
                        name = "Jaccard", limits = c(0, 1)) +
  theme_minimal(base_size = 5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 2),
        axis.text.y = element_text(size = 2),
        panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
        title = "Cross-contrast nhood-membership correspondence",
        subtitle = sprintf("%d entities (%d L2 + %d viable NhoodGroups across %d contrasts). Hierarchical clustering on 1-Jaccard, average linkage.",
                              nrow(ent), sum(ent$entity_type == "L2"),
                              sum(ent$entity_type == "NhoodGroup"), length(nhood_maps)))

out_pdf_global <- file.path(fig_dir, "correspondence_clustered_global.pdf")
ggsave(out_pdf_global, p_global, width = 16, height = 14)
cat(sprintf("Wrote: %s\n", out_pdf_global))

# ------------------------------------------------------------------------
# 5. Per-parent-L2 facet heatmap (clearer reading)
# ------------------------------------------------------------------------
# For each parent_L2: pull its L2 entity + all its viable NhoodGroups across
# contrasts. Compute submatrix Jaccard. Tile.
per_l2_rows <- list()
for (l2 in unique(ent$entity_parent_L2)) {
  if (is.na(l2)) next
  rows <- which(ent$entity_parent_L2 == l2)
  if (length(rows) < 2) next
  sub <- J[rows, rows, drop = FALSE]
  ids <- ent$entity_id[rows]
  contrasts <- ifelse(is.na(ent$entity_contrast[rows]), "L2_anchor",
                       ent$entity_contrast[rows])
  # Cluster rows within this L2
  hc_sub <- hclust(as.dist(1 - sub), method = "average")
  ord_sub <- hc_sub$order
  sub_df <- as.data.frame(sub[ord_sub, ord_sub, drop = FALSE]) %>%
    tibble::rownames_to_column("a") %>%
    pivot_longer(-a, names_to = "b", values_to = "jaccard") %>%
    mutate(a = factor(a, levels = ids[ord_sub]),
            b = factor(b, levels = ids[ord_sub]),
            parent_L2 = l2,
            contrast_a = contrasts[ord_sub][match(a, ids[ord_sub])],
            contrast_b = contrasts[ord_sub][match(b, ids[ord_sub])])
  per_l2_rows[[l2]] <- sub_df
}
per_l2 <- bind_rows(per_l2_rows)

# Limit to top L2s (most entities) for the facet PDF — readability
top_l2s <- ent %>%
  filter(!is.na(entity_parent_L2)) %>%
  count(entity_parent_L2, sort = TRUE) %>%
  filter(n >= 3) %>%   # at least L2 + 2 NhoodGroups
  pull(entity_parent_L2)
per_l2_plot <- per_l2 %>% filter(parent_L2 %in% top_l2s)

p_facet <- ggplot(per_l2_plot, aes(x = b, y = a, fill = jaccard)) +
  geom_tile() +
  facet_wrap(~ parent_L2, scales = "free", ncol = 4) +
  scale_fill_gradient(low = "white", high = "#08306B",
                        name = "Jaccard", limits = c(0, 1)) +
  theme_minimal(base_size = 5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4),
        axis.text.y = element_text(size = 4),
        panel.grid = element_blank(),
        strip.text = element_text(size = 7, face = "bold")) +
  labs(x = NULL, y = NULL,
        title = "Per-parent-L2 cross-contrast correspondence",
        subtitle = sprintf("%d L2s with >=3 entities (L2 anchor + viable NhoodGroups across contrasts). Within-L2 hierarchical clustering on 1-Jaccard.",
                              length(top_l2s)))

out_pdf_facet <- file.path(fig_dir, "correspondence_per_parent_L2.pdf")
n_l2 <- length(top_l2s)
n_cols <- 4
n_rows <- ceiling(n_l2 / n_cols)
ggsave(out_pdf_facet, p_facet, width = 16, height = max(4, n_rows * 3.5),
        limitsize = FALSE)
cat(sprintf("Wrote: %s (%d L2 facets)\n", out_pdf_facet, n_l2))

cat("\n=== cross-contrast correspondence done ===\n")
