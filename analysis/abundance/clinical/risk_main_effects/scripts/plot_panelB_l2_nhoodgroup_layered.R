#!/usr/bin/env Rscript
# plot_panelB_l2_nhoodgroup_layered.R
#
# Per-L2 panel: overlapping NhoodGroup-grain density distributions of
# per-nhood logFC, filtered to limma-viable populations (L2s with at least
# one viable NhoodGroup, n_nhoods >= 10).
#
# Display unit = NhoodGroup density curve over its members' logFC values.
# No per-nhood beeswarm, no backdrop. Color = NhoodGroup direction (sig_up
# red / sig_dn blue / mixed grey). Density area scaled by n_nhoods_in_group.
#
# Layout:
#   y = composite (L2_joint x contrast), each L2 stacked across contrasts
#   x = logFC
#   facet rows = compartment, ordered by parity_in_AR median lfc within compartment
#
# Inputs:
#   stageD/<contrast>/da_results.csv          per-nhood logFC + SpatialFDR
#   stageF1_nhoodgroups/<contrast>/nhood_groups.csv   nhood -> NhoodGroup map
#   nhoodgroup_renaming/lookup.csv            viable group table + direction

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(tidyr); library(yaml)
})
have_ridges <- requireNamespace("ggridges", quietly = TRUE)

`%||%` <- function(a, b) if (is.null(a)) b else a

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "plot_panelB_l2_nhoodgroup_layered.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths

CONTRASTS <- c("parity_in_AR", "parity_x_HR_BRCA1")
contrast_short <- c(parity_in_AR = "AR", parity_x_HR_BRCA1 = "BR1")

SIG_PALETTE <- c(sig_up = "#d62728", sig_dn = "#1f77b4", mixed = "grey70")

# ---- Load lookup (viable NhoodGroups only, restrict to target contrasts) ----
lookup <- read_csv(file.path(paths$inquiry_root, "outputs",
                              "nhoodgroup_renaming", "lookup.csv"),
                    show_col_types = FALSE) %>%
  filter(viable, contrast %in% CONTRASTS) %>%
  mutate(NhoodGroup = as.character(NhoodGroup),
         group_dir = case_when(
           group_pct_up >= 60 ~ "sig_up",
           group_pct_up <= 40 ~ "sig_dn",
           TRUE                ~ "mixed"
         )) %>%
  select(contrast, parent_L2_joint, parent_compartment, NhoodGroup,
         n_nhoods_in_group, group_med_lfc, group_pct_up, group_dir,
         size_rank_within_L2, name_short, name_full)

cat(sprintf("Viable NhoodGroups in target contrasts: %d (across %d L2s)\n",
            nrow(lookup), length(unique(lookup$parent_L2_joint))))

viable_l2s <- unique(lookup$parent_L2_joint)

# ---- Load per-nhood DA + group membership; restrict to viable groups ----
all_nhoods <- list()
for (cn in CONTRASTS) {
  da_csv <- file.path(paths$outputs$stageD, cn, "da_results.csv")
  ng_csv <- file.path(paths$outputs$stageF1, cn, "nhood_groups.csv")
  if (!file.exists(da_csv) || !file.exists(ng_csv)) {
    cat(sprintf("MISSING substrate for %s\n", cn)); next
  }
  da <- read_csv(da_csv, show_col_types = FALSE)
  ng <- read_csv(ng_csv, show_col_types = FALSE)
  if ("NhoodGroup" %in% colnames(ng)) ng$NhoodGroup <- as.character(ng$NhoodGroup)
  if (!"compartment" %in% colnames(da) && "L2_compartment" %in% colnames(da))
    da$compartment <- da$L2_compartment
  if (!"label" %in% colnames(da) && "L2_label" %in% colnames(da))
    da$label <- da$L2_label
  da$L2_joint <- paste(da$compartment, da$label, sep = "::")
  da$contrast <- cn

  ng_key <- intersect(c("Nhood", "NhoodIndex", "nhood_index"), colnames(ng))[1]
  da_key <- intersect(c("Nhood", "NhoodIndex", "nhood_index"), colnames(da))[1]
  if (is.na(ng_key) || is.na(da_key)) {
    cat(sprintf("WARN no matching NhoodGroup key for %s\n", cn)); next
  }
  da <- da %>%
    left_join(ng %>% select(all_of(ng_key), NhoodGroup),
              by = setNames(ng_key, da_key)) %>%
    filter(!is.na(NhoodGroup))
  da$NhoodGroup <- as.character(da$NhoodGroup)

  da <- da %>% filter(L2_joint %in% viable_l2s)

  da <- da %>% inner_join(
    lookup %>% filter(contrast == cn) %>%
      select(parent_L2_joint, NhoodGroup, group_dir,
             n_nhoods_in_group, size_rank_within_L2, name_short, name_full),
    by = c("L2_joint" = "parent_L2_joint", "NhoodGroup" = "NhoodGroup")
  )
  all_nhoods[[cn]] <- da
}
df <- bind_rows(all_nhoods)
cat(sprintf("Total per-nhood rows (viable groups, target contrasts): %d\n",
            nrow(df)))
cat(sprintf("Per-NhoodGroup rows: %d\n",
            n_distinct(df$contrast, df$L2_joint, df$NhoodGroup)))

# ---- Composite y: L2 stacked across contrasts ----
order_l2 <- df %>%
  filter(contrast == "parity_in_AR") %>%
  group_by(compartment, L2_joint) %>%
  summarise(med = median(logFC, na.rm = TRUE), .groups = "drop") %>%
  arrange(compartment, desc(med))

# Include any L2 that's viable in BR1 but not AR (so both contrasts represented)
extra_l2 <- setdiff(unique(df$L2_joint), order_l2$L2_joint)
if (length(extra_l2) > 0) {
  ex_rows <- df %>% filter(L2_joint %in% extra_l2) %>%
    distinct(compartment, L2_joint) %>%
    mutate(med = 0)
  order_l2 <- bind_rows(order_l2, ex_rows) %>%
    arrange(compartment, desc(med))
}

df$contrast_short <- contrast_short[as.character(df$contrast)]
df$L2_contrast <- paste0(df$L2_joint, "  [", df$contrast_short, "]")
lc_levels <- unlist(lapply(order_l2$L2_joint, function(l)
  paste0(l, "  [", contrast_short[CONTRASTS], "]")))
df$L2_contrast <- factor(df$L2_contrast, levels = rev(lc_levels))
df$compartment <- factor(df$compartment, levels = c("epi", "imm", "str"))

# Drawing order: smallest groups first so largest paint last (largest most visible
# at high alpha; smaller groups still legible through the alpha overlap).
df <- df %>% arrange(contrast, L2_joint, n_nhoods_in_group)

# ---- Plot ----
if (!have_ridges) {
  stop("ggridges not available; cannot render overlapping density panel.")
}

p <- ggplot(df,
            aes(x = logFC, y = L2_contrast,
                group = interaction(L2_contrast, NhoodGroup),
                fill = group_dir, height = after_stat(density))) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.3) +
  ggridges::geom_density_ridges(
    stat = "density", scale = 1.0, alpha = 0.55, color = "white",
    linewidth = 0.18, rel_min_height = 0.01) +
  scale_fill_manual("NhoodGroup direction",
                    values = SIG_PALETTE,
                    breaks = c("sig_up", "sig_dn", "mixed"),
                    labels = c("sig up", "sig down", "mixed")) +
  facet_grid(compartment ~ ., scales = "free_y", space = "free_y") +
  labs(x = "per-nhood logFC", y = NULL) +
  theme_classic(base_size = 8) +
  theme(strip.text = element_text(face = "bold", size = 8.5),
        strip.background = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(0.3, "cm"),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 8),
        panel.grid = element_blank(),
        panel.background = element_blank())

# ---- Render ----
out_dir <- file.path(paths$outputs$reports, "figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out_pdf <- file.path(out_dir, "panelB_l2_nhoodgroup_layered.pdf")

n_rows <- length(unique(df$L2_contrast))
ggsave(out_pdf, p,
       width  = 7,
       height = max(5, 0.18 * n_rows + 1.5),
       limitsize = FALSE)
cat(sprintf("Wrote %s (%d L2_contrast rows)\n", out_pdf, n_rows))

# ---- Sidecar summary ----
ng_summary <- df %>%
  group_by(contrast, compartment, L2_joint, NhoodGroup,
           name_short, name_full, n_nhoods_in_group, size_rank_within_L2,
           group_dir) %>%
  summarise(median_logFC = median(logFC, na.rm = TRUE),
            q25 = quantile(logFC, 0.25, na.rm = TRUE),
            q75 = quantile(logFC, 0.75, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(contrast, compartment, L2_joint, size_rank_within_L2)

write_csv(ng_summary,
          file.path(out_dir, "panelB_nhoodgroup_summary.csv"))
cat(sprintf("Wrote panelB_nhoodgroup_summary.csv (%d rows)\n", nrow(ng_summary)))
