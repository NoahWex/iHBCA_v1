#!/usr/bin/env Rscript
# plot_panelB_split_forest_absolute.R
#
# VISUAL AID — NOT canonical. Same layout as plot_panelB_split_forest.R, but
# right-half violin + right-side forest dots use the DERIVED absolute parity
# effect in BR1 (logFC_AR + logFC_int) instead of the interaction logFC.
# NhoodGroup carving is still from the parity_x_HR_BRCA1 (interaction) F.1.
# Output: panelB_split_forest_absolute.pdf

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(tidyr); library(yaml)
})
`%||%` <- function(a, b) if (is.null(a)) b else a

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "plot_panelB_split_forest_absolute.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths

CONTRASTS <- c("parity_in_AR", "parity_x_HR_BRCA1")
contrast_short <- c(parity_in_AR = "AR", parity_x_HR_BRCA1 = "BR1 (abs)")
SIDE <- c(parity_in_AR = -1L, parity_x_HR_BRCA1 = +1L)
DODGE <- 0.22

RANK_PALETTE <- c(
  "1"  = "#332288", "2"  = "#117733", "3"  = "#44AA99", "4"  = "#88CCEE",
  "5"  = "#DDCC77", "6"  = "#CC6677", "7"  = "#AA4499", "8"  = "#882255",
  "9"  = "#999933", "10" = "#661100"
)

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

# Artifact filter
labels_path <- paths$inputs$labels
artifact_l2s <- character(0)
if (file.exists(labels_path)) {
  art <- read_csv(labels_path, show_col_types = FALSE,
                   col_select = c("compartment","label","is_artifact")) %>%
    mutate(is_artifact = as.logical(is_artifact)) %>%
    group_by(compartment, label) %>%
    summarise(p = mean(is_artifact, na.rm = TRUE), .groups = "drop") %>%
    filter(p >= 0.5) %>%
    mutate(L2_joint = paste(compartment, label, sep = "::"))
  artifact_l2s <- art$L2_joint
  cat(sprintf("Dropping %d artifact L2s\n", length(artifact_l2s)))
}

# Load AR contrast first to get the absolute baseline
da_AR <- read_csv(file.path(paths$outputs$stageD, "parity_in_AR", "da_results.csv"),
                   show_col_types = FALSE) %>%
  select(Nhood, logFC_AR = logFC, compartment, label)

# Per-contrast loop
all_nhoods <- list()
for (cn in CONTRASTS) {
  da_csv <- file.path(paths$outputs$stageD, cn, "da_results.csv")
  ng_csv <- file.path(paths$outputs$stageF1, cn, "nhood_groups.csv")
  if (!file.exists(da_csv) || !file.exists(ng_csv)) next
  da <- read_csv(da_csv, show_col_types = FALSE)
  ng <- read_csv(ng_csv, show_col_types = FALSE)
  if ("NhoodGroup" %in% colnames(ng)) ng$NhoodGroup <- as.character(ng$NhoodGroup)
  if (!"compartment" %in% colnames(da) && "L2_compartment" %in% colnames(da))
    da$compartment <- da$L2_compartment
  if (!"label" %in% colnames(da) && "L2_label" %in% colnames(da))
    da$label <- da$L2_label
  da$L2_joint <- paste(da$compartment, da$label, sep = "::")
  da$contrast <- cn

  # KEY SUBSTITUTION: for the BRCA1 contrast, replace logFC with absolute
  # parity-in-BR1 = logFC_AR + logFC_int.
  if (cn == "parity_x_HR_BRCA1") {
    da <- da %>%
      left_join(da_AR %>% select(Nhood, logFC_AR), by = "Nhood") %>%
      mutate(logFC_int = logFC,
             logFC = logFC_AR + logFC_int)
    cat(sprintf("Substituted BR1 logFC = AR + interaction. New median: %.3f\n",
                median(da$logFC, na.rm = TRUE)))
  }

  ng_key <- intersect(c("Nhood","NhoodIndex","nhood_index"), colnames(ng))[1]
  da_key <- intersect(c("Nhood","NhoodIndex","nhood_index"), colnames(da))[1]
  da <- da %>%
    left_join(ng %>% select(all_of(ng_key), NhoodGroup),
              by = setNames(ng_key, da_key)) %>%
    filter(!L2_joint %in% artifact_l2s)
  da$NhoodGroup <- as.character(da$NhoodGroup)
  da <- da %>% left_join(
    lookup %>% filter(contrast == cn) %>%
      select(parent_L2_joint, NhoodGroup, group_dir,
             n_nhoods_in_group, size_rank_within_L2, name_short, name_full),
    by = c("L2_joint" = "parent_L2_joint", "NhoodGroup" = "NhoodGroup"))
  all_nhoods[[cn]] <- da
}
df <- bind_rows(all_nhoods)
cat(sprintf("Per-nhood rows: %d\n", nrow(df)))

# L2 ordering
order_l2 <- df %>%
  filter(contrast == "parity_in_AR") %>%
  group_by(compartment, L2_joint) %>%
  summarise(med = median(logFC, na.rm = TRUE), .groups = "drop") %>%
  arrange(compartment, desc(med))
extra_l2 <- setdiff(unique(df$L2_joint), order_l2$L2_joint)
if (length(extra_l2) > 0) {
  ex_rows <- df %>% filter(L2_joint %in% extra_l2) %>%
    distinct(compartment, L2_joint) %>% mutate(med = 0)
  order_l2 <- bind_rows(order_l2, ex_rows) %>% arrange(compartment, desc(med))
}
df$compartment <- factor(df$compartment, levels = c("epi","imm","str"))
df$L2_joint <- factor(df$L2_joint, levels = order_l2$L2_joint)

# Half-violin polygon builder
build_half <- function(d, side, max_width = 0.42) {
  empty <- data.frame(x_off = numeric(0), y = numeric(0))
  if (nrow(d) < 5) return(empty)
  dn <- density(d$logFC, n = 128, na.rm = TRUE,
                from = min(d$logFC, na.rm = TRUE),
                to   = max(d$logFC, na.rm = TRUE))
  scale <- max(dn$y, na.rm = TRUE)
  if (!is.finite(scale) || scale <= 0) return(empty)
  width <- (dn$y / scale) * max_width
  data.frame(x_off = c(side * width, rep(0, length(dn$x))),
             y     = c(dn$x, rev(dn$x)))
}

violin_df <- df %>%
  group_by(compartment, L2_joint, contrast) %>%
  group_modify(~ build_half(.x, SIDE[as.character(.y$contrast)])) %>%
  ungroup()

l2_pos <- order_l2 %>% ungroup() %>%
  arrange(compartment, desc(med)) %>%
  mutate(x_pos = row_number()) %>%
  select(compartment, L2_joint, x_pos)
violin_df <- violin_df %>%
  left_join(l2_pos, by = c("compartment","L2_joint")) %>%
  mutate(x = x_pos + x_off)

forest_df <- df %>%
  filter(!is.na(size_rank_within_L2)) %>%
  group_by(compartment, L2_joint, contrast, NhoodGroup, name_short,
           group_dir, n_nhoods_in_group, size_rank_within_L2) %>%
  summarise(median = median(logFC, na.rm = TRUE),
            q25    = quantile(logFC, 0.25, na.rm = TRUE),
            q75    = quantile(logFC, 0.75, na.rm = TRUE),
            .groups = "drop") %>%
  left_join(l2_pos, by = c("compartment","L2_joint")) %>%
  mutate(x = x_pos + ifelse(contrast == "parity_in_AR", -DODGE, DODGE),
         rank_str = ifelse(size_rank_within_L2 %in% 1:10,
                            as.character(size_rank_within_L2), "10"))

breaks_df <- l2_pos %>% mutate(label = sub("^[^:]*::", "", L2_joint))

p <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.3) +
  geom_polygon(data = violin_df,
               aes(x = x, y = y, group = interaction(L2_joint, contrast)),
               fill = NA, color = "grey35", linewidth = 0.3) +
  geom_linerange(data = forest_df,
                 aes(x = x, ymin = q25, ymax = q75, color = rank_str),
                 linewidth = 0.45, alpha = 0.9) +
  geom_point(data = forest_df,
             aes(x = x, y = median, color = rank_str,
                 size = n_nhoods_in_group),
             shape = 16, alpha = 0.95) +
  scale_x_continuous(breaks = breaks_df$x_pos, labels = breaks_df$label,
                     expand = expansion(add = 0.05),
                     guide = guide_axis(check.overlap = FALSE)) +
  coord_cartesian(clip = "off") +
  scale_color_manual("NhoodGroup rank (within L2)",
                      values = RANK_PALETTE,
                      breaks = as.character(1:10),
                      na.value = "grey60") +
  scale_size_continuous("NhoodGroup size (n_nhoods)", range = c(0.8, 3.2)) +
  facet_grid(. ~ compartment, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "per-nhood logFC",
        caption = "VISUAL AID: left = parity in AR (baseline), right = parity in BR1 (DERIVED absolute = AR + interaction)") +
  theme_classic(base_size = 8) +
  theme(strip.text = element_text(face = "bold", size = 8.5),
        strip.background = element_blank(),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.key.size = unit(0.3, "cm"),
        axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5),
        axis.title.y = element_text(size = 8),
        plot.caption = element_text(size = 7, color = "grey30"),
        panel.grid = element_blank(),
        panel.background = element_blank())

out_dir <- file.path(paths$outputs$reports, "figures")
out_pdf <- file.path(out_dir, "panelB_split_forest_absolute.pdf")
ggsave(out_pdf, p,
       width = max(7, 0.32 * nrow(l2_pos) + 1.5), height = 6.5, limitsize = FALSE)
cat(sprintf("Wrote %s\n", out_pdf))
