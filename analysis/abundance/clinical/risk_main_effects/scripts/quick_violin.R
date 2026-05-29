suppressPackageStartupMessages({ library(dplyr); library(readr); library(ggplot2); library(tidyr) })

df <- read_csv("/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/da_pipeline/inquiries/parity_x_risk_20260504/outputs/per_nhood_BR1_absolute_derived.csv",
                show_col_types = FALSE)

# Drop artifacts at L2 grain (same filter as Panel B v11)
labels_path <- "/share/crsp/lab/dalawson/nwechter/iHBCA_publication/publication/analysis/annotation/labels_full.csv"
art_l2 <- read_csv(labels_path, show_col_types = FALSE,
                    col_select = c("compartment","label","is_artifact")) %>%
  mutate(is_artifact = as.logical(is_artifact)) %>%
  group_by(compartment, label) %>%
  summarise(p = mean(is_artifact, na.rm = TRUE), .groups = "drop") %>%
  filter(p >= 0.5) %>% mutate(L2_joint = paste(compartment, label, sep = "::")) %>%
  pull(L2_joint)
df <- df %>% filter(!L2_joint %in% art_l2)

long <- df %>%
  pivot_longer(c(logFC_AR, logFC_int, logFC_BR1_abs),
               names_to = "contrast", values_to = "logFC") %>%
  mutate(contrast = factor(contrast,
                            levels = c("logFC_AR", "logFC_int", "logFC_BR1_abs"),
                            labels = c("parity in AR\n(baseline)",
                                       "parity x BRCA1\n(modification)",
                                       "parity in BR1\n(derived absolute)")))

p <- ggplot(long, aes(x = contrast, y = logFC, fill = contrast)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  geom_violin(alpha = 0.55, color = "grey25", linewidth = 0.3) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.8,
               color = "grey25", linewidth = 0.3, fill = "white") +
  scale_fill_manual(values = c("#1f77b4", "grey60", "#d62728"), guide = "none") +
  labs(x = NULL, y = "per-nhood logFC",
        caption = sprintf("n = %d nhoods (artifact L2s dropped). BR1_abs = AR + interaction.",
                           nrow(df))) +
  theme_classic(base_size = 9) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 8))

out <- "/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/da_pipeline/inquiries/parity_x_risk_20260504/reports/figures/quick_violin_AR_int_BR1abs.pdf"
ggsave(out, p, width = 5.5, height = 4.5)
cat(sprintf("Wrote %s\n", out))
