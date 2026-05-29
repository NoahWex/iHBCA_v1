#!/usr/bin/env Rscript
# Quick viability check: derive absolute parity-in-BR1 per nhood from
# parity_in_AR + parity_x_HR_BRCA1 Stage D outputs.
suppressPackageStartupMessages({
  library(dplyr); library(readr)
})

base <- "/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/da_pipeline/inquiries/parity_x_risk_20260504/outputs/stageD_da_results"
ar <- read_csv(file.path(base, "parity_in_AR/da_results.csv"),
                show_col_types = FALSE)
int <- read_csv(file.path(base, "parity_x_HR_BRCA1/da_results.csv"),
                 show_col_types = FALSE)

cat(sprintf("AR  rows=%d  cols=%s\n", nrow(ar),  paste(colnames(ar)[1:6], collapse=",")))
cat(sprintf("INT rows=%d  cols=%s\n\n", nrow(int), paste(colnames(int)[1:6], collapse=",")))

# Join by Nhood index
both <- inner_join(ar  %>% select(Nhood, logFC_AR  = logFC,  SpatialFDR_AR  = SpatialFDR,
                                   compartment, label),
                   int %>% select(Nhood, logFC_int = logFC,  SpatialFDR_int = SpatialFDR),
                   by = "Nhood")
both$L2_joint <- paste(both$compartment, both$label, sep = "::")
both$logFC_BR1_abs <- both$logFC_AR + both$logFC_int

cat(sprintf("Matched nhoods: %d / %d\n\n", nrow(both),
            min(nrow(ar), nrow(int))))

# Distribution summaries
cat("=== Per-nhood logFC summaries ===\n")
sumstat <- function(x, name) {
  q <- quantile(x, c(0.01, 0.25, 0.5, 0.75, 0.99), na.rm = TRUE)
  cat(sprintf("%-22s  min=%6.2f  q01=%6.2f  q25=%6.2f  med=%6.2f  q75=%6.2f  q99=%6.2f  max=%6.2f  sd=%5.2f\n",
              name, min(x, na.rm = TRUE), q[1], q[2], q[3], q[4], q[5],
              max(x, na.rm = TRUE), sd(x, na.rm = TRUE)))
}
sumstat(both$logFC_AR,      "parity_in_AR")
sumstat(both$logFC_int,     "parity_x_HR_BRCA1 (mod)")
sumstat(both$logFC_BR1_abs, "parity_in_BR1 (derived)")

# Concordance: where does AR != BR1_abs (i.e., BRCA1 modifies parity effect)
abs_delta <- abs(both$logFC_BR1_abs - both$logFC_AR)
cat(sprintf("\nMod magnitude |L_int|: median=%.3f, q90=%.3f, q99=%.3f\n",
            median(abs_delta), quantile(abs_delta, 0.90, na.rm=TRUE),
            quantile(abs_delta, 0.99, na.rm=TRUE)))

# Sign-flips: AR > 0 and BR1_abs < 0 (or reverse)
flip <- (sign(both$logFC_AR) != sign(both$logFC_BR1_abs)) &
         abs(both$logFC_AR) > 0.1 & abs(both$logFC_BR1_abs) > 0.1
cat(sprintf("Sign-flips between AR and BR1_abs (both |logFC|>0.1): %d (%.1f%%)\n",
            sum(flip, na.rm = TRUE),
            100 * mean(flip, na.rm = TRUE)))

# Per-L2 medians: how does the absolute look at L2 grain?
per_l2 <- both %>%
  group_by(compartment, L2_joint) %>%
  summarise(n = dplyr::n(),
            med_AR  = median(logFC_AR, na.rm = TRUE),
            med_int = median(logFC_int, na.rm = TRUE),
            med_BR1_abs = median(logFC_BR1_abs, na.rm = TRUE),
            sig_AR  = sum(SpatialFDR_AR  < 0.05, na.rm = TRUE),
            sig_int = sum(SpatialFDR_int < 0.05, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(desc(abs(med_BR1_abs - med_AR)))

cat("\n=== Top 20 L2s by |BR1_abs − AR| (largest BRCA1 modifications) ===\n")
print(head(per_l2, 20), n = 20)

cat("\n=== L2s where direction reverses (med_AR vs med_BR1_abs) ===\n")
flip_l2 <- per_l2 %>%
  filter(sign(med_AR) != sign(med_BR1_abs),
         abs(med_AR) > 0.05, abs(med_BR1_abs) > 0.05)
print(flip_l2, n = Inf)

# Save the per-nhood derived table for downstream use
out_csv <- "/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/da_pipeline/inquiries/parity_x_risk_20260504/outputs/per_nhood_BR1_absolute_derived.csv"
write_csv(both %>% select(Nhood, compartment, label, L2_joint,
                           logFC_AR, logFC_int, logFC_BR1_abs,
                           SpatialFDR_AR, SpatialFDR_int),
          out_csv)
cat(sprintf("\nWrote derived per-nhood table: %s\n", out_csv))
