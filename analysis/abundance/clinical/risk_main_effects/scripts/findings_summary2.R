suppressPackageStartupMessages({ library(dplyr); library(readr) })

cat("====================== PARITY x RISK ======================\n")
ft <- read_csv("/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/da_pipeline/inquiries/parity_x_risk_20260504/outputs/stageH_findings_table.csv",
                show_col_types = FALSE)
ftL <- ft %>% filter(row_type == "L2") %>% distinct(L2_joint, .keep_all = TRUE)

cat("\n--- Top 12 parity_in_AR by |median lfc| (n_sig >= 30) ---\n")
ftL %>% filter(parity_in_AR__n_sig >= 30, abs(parity_in_AR__med_lfc) > 0.10) %>%
  arrange(desc(abs(parity_in_AR__med_lfc))) %>%
  head(12) %>%
  transmute(L2_joint, ar_n = parity_in_AR__n_sig,
            ar_med = round(parity_in_AR__med_lfc, 2),
            ar_pctup = parity_in_AR__pct_up,
            br1_n = parity_x_HR_BRCA1__n_sig,
            br1_int = round(parity_x_HR_BRCA1__med_lfc, 2)) %>% as.data.frame() %>% print()

cat("\n--- Top 12 BR1 modifications by |interaction lfc| (n_sig >= 30) ---\n")
ftL %>% filter(parity_x_HR_BRCA1__n_sig >= 30, abs(parity_x_HR_BRCA1__med_lfc) > 0.20) %>%
  arrange(desc(abs(parity_x_HR_BRCA1__med_lfc))) %>%
  head(12) %>%
  transmute(L2_joint,
            ar_med = round(parity_in_AR__med_lfc, 2),
            br1_int = round(parity_x_HR_BRCA1__med_lfc, 2),
            br1_abs = round(parity_in_AR__med_lfc + parity_x_HR_BRCA1__med_lfc, 2),
            br1_n = parity_x_HR_BRCA1__n_sig) %>% as.data.frame() %>% print()

cat("\n--- L2s where BR1 reverses parity direction ---\n")
ftL %>% filter(parity_in_AR__n_sig >= 30, parity_x_HR_BRCA1__n_sig >= 5,
                sign(parity_in_AR__med_lfc) != sign(parity_x_HR_BRCA1__med_lfc),
                abs(parity_in_AR__med_lfc) > 0.10,
                abs(parity_x_HR_BRCA1__med_lfc) > 0.10) %>%
  transmute(L2_joint, ar_med = round(parity_in_AR__med_lfc, 2),
            br1_int = round(parity_x_HR_BRCA1__med_lfc, 2),
            br1_abs = round(parity_in_AR__med_lfc + parity_x_HR_BRCA1__med_lfc, 2)) %>%
  as.data.frame() %>% print()

cat("\n\n====================== RISK MAIN EFFECTS ======================\n")
rt <- read_csv("/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/da_pipeline/inquiries/risk_main_effects_20260507/outputs/stageH_findings_table.csv",
                show_col_types = FALSE)
rtL <- rt %>% filter(row_type == "L2") %>% distinct(L2_joint, .keep_all = TRUE)
cat(sprintf("L2 rows: %d\n", nrow(rtL)))

cat("\n--- Top 12 BR1_vs_AR (n_sig >= 30, |med| > 0.10) ---\n")
rtL %>% filter(BR1_vs_AR__n_sig >= 30, abs(BR1_vs_AR__med_lfc) > 0.10) %>%
  arrange(desc(abs(BR1_vs_AR__med_lfc))) %>%
  head(12) %>%
  transmute(L2_joint, n_sig = BR1_vs_AR__n_sig,
            med = round(BR1_vs_AR__med_lfc, 2),
            pct_up = BR1_vs_AR__pct_up) %>% as.data.frame() %>% print()

cat("\n--- Top 12 HRS_vs_AR (n_sig >= 30, |med| > 0.10) ---\n")
rtL %>% filter(HRS_vs_AR__n_sig >= 30, abs(HRS_vs_AR__med_lfc) > 0.10) %>%
  arrange(desc(abs(HRS_vs_AR__med_lfc))) %>%
  head(12) %>%
  transmute(L2_joint, n_sig = HRS_vs_AR__n_sig,
            med = round(HRS_vs_AR__med_lfc, 2),
            pct_up = HRS_vs_AR__pct_up) %>% as.data.frame() %>% print()
