suppressPackageStartupMessages({ library(dplyr); library(readr) })

cat("====================== PARITY x RISK ======================\n\n")
ft <- read_csv("/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/da_pipeline/inquiries/parity_x_risk_20260504/outputs/stageH_findings_table.csv",
                show_col_types = FALSE)
cat(sprintf("findings_table: %d rows, %d cols\n", nrow(ft), ncol(ft)))
cat("row_type:\n"); print(ft %>% count(row_type))

cat("\n--- L2s with parity_in_AR signal (n_sig >= 30, |med_lfc| > 0.10) ---\n")
ar_l2 <- ft %>% filter(row_type == "L2",
                        parity_in_AR__n_sig >= 30,
                        abs(parity_in_AR__med_lfc) > 0.10) %>%
  arrange(desc(abs(parity_in_AR__med_lfc))) %>%
  transmute(L2_joint,
            ar_n = parity_in_AR__n_sig,
            ar_med = round(parity_in_AR__med_lfc, 2),
            ar_pctup = parity_in_AR__pct_up,
            br1_n = parity_x_HR_BRCA1__n_sig,
            br1_med = round(parity_x_HR_BRCA1__med_lfc, 2),
            br1_pctup = parity_x_HR_BRCA1__pct_up)
print(ar_l2, n = Inf)

cat("\n--- BR1 modifications (parity_x_HR_BRCA1, n_sig >= 30, |med_lfc| > 0.20) ---\n")
br1_mods <- ft %>% filter(row_type == "L2",
                           parity_x_HR_BRCA1__n_sig >= 30,
                           abs(parity_x_HR_BRCA1__med_lfc) > 0.20) %>%
  arrange(desc(abs(parity_x_HR_BRCA1__med_lfc))) %>%
  transmute(L2_joint,
            ar_med = round(parity_in_AR__med_lfc, 2),
            br1_int = round(parity_x_HR_BRCA1__med_lfc, 2),
            br1_abs = round(parity_in_AR__med_lfc + parity_x_HR_BRCA1__med_lfc, 2),
            br1_n = parity_x_HR_BRCA1__n_sig)
print(br1_mods, n = Inf)

cat("\n--- L2s where BR1 reverses parity direction ---\n")
flips <- ft %>% filter(row_type == "L2",
                        parity_in_AR__n_sig >= 30,
                        parity_x_HR_BRCA1__n_sig >= 5,
                        sign(parity_in_AR__med_lfc) != sign(parity_x_HR_BRCA1__med_lfc),
                        abs(parity_in_AR__med_lfc) > 0.10,
                        abs(parity_x_HR_BRCA1__med_lfc) > 0.10) %>%
  transmute(L2_joint,
            ar_med = round(parity_in_AR__med_lfc, 2),
            br1_int = round(parity_x_HR_BRCA1__med_lfc, 2),
            br1_abs = round(parity_in_AR__med_lfc + parity_x_HR_BRCA1__med_lfc, 2))
print(flips, n = Inf)

cat("\n\n====================== RISK MAIN EFFECTS ======================\n\n")
rt <- read_csv("/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/da_pipeline/inquiries/risk_main_effects_20260507/outputs/stageH_findings_table.csv",
                show_col_types = FALSE)
cat(sprintf("findings_table: %d rows, %d cols\n", nrow(rt), ncol(rt)))
cat("row_type:\n"); print(rt %>% count(row_type))

cat("\n--- BR1 main effect (BR1_vs_AR, n_sig >= 30, |med_lfc| > 0.10) ---\n")
if ("BR1_vs_AR__n_sig" %in% colnames(rt)) {
  br1_main <- rt %>% filter(row_type == "L2",
                             BR1_vs_AR__n_sig >= 30,
                             abs(BR1_vs_AR__med_lfc) > 0.10) %>%
    arrange(desc(abs(BR1_vs_AR__med_lfc))) %>%
    transmute(L2_joint,
              n_sig = BR1_vs_AR__n_sig,
              med = round(BR1_vs_AR__med_lfc, 2),
              pct_up = BR1_vs_AR__pct_up)
  print(br1_main, n = Inf)
}

cat("\n--- HRS main effect (HRS_vs_AR, n_sig >= 30, |med_lfc| > 0.10) ---\n")
if ("HRS_vs_AR__n_sig" %in% colnames(rt)) {
  hrs_main <- rt %>% filter(row_type == "L2",
                             HRS_vs_AR__n_sig >= 30,
                             abs(HRS_vs_AR__med_lfc) > 0.10) %>%
    arrange(desc(abs(HRS_vs_AR__med_lfc))) %>%
    transmute(L2_joint,
              n_sig = HRS_vs_AR__n_sig,
              med = round(HRS_vs_AR__med_lfc, 2),
              pct_up = HRS_vs_AR__pct_up)
  print(hrs_main, n = Inf)
}
