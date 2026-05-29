suppressPackageStartupMessages({ library(dplyr); library(readr) })

ft <- read_csv("/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/da_pipeline/inquiries/parity_x_risk_20260504/outputs/stageH_findings_table.csv",
                show_col_types = FALSE)

cat("====================== LASP family — interaction zoom ======================\n")
for (lasp in c("LASP-major", "LASP-basal", "LASP-KIT")) {
  cat(sprintf("\n--- epi::%s in parity_x_HR_BRCA1 (NG-grain, viable only) ---\n", lasp))
  ng <- ft %>% filter(row_type == "nhoodgroup",
                       contrast == "parity_x_HR_BRCA1",
                       grepl(paste0("epi::", lasp, "_\\d+$"), NhoodGroup_renamed)) %>%
    arrange(desc(n_nhoods_in_group)) %>%
    transmute(NG = NhoodGroup_renamed,
              n_nhoods = n_nhoods_in_group,
              n_sig = n_sig,
              med_lfc = round(group_med_lfc, 2),
              pct_up = round(group_pct_up, 1)) %>%
    head(8)
  print(ng %>% as.data.frame())
}

# Markers per top NG for LASP-major and LASP-basal
marker_dir <- "/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/da_pipeline/inquiries/parity_x_risk_20260504/outputs/stageF3_markers/parity_x_HR_BRCA1"

for (lasp_safe in c("LASP_major", "LASP_basal", "LASP_KIT")) {
  cat(sprintf("\n\n=== epi::%s — top markers per BR1 NG (FDR<0.05, |logFC|>0.5) ===\n",
              gsub("_", "-", lasp_safe)))
  files <- list.files(marker_dir,
                       pattern = sprintf("^epi__%s_\\d+_vs_parent\\.csv$", lasp_safe),
                       full.names = TRUE)
  if (length(files) == 0) {
    cat("  (no marker files)\n"); next
  }
  for (fp in files) {
    ng_id <- gsub(sprintf("^epi__%s_(\\d+)_vs_parent\\.csv$", lasp_safe),
                   "NG_\\1", basename(fp))
    m <- read_csv(fp, show_col_types = FALSE) %>%
      filter(!is_dissoc_stress, !is.na(symbol), symbol != "",
             FDR < 0.05, abs(logFC) > 0.5)
    if (nrow(m) == 0) next
    cat(sprintf("\n%s — %d sig markers\n", ng_id, nrow(m)))
    cat("  UP top8: ");
    cat(m %>% filter(logFC > 0) %>% arrange(desc(logFC)) %>% head(8) %>% pull(symbol) %>%
          paste(collapse=", "))
    cat("\n  DN top8: ");
    cat(m %>% filter(logFC < 0) %>% arrange(logFC) %>% head(8) %>% pull(symbol) %>%
          paste(collapse=", "))
    cat("\n")
  }
}

# Also show: cross-context summary for LASP family at L2 grain
cat("\n\n====================== LASP family — L2-grain cross-context ======================\n")
ftL <- ft %>% filter(row_type == "L2") %>% distinct(L2_joint, .keep_all = TRUE)
lasp_l2 <- ftL %>% filter(L2_joint %in% c("epi::LASP-major","epi::LASP-basal","epi::LASP-KIT")) %>%
  transmute(L2_joint,
            ar_n = parity_in_AR__n_sig,
            ar_med = round(parity_in_AR__med_lfc, 2),
            br1_n = parity_x_HR_BRCA1__n_sig,
            br1_int = round(parity_x_HR_BRCA1__med_lfc, 2),
            br1_abs = round(parity_in_AR__med_lfc + parity_x_HR_BRCA1__med_lfc, 2),
            sporadic_n = parity_x_HR_sporadic__n_sig,
            sporadic_int = round(parity_x_HR_sporadic__med_lfc, 2))
print(lasp_l2 %>% as.data.frame())
