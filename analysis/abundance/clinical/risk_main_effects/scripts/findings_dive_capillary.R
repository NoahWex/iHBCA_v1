suppressPackageStartupMessages({ library(dplyr); library(readr); library(tidyr) })

ft <- read_csv("/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/da_pipeline/inquiries/parity_x_risk_20260504/outputs/stageH_findings_table.csv",
                show_col_types = FALSE)

cat("====================== Vas-capillary cross-context ======================\n")
cap_l2 <- ft %>% filter(row_type == "L2") %>%
  filter(L2_joint == "str::Vas-capillary") %>%
  distinct(L2_joint, .keep_all = TRUE)
print(cap_l2 %>% transmute(L2_joint,
                            ar_n = parity_in_AR__n_sig,
                            ar_med = round(parity_in_AR__med_lfc, 2),
                            br1_n = parity_x_HR_BRCA1__n_sig,
                            br1_int = round(parity_x_HR_BRCA1__med_lfc, 2),
                            br1_abs = round(parity_in_AR__med_lfc + parity_x_HR_BRCA1__med_lfc, 2),
                            br2_int = round(parity_x_HR_BRCA2__med_lfc, 2),
                            sp_int = round(parity_x_HR_sporadic__med_lfc, 2)) %>%
        as.data.frame())

cat("\n--- str::Vas-capillary in parity_x_HR_BRCA1 (NG-grain, viable only) ---\n")
ng <- ft %>% filter(row_type == "nhoodgroup",
                     contrast == "parity_x_HR_BRCA1",
                     grepl("Vas-capillary_\\d+$", NhoodGroup_renamed)) %>%
  arrange(desc(n_nhoods_in_group)) %>%
  transmute(NG = NhoodGroup_renamed,
            n_nhoods = n_nhoods_in_group,
            n_sig = n_sig,
            med_lfc = round(group_med_lfc, 2),
            pct_up = round(group_pct_up, 1)) %>%
  head(10)
print(ng %>% as.data.frame())

# Top markers per Vas-capillary NG
marker_dir <- "/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/da_pipeline/inquiries/parity_x_risk_20260504/outputs/stageF3_markers/parity_x_HR_BRCA1"
files <- list.files(marker_dir,
                     pattern = "^str__Vas_capillary_\\d+_vs_parent\\.csv$",
                     full.names = TRUE)
cat(sprintf("\nVas-capillary BR1 marker files: %d\n", length(files)))
for (fp in files) {
  ng_id <- gsub("^str__Vas_capillary_(\\d+)_vs_parent\\.csv$",
                 "NG_\\1", basename(fp))
  m <- read_csv(fp, show_col_types = FALSE) %>%
    filter(!is_dissoc_stress, !is.na(symbol), symbol != "",
           FDR < 0.05, abs(logFC) > 0.5)
  cat(sprintf("\n%s — %d sig markers\n", ng_id, nrow(m)))
  if (nrow(m) == 0) next
  cat("  UP top10: ");
  cat(m %>% filter(logFC > 0) %>% arrange(desc(logFC)) %>% head(10) %>% pull(symbol) %>%
        paste(collapse=", "))
  cat("\n  DN top10: ");
  cat(m %>% filter(logFC < 0) %>% arrange(logFC) %>% head(10) %>% pull(symbol) %>%
        paste(collapse=", "))
  cat("\n")
}

# PTN axis search: PTN ligand + receptors (PTPRZ1, SDC1, SDC2, SDC3, SDC4, NCL, ALK, GPC3)
cat("\n\n====================== PTN signaling axis search ======================\n")
ptn_axis <- c("PTN", "MK", "MDK", "PTPRZ1", "SDC1", "SDC2", "SDC3", "SDC4",
               "NCL", "ALK", "GPC1", "GPC2", "GPC3")
hits <- list()
for (fp in files) {
  ng_id <- gsub("^str__Vas_capillary_(\\d+)_vs_parent\\.csv$",
                 "NG_\\1", basename(fp))
  m <- read_csv(fp, show_col_types = FALSE) %>%
    filter(symbol %in% ptn_axis)
  if (nrow(m) > 0) {
    m$ng <- ng_id
    hits[[ng_id]] <- m
  }
}
ptn_df <- bind_rows(hits)
if (nrow(ptn_df) > 0) {
  cat("\nPTN-axis genes with sig signal in Vas-capillary BR1 NGs (any FDR):\n")
  print(ptn_df %>%
          transmute(ng, symbol, logFC = round(logFC, 2),
                    FDR = signif(FDR, 2),
                    AveExpr = round(AveExpr, 2)) %>%
          arrange(symbol, ng) %>% as.data.frame())
} else {
  cat("\nNo PTN-axis hits at any FDR in Vas-capillary BR1 NG markers.\n")
}

# Also check across LASP-major NGs (PTN appeared there earlier)
cat("\n--- PTN axis in LASP-major BR1 NGs (cross-reference) ---\n")
lasp_files <- list.files(marker_dir,
                          pattern = "^epi__LASP_major_\\d+_vs_parent\\.csv$",
                          full.names = TRUE)
hits2 <- list()
for (fp in lasp_files) {
  ng_id <- gsub("^epi__LASP_major_(\\d+)_vs_parent\\.csv$",
                 "NG_\\1", basename(fp))
  m <- read_csv(fp, show_col_types = FALSE) %>%
    filter(symbol %in% ptn_axis)
  if (nrow(m) > 0) {
    m$ng <- ng_id; hits2[[ng_id]] <- m
  }
}
ptn_lasp <- bind_rows(hits2)
if (nrow(ptn_lasp) > 0) {
  print(ptn_lasp %>%
          transmute(ng, symbol, logFC = round(logFC, 2), FDR = signif(FDR, 2)) %>%
          arrange(symbol, ng) %>% as.data.frame())
}

# And LHS-major NGs
cat("\n--- PTN axis in LHS-major BR1 NGs (cross-reference) ---\n")
lhs_files <- list.files(marker_dir,
                         pattern = "^epi__LHS_major_\\d+_vs_parent\\.csv$",
                         full.names = TRUE)
hits3 <- list()
for (fp in lhs_files) {
  ng_id <- gsub("^epi__LHS_major_(\\d+)_vs_parent\\.csv$",
                 "NG_\\1", basename(fp))
  m <- read_csv(fp, show_col_types = FALSE) %>%
    filter(symbol %in% ptn_axis)
  if (nrow(m) > 0) {
    m$ng <- ng_id; hits3[[ng_id]] <- m
  }
}
ptn_lhs <- bind_rows(hits3)
if (nrow(ptn_lhs) > 0) {
  print(ptn_lhs %>%
          transmute(ng, symbol, logFC = round(logFC, 2), FDR = signif(FDR, 2)) %>%
          arrange(symbol, ng) %>% as.data.frame())
}
