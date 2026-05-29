suppressPackageStartupMessages({ library(dplyr); library(readr) })

cat("====================== HRS — CD4 / Treg check ======================\n")
rt <- read_csv("/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/da_pipeline/inquiries/risk_main_effects_20260507/outputs/stageH_findings_table.csv",
                show_col_types = FALSE)
rtL <- rt %>% filter(row_type == "L2") %>% distinct(L2_joint, .keep_all = TRUE)

cd4_lineage <- c("imm::CD4_Th_like", "imm::Treg", "imm::Th17",
                  "imm::CD8_Trm", "imm::CD8_Tem", "imm::CD8_Resting",
                  "imm::ZNF683_T", "imm::IFNg_T", "imm::SS")
cd4 <- rtL %>%
  filter(L2_joint %in% cd4_lineage) %>%
  transmute(L2_joint,
            BR1_n = BR1_vs_AR__n_sig,
            BR1_med = round(BR1_vs_AR__med_lfc, 2),
            BR1_pctup = BR1_vs_AR__pct_up,
            HRS_n = HRS_vs_AR__n_sig,
            HRS_med = round(HRS_vs_AR__med_lfc, 2),
            HRS_pctup = HRS_vs_AR__pct_up,
            BR2_n = BR2_vs_AR__n_sig,
            BR2_med = round(BR2_vs_AR__med_lfc, 2),
            BR2_pctup = BR2_vs_AR__pct_up) %>%
  arrange(HRS_med)
print(cd4 %>% as.data.frame())

cat("\n\n====================== LHS-major in BR1 — NG-level zoom ======================\n")
ft <- read_csv("/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/da_pipeline/inquiries/parity_x_risk_20260504/outputs/stageH_findings_table.csv",
                show_col_types = FALSE)
lhs <- ft %>% filter(row_type == "nhoodgroup",
                      grepl("LHS-major", NhoodGroup_renamed, fixed = TRUE),
                      contrast == "parity_x_HR_BRCA1") %>%
  transmute(NG = NhoodGroup_renamed,
            n_nhoods = n_nhoods_in_group,
            n_sig = n_sig,
            med_lfc = round(group_med_lfc, 2),
            pct_up = round(group_pct_up, 1)) %>%
  arrange(desc(n_nhoods))
cat("\nLHS-major NhoodGroups in parity_x_HR_BRCA1 (interaction):\n")
print(lhs %>% as.data.frame())

# Pull top up + down markers per LHS-major NG in BR1
cat("\nLHS-major BR1 — top markers per NG (limma-voom F.3, FDR<0.05, |logFC|>0.5):\n")
marker_dir <- "/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/da_pipeline/inquiries/parity_x_risk_20260504/outputs/stageF3_markers/parity_x_HR_BRCA1"
files <- list.files(marker_dir, pattern = "^epi__LHS_major_.*_vs_parent\\.csv$", full.names = TRUE)
for (fp in head(files, 5)) {
  ng <- gsub("^epi__LHS_major_(\\d+)_vs_parent\\.csv$", "NG_\\1", basename(fp))
  m <- read_csv(fp, show_col_types = FALSE) %>%
    filter(!is_dissoc_stress, !is.na(symbol), symbol != "",
           FDR < 0.05, abs(logFC) > 0.5)
  cat(sprintf("\n--- %s (%d sig markers, |lfc|>0.5) ---\n", ng, nrow(m)))
  if (nrow(m) > 0) {
    cat("UP top10:  ");
    print(m %>% filter(logFC > 0) %>% arrange(desc(logFC)) %>% head(10) %>% pull(symbol))
    cat("DN top10:  ");
    print(m %>% filter(logFC < 0) %>% arrange(logFC) %>% head(10) %>% pull(symbol))
  }
}

# GSEA NES top per LHS-major NG
gsea_path <- "/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/da_pipeline/inquiries/parity_x_risk_20260504/outputs/stageI3_gsea/per_nhoodgroup_nes.csv"
if (file.exists(gsea_path)) {
  cat("\n\nLHS-major BR1 — top Hallmark/Reactome pathways (NES, padj<0.05):\n")
  g <- read_csv(gsea_path, show_col_types = FALSE)
  cat(sprintf("GSEA file: %d rows, cols: %s\n", nrow(g), paste(colnames(g), collapse=",")))
  if ("group" %in% colnames(g) && "contrast" %in% colnames(g)) {
    g_lhs <- g %>%
      filter(contrast == "parity_x_HR_BRCA1",
             grepl("LHS-major", group, fixed = TRUE),
             padj < 0.05) %>%
      arrange(desc(abs(NES))) %>% head(20)
    print(g_lhs %>% select(any_of(c("group", "pathway", "NES", "padj", "size"))) %>% as.data.frame())
  }
}
