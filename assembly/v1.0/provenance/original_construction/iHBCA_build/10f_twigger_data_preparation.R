#R

#Prepare the Twigger et al dataset for joining and integration in iHBCA.

#conda env milo




#libraries
library(scran)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(Seurat)

library(plyr)
library(dplyr)
library(Matrix)



##Read in data
rawdata_dir <- '/home/adr44/rds/hpc-work/hbca/data/integration/iHBCA_raw_data/'

#load data
sce <- readRDS(paste0(rawdata_dir, '/twigger_data/original/sce_all_nospike_2.rds'))


#make seperate counts, rowdata and coldata of set formatting
dir.create(paste0(rawdata_dir, '/twigger_data/formatted'))

#rowdata
twigger_rowData <- data.frame(#row.names = rownames(sce),
                              'symbol' = rownames(sce))
twigger_rowData$gene_id <- mapIds(org.Hs.eg.db,
                              keys=twigger_rowData$symbol,
                              column="ENSEMBL",
                              keytype="SYMBOL",
                              multiVals="first")

#identify the duplicate gene symbols
n_occur <- data.frame(table(twigger_rowData$symbol))
# twigger_rowData[twigger_rowData$symbol %in% n_occur$Var1[n_occur$Freq > 1],]
dupped_symbols <- unique(twigger_rowData[twigger_rowData$symbol %in% n_occur$Var1[n_occur$Freq > 1], 'symbol'])
write.csv(dupped_symbols, paste0(rawdata_dir, '/twigger_data/duplicated_symbols.csv'))

temp_counts <- counts(sce)
#temp_sub <- temp_counts[rownames(temp_counts) %in% dupped_symbols,]
temp_sums <- data.frame('symbols' = names(rowSums(temp_counts)),
                            'sums' = rowSums(temp_counts),
                            'keep' = TRUE)
for (symbol in dupped_symbols){
  temp <- temp_sums[temp_sums$symbols %in% symbol,]
  max_val <- max(temp$sums)
  temp$keep <- temp$sums == max_val
  temp_sums[temp_sums$symbols %in% symbol,] <- temp
}

#use the keep column to remove duplicated genes.
twigger_rowData <- twigger_rowData[temp_sums$keep,]
rownames(twigger_rowData) <- twigger_rowData$symbol

write.csv(twigger_rowData, paste0(rawdata_dir, '/twigger_data/formatted/twigger_features.csv'))


#coldata
twigger_meta <- colData(sce)
twigger_patient_meta <- read.csv(paste0(rawdata_dir, '/twigger_data/original/patient_meta.csv'))
twigger_patient_meta$sampleID <- twigger_patient_meta$patientID
twigger_patient_meta$patientID <- NULL


#add/fix extra metadata
twigger_meta$sampleID <- mapvalues(twigger_meta$Sample, 
                                    from=c('HMC1', 'HMC2', 'HMC2B', 'HMC3', 
                                           'HMC4', 'HMC5', 'HMC6', 'HMC7', 
                                           'HMC8',  'HMC9', 
                                           'RB1', 'RB2', 'RB3', 'RB4', 'RB5', 
                                           'RB6', 'RB7', 'RB8'),
                                    to=c('LMC1', 'LMC2', 'LMC2B', 'LMC3', 
                                         'LMC4', 'LMC5', 'LMC6', 'LMC7', 
                                         'LMC8',  'LMC9', 
                                         'NMC1', 'NMC2', 'NMC3', 'NMC4', 'NMC5', 
                                         'NMC6', 'NMC7', 'NMC1B'))
twigger_meta$patientID <- mapvalues(twigger_meta$Sample, 
                                    from=c('HMC1', 'HMC2', 'HMC2B', 'HMC3', 
                                           'HMC4', 'HMC5', 'HMC6', 'HMC7', 
                                           'HMC8',  'HMC9', 
                                           'RB1', 'RB2', 'RB3', 'RB4', 'RB5', 
                                           'RB6', 'RB7', 'RB8'),
                                    to=c('LMC1', 'LMC2', 'LMC2', 'LMC3', 
                                         'LMC4', 'LMC5', 'LMC6', 'LMC7', 
                                         'LMC8',  'LMC9', 
                                         'NMC1', 'NMC2', 'NMC3', 'NMC4', 'NMC5', 
                                         'NMC6', 'NMC7', 'NMC1'))
twigger_meta <- merge(twigger_meta, unique(twigger_patient_meta), by='sampleID', all.x=T, sort=F)
twigger_meta$FACS_status <- 'no_sort'
twigger_meta$sample_type <- 'mixed'
twigger_meta$risk_status <- 'AR'
twigger_meta$tissue_origin <- 'frozen'



twigger_colData <- data.frame(row.names = twigger_meta$Barcode,
                              'cellID' = twigger_meta$Barcode,
                              'patientID' = twigger_meta$Sample,
                              'sampleID' = twigger_meta$sampleID,
                              'age' = twigger_meta$patient_age,
                              'parous' = twigger_meta$parity > 0,
                              'parity' = twigger_meta$parity,
                              'risk_status' = twigger_meta$risk_status,
                              'tissue_origin' = twigger_meta$tissue_origin,
                              'FACS_status' = twigger_meta$FACS_status,
                              'sample_type' = twigger_meta$sample_type,
                              'batch' = twigger_meta$Batches,
                              'level0' = mapvalues(twigger_meta$Identity, 
                                                   from = c("FB", "LP", "LC1", "LC2", "IM", 
                                                            "EN", "BA", "VA", "HR"),
                                                   to = c("Stroma", "Epithelial", "Epithelial", "Epithelial", "Immune",
                                                          "Stroma", "Epithelial", "Stroma", "Epithelial")),
                              'level1' = twigger_meta$Identity,
                              'level2' = NA)

print(head(twigger_colData))
write.csv(twigger_meta, paste0(rawdata_dir, '/twigger_data/formatted/twigger_full_phenodata.csv'))
write.csv(twigger_colData, paste0(rawdata_dir, '/twigger_data/formatted/twigger_phenodata.csv'))

#counts
twigger_counts <- counts(sce)
twigger_counts <- twigger_counts[temp_sums$keep, ] #remove dupped genes
writeMM(twigger_counts, file = paste0(rawdata_dir, '/twigger_data/formatted/twigger_counts.mtx'))
