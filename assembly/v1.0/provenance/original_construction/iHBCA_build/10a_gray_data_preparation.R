#R

#Prepare the Gray et al dataset for joining and integration in iHBCA.

#conda env milo



#libraries
library(scran)
library(AnnotationDbi)
library(org.Hs.eg.db)

library(plyr)
library(dplyr)
library(stringr)
library(Matrix)

#read dir
rawdata_dir <- '/home/adr44/rds/hpc-work/hbca/data/integration/iHBCA_raw_data/'

#load data
gray_counts <- read.table(paste0(rawdata_dir, '/gray_data/original/GSE180878_Li_Brugge_10XscRNAseq_GeneCellMatrix_RNAcounts_human.csv'),
                          sep=',', header = TRUE, row.names = 1)
gray_meta <- read.table(paste0(rawdata_dir, '/gray_data/original/GSE180878_Li_Brugge_10XscRNAseq_Metadata_human.csv'), 
                        sep=',', header=TRUE)
gray_patient_meta <- read.csv(paste0(rawdata_dir, '/gray_data/original/patient_meta.csv'))

#make seperate counts, rowdata and coldata of set formatting
dir.create(paste0(rawdata_dir, '/gray_data/formatted'))


#rowdata
gray_rowData <- data.frame(row.names = rownames(gray_counts),
                           'symbol' = rownames(gray_counts))
gray_rowData$gene_id <- mapIds(org.Hs.eg.db,
                                  keys=gray_rowData$symbol,
                                  column="ENSEMBL",
                                  keytype="SYMBOL",
                                  multiVals="first")
#rownames(gray_rowdata) <- gray_rowdata$gene_id #fails as some ensemble IDs are not found
print(head(gray_rowData))
write.csv(gray_rowData, paste0(rawdata_dir, '/gray_data/formatted/gray_features.csv'))

#counts
colnames(gray_counts) <- sapply(strsplit(colnames(gray_counts), '\\.'), function(x) paste0(x[1], '-', x[2]))
gray_counts <- as.matrix(gray_counts)
gray_counts <- as(gray_counts, "sparseMatrix")
writeMM(gray_counts, file = paste0(rawdata_dir, '/gray_data/formatted/gray_counts.mtx'))


#coldata

#Grab and format the coldata info
rownames(gray_meta) <- gray_meta$cellID
gray_meta <- gray_meta[colnames(gray_counts), ]

#add extra metadata
patients <- sapply(strsplit(gray_meta$cellID, '_'), function(x) x[1])
gray_patient_meta$tissue_condition <- paste0(gray_patient_meta$Surgery, ' ', word(gray_patient_meta$Genotype..mutation., 1))
gray_patient_meta$risk_status <- mapvalues(gray_patient_meta$tissue_condition, 
                                           from=c('Contralateral prophylactic mastectomy BRCA1',
                                                  'Contralateral prophylactic mastectomy BRCA2',
                                                  'Prophylactic mastectomy BRCA1',
                                                  'Prophylactic mastectomy BRCA2',
                                                  'Reductive mammoplasty RAD51C',
                                                  'Reductive mammoplasty WT'),
                                           to=c('HR-cBR1', 'HR-cBR2', 'HR-BR1', 'HR-BR2', 'HR-RAD', 'AR'))
gray_patient_meta_sub <- gray_patient_meta[,c('scRNA.seq', 'Age', 'Births', 'risk_status')]
names(gray_patient_meta_sub) <- c('patientID', 'patient_age', 'parity', 'risk_status')
gray_meta$patientID <- patients
gray_meta <- merge(gray_meta, unique(gray_patient_meta_sub), by='patientID', all.x=T, sort=F)
gray_meta$FACS_status <- 'not_sorted'
gray_meta$sample_type <- 'mixed'
gray_meta$tissue_origin <- 'fresh'

gray_colData <- data.frame(row.names = gray_meta$cellID,
                           'cellID' = gray_meta$cellID,
                           'patientID' = gray_meta$patientID,
                           'sampleID' = gray_meta$patientID,
                           'age' = gray_meta$patient_age,
                           'parous' = gray_meta$parity > 0,
                           'parity' = gray_meta$parity,
                           'risk_status' = gray_meta$risk_status,
                           'tissue_origin' = gray_meta$tissue_origin,
                           'FACS_status' = gray_meta$FACS_status,
                           'sample_type' = gray_meta$sample_type,
                           'batch' = gray_meta$patientID,
                           'level0' = mapvalues(gray_meta$Major.subtype,
                                                from = c("AV", "BA", "HS", "Fibroblast", "VasLymph", "Immune"),
                                                to = c("Epithelial", "Epithelial", "Epithelial", "Stroma", "Stroma", "Immune")),
                           'level1' = gray_meta$Major.subtype,
                           'level2' = gray_meta$Cell.subtype)

print(head(gray_colData))
write.csv(gray_colData, paste0(rawdata_dir, '/gray_data/formatted/gray_phenodata.csv'))
