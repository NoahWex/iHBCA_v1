#R

#Prepare the Murrow et al dataset for joining and integration in iHBCA.

#conda env milo


#libraries
library(scran)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(Seurat)

library(plyr)
library(dplyr)
library(Matrix)




#read dir
rawdata_dir <- '/home/adr44/rds/hpc-work/hbca/data/integration/iHBCA_raw_data/'


#load data
srt <- readRDS(paste0(rawdata_dir, '/murrow_data/original/GSE198732_breast.data.rds'))
srt <- UpdateSeuratObject(srt)

#make seperate counts, rowdata and coldata of set formatting
dir.create(paste0(rawdata_dir, '/murrow_data/formatted'))

#rowdata
murrow_rowData <- data.frame(row.names = rownames(srt@assays$RNA@counts),
                           'symbol' = rownames(srt@assays$RNA@counts))
murrow_rowData$gene_id <- mapIds(org.Hs.eg.db,
                               keys=murrow_rowData$symbol,
                               column="ENSEMBL",
                               keytype="SYMBOL",
                               multiVals="first")
#rownames(murrow_rowdata) <- murrow_rowdata$gene_id #fails as some ensemble IDs are not found
print(head(murrow_rowData))
write.csv(murrow_rowData, paste0(rawdata_dir, '/murrow_data/formatted/murrow_features.csv'))

#coldata

#Grab and format the coldata info
murrow_colData <- srt@meta.data

#add extra metadata
murrow_colData$sampleID <- paste(murrow_colData$Sample, murrow_colData$orig.ident, sep = '_')
murrow_colData$FACS_status <- mapvalues(murrow_colData$Sort, 
                                                   from=c('Basal', 'Epithelial', 'Live_singlet', 'Luminal'), 
                                                   to = c('cell_type_sorted', 'cell_type_sorted', 'live_sorted', 'cell_type_sorted'))
murrow_colData$sample_type <- 'mixed'
murrow_colData$risk_status <- 'AR'
murrow_colData$tissue_origin <- 'frozen'

#save a raw colData
write.csv(murrow_colData, paste0(rawdata_dir, '/murrow_data/formatted/murrow_full_phenodata.csv'))

#subset and save common colData
murrow_colData_short <-   data.frame(row.names = rownames(murrow_colData),
                                     'cellID' = rownames(murrow_colData),
                                     'patientID' = murrow_colData$Sample,
                                     'sampleID' = murrow_colData$sampleID,
                                     'age' = murrow_colData$Age,
                                     'parous' = murrow_colData$Parity > 0,
                                     'parity' = murrow_colData$Parity,
                                     'risk_status' = murrow_colData$risk_status,
                                     'tissue_origin' = murrow_colData$tissue_origin,
                                     'FACS_status' = murrow_colData$FACS_status,
                                     'sample_type' = murrow_colData$sample_type,
                                     'batch' = murrow_colData$Batch,
                                     'level0' = mapvalues(murrow_colData$Type,
                                                          from = c("HRpos_Luminal", "Basal", "Secretory_Luminal", 
                                                                   "Fibroblast",  "Vascular_Endothelial", "Vascular_Accessory",
                                                                   "Lymphocyte", "Lymphatic_Endothelial", "Macrophage"),
                                                          to = c("Epithelial", "Epithelial", "Epithelial", "Stroma", "Stroma", "Stroma", "Immune", "Stroma", "Immune")),
                                     'level1' = murrow_colData$Type,
                                     'level2' = NA)

print(head(murrow_colData_short))
write.csv(murrow_colData_short, paste0(rawdata_dir, '/murrow_data/formatted/murrow_phenodata.csv'))

#counts
murrow_counts <- srt@assays$RNA@counts
writeMM(murrow_counts, file = paste0(rawdata_dir, '/murrow_data/formatted/murrow_counts.mtx'))
