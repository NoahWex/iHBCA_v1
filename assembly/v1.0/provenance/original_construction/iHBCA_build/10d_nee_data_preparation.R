#R

#Prepare the Nee et al dataset for joining and integration in iHBCA.

#conda env milo




#libraries
library(scran)
library(AnnotationDbi)
library(org.Hs.eg.db)

library(plyr)
library(dplyr)
library(Matrix)


#read dir
rawdata_dir <- '/home/adr44/rds/hpc-work/hbca/data/integration/iHBCA_raw_data/'


#load data
mat_list <- list()
mat_list[['ctrl1']] <- read.table(paste0(rawdata_dir, '/nee_data/original/GSM5320152_scRNA_ctrl1_matrix.txt'), sep=' ', header=TRUE, row.names=1)
mat_list[['ctrl2']] <- read.table(paste0(rawdata_dir, '/nee_data/original/GSM5320153_scRNA_ctrl2_matrix.txt'), sep=' ', header=TRUE, row.names=1)
mat_list[['ctrl3']] <- read.table(paste0(rawdata_dir, '/nee_data/original/GSM5320154_scRNA_ctrl3_matrix.txt'), sep=' ', header=TRUE, row.names=1)
mat_list[['ctrl4']] <- read.table(paste0(rawdata_dir, '/nee_data/original/GSM5320155_scRNA_ctrl4_matrix.txt'), sep=' ', header=TRUE, row.names=1)
mat_list[['ctrl5']] <- read.table(paste0(rawdata_dir, '/nee_data/original/GSM5320156_scRNA_ctrl5_matrix.txt'), sep=' ', header=TRUE, row.names=1)
mat_list[['ctrl6']] <- read.table(paste0(rawdata_dir, '/nee_data/original/GSM5320157_scRNA_ctrl6_matrix.txt'), sep=' ', header=TRUE, row.names=1)
mat_list[['ctrl7']] <- read.table(paste0(rawdata_dir, '/nee_data/original/GSM5320158_scRNA_ctrl7_matrix.txt'), sep=' ', header=TRUE, row.names=1)
mat_list[['ctrl8']] <- read.table(paste0(rawdata_dir, '/nee_data/original/GSM5320159_scRNA_ctrl8_matrix.txt'), sep=' ', header=TRUE, row.names=1)
mat_list[['ctrl9']] <- read.table(paste0(rawdata_dir, '/nee_data/original/GSM5320160_scRNA_ctrl9_matrix.txt'), sep=' ', header=TRUE, row.names=1)
mat_list[['ctrl10']] <- read.table(paste0(rawdata_dir, '/nee_data/original/GSM5320161_scRNA_ctrl10_matrix.txt'), sep=' ', header=TRUE, row.names=1)
mat_list[['ctrl11']] <- read.table(paste0(rawdata_dir, '/nee_data/original/GSM5320162_scRNA_ctrl11_matrix.txt'), sep=' ', header=TRUE, row.names=1)
mat_list[['brca1']] <- read.table(paste0(rawdata_dir, '/nee_data/original/GSM5320163_scRNA_brca1_matrix.txt'), sep=' ', header=TRUE, row.names=1)
mat_list[['brca2']] <- read.table(paste0(rawdata_dir, '/nee_data/original/GSM5320164_scRNA_brca2_matrix.txt'), sep=' ', header=TRUE, row.names=1)
mat_list[['brca3']] <- read.table(paste0(rawdata_dir, '/nee_data/original/GSM5320165_scRNA_brca3_matrix.txt'), sep=' ', header=TRUE, row.names=1)
mat_list[['brca4']] <- read.table(paste0(rawdata_dir, '/nee_data/original/GSM5320166_scRNA_brca4_matrix.txt'), sep=' ', header=TRUE, row.names=1)
mat_list[['brca5']] <- read.table(paste0(rawdata_dir, '/nee_data/original/GSM5320167_scRNA_brca5_matrix.txt'), sep=' ', header=TRUE, row.names=1)
mat_list[['brca6']] <- read.table(paste0(rawdata_dir, '/nee_data/original/GSM5320168_scRNA_brca6_matrix.txt'), sep=' ', header=TRUE, row.names=1)
mat_list[['brca7']] <- read.table(paste0(rawdata_dir, '/nee_data/original/GSM5320169_scRNA_brca7_matrix.txt'), sep=' ', header=TRUE, row.names=1)
mat_list[['brca8']] <- read.table(paste0(rawdata_dir, '/nee_data/original/GSM5320170_scRNA_brca8_matrix.txt'), sep=' ', header=TRUE, row.names=1)
mat_list[['brca9']] <- read.table(paste0(rawdata_dir, '/nee_data/original/GSM5320171_scRNA_brca9_matrix.txt'), sep=' ', header=TRUE, row.names=1)
mat_list[['brca10']] <- read.table(paste0(rawdata_dir, '/nee_data/original/GSM5320172_scRNA_brca10_matrix.txt'), sep=' ', header=TRUE, row.names=1)
mat_list[['brca11']] <- read.table(paste0(rawdata_dir, '/nee_data/original/GSM5320173_scRNA_brca11_matrix.txt'), sep=' ', header=TRUE, row.names=1)



#make seperate counts, rowdata and coldata of set formatting
dir.create(paste0(rawdata_dir, '/nee_data/formatted'))


#counts
nee_counts <- do.call(cbind, mat_list)

#fix barcodes
bcs <- sapply(strsplit(colnames(nee_counts), '\\.'), function(x) {
  if (x[2] == 'UCI_Patients'){
    return(paste0(x[2], ' ', x[3], '-', x[4]))
  } else {
    return(paste0(x[2], '-', x[3]))
  }
})

nee_counts <- as.matrix(nee_counts)
nee_counts <- as(nee_counts, "sparseMatrix")

colnames(nee_counts) <- bcs
writeMM(nee_counts, file = paste0(rawdata_dir, '/nee_data/formatted/nee_counts.mtx'))
# nee_counts <- readMM(file = paste0(rawdata_dir, '/nee_data/formatted/nee_counts.mtx')) 

#rowdata

nee_rowData <- data.frame(row.names = rownames(nee_counts),
                             'symbol' = rownames(nee_counts))
nee_rowData$gene_id <- mapIds(org.Hs.eg.db,
                                 keys=nee_rowData$symbol,
                                 column="ENSEMBL",
                                 keytype="SYMBOL",
                                 multiVals="first")
#rownames(nee_rowdata) <- nee_rowdata$gene_id #fails as some ensemble IDs are not found
print(head(nee_rowData))
write.csv(nee_rowData, paste0(rawdata_dir, '/nee_data/formatted/nee_features.csv'))

#coldata

#Grab and format the coldata info
# nee_colData_short2 <- read.csv(paste0(rawdata_dir, '/nee_data/formatted/nee_phenodata_old.csv'), row.names=1)
nee_colData <- read.table(paste0(rawdata_dir, '/nee_data/original/nee_etal_metadata.txt'), sep=' ', header=TRUE, row.names=1)
nee_colData <- nee_colData[colnames(nee_counts), ]
# nee_colData <- nee_colData[rownames(nee_colData_short2), ]

#add extra metadata
nee_colData$FACS_status <- 'cell_type_sorted' #all cells FACs sorted
nee_colData$sample_type <- 'mixed'
nee_colData$risk_status <- mapvalues(nee_colData$tissue, 
                                        from = c('Contralateral', 'Prophylactic Mastectomy', 'Prophylatctic Mastectomy',
                                                 'Reduction Mammoplasty'),
                                        to = c('HR-cBR1', 'HR-BR1', 'HR-BR1', 'AR'))
nee_colData$tissue_origin <- 'frozen'


#save a raw colData
write.csv(nee_colData, paste0(rawdata_dir, '/nee_data/formatted/nee_full_phenodata.csv'))

#subset and save common colData
nee_colData_short <- data.frame(row.names = rownames(nee_colData),
                                'cellID' = rownames(nee_colData),
                                'patientID' = nee_colData$patient_id,
                                'sampleID' = nee_colData$orig.ident,
                                'age' = nee_colData$age,
                                'parous' = mapvalues(nee_colData$parity, 
                                                     from = c('Nulliarous', 'Nulliparous', 'Parous', 'Unknown'),
                                                     to = c(FALSE, FALSE, TRUE, NA)),
                                'parity' = mapvalues(nee_colData$parity, 
                                                     from = c('Nulliarous', 'Nulliparous', 'Parous', 'Unknown'),
                                                     to = c(0, 0, NA, NA)),
                                'risk_status' = nee_colData$risk_status,
                                'tissue_origin' = nee_colData$tissue_origin,
                                'FACS_status' = nee_colData$FACS_status,
                                'sample_type' = nee_colData$sample_type,
                                
                                'batch' = nee_colData$patient_id,
                                'level0' = mapvalues(nee_colData$cell_type_final,
                                                      from = c("Basal", "Endothelial", "Fibroblasts", "Immune", "Luminal1", "Luminal2", 
                                                              "Lymphatic", "Pericytes"),
                                                      to = c("Epithelial", "Stroma", "Stroma", "Immune", "Epithelial", "Epithelial", 
                                                            "Stroma", "Stroma")),
                                'level1' = nee_colData$cell_type_final,
                                'level2' = nee_colData$cell_state)

print(head(nee_colData_short))
write.csv(nee_colData_short, paste0(rawdata_dir, '/nee_data/formatted/nee_phenodata.csv'))
# nee_colData_short2 <- read.csv(paste0(rawdata_dir, '/nee_data/formatted/nee_phenodata_old.csv'), row.names=1)
# write.csv(nee_colData_short2, paste0(rawdata_dir, '/nee_data/formatted/nee_phenodata_old.csv'))

