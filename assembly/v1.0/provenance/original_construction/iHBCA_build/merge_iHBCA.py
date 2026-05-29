#Code to merge and harmonise iHBCA data for integration


# conda env py-scvi

#libraries
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
import anndata as ad

import glob
import os
from collections import Counter
import matplotlib.pyplot as plt






##### Functions #####

#data loading
def get_adata_Gray():
    data_dir = ihbca_raw_data_dir + "/gray_data/formatted/"
    adata = ad.read_mtx(data_dir+'gray_counts.mtx').T
    metadata = pd.read_csv(data_dir+'gray_phenodata.csv', index_col=0)
    genemap = pd.read_csv(data_dir+'gray_features.csv', index_col=0)
    adata.obs = metadata
    adata.var = genemap
    adata.var['gene_id'] = adata.var['gene_id'].astype(str)
    #
    adata.layers['raw'] = adata.X.copy()
    #
    return(adata)


def get_adata_Murrow():
    data_dir = ihbca_raw_data_dir + "/murrow_data/formatted/"
    adata = ad.read_mtx(data_dir+'murrow_counts.mtx').T
    metadata = pd.read_csv(data_dir+'murrow_phenodata.csv', index_col=0)
    genemap = pd.read_csv(data_dir+'murrow_features.csv', index_col=0)
    adata.obs = metadata
    adata.var = genemap
    adata.var['gene_id'] = adata.var['gene_id'].astype(str)
    #
    # adata.obs['batch'] = adata.obs['Batch']
    # adata.obs['patientID'] = adata.obs['Sample']
    #
    adata.layers['raw'] = adata.X.copy()
    
    return(adata)


def get_adata_Nee():
    data_dir = ihbca_raw_data_dir + "/nee_data/formatted/"
    adata = ad.read_mtx(data_dir+'nee_counts.mtx').T
    metadata = pd.read_csv(data_dir+'nee_phenodata.csv', index_col=0, engine='python')
    genemap = pd.read_csv(data_dir+'nee_features.csv', index_col=0)
    adata.obs = metadata
    adata.var = genemap
    adata.var['gene_id'] = adata.var['gene_id'].astype(str)
    #
    adata.layers['raw'] = adata.X.copy()
    #
    return(adata)


def get_adata_Pal():
    data_dir = ihbca_raw_data_dir + "/pal_data/formatted/"
    adata = ad.read_mtx(data_dir+'pal_counts.mtx').T
    metadata = pd.read_csv(data_dir+'pal_phenodata.csv', index_col=0, low_memory=False)
    genemap = pd.read_csv(data_dir+'pal_features.csv', index_col=0)
    adata.obs = metadata
    adata.var = genemap
    adata.var['gene_id'] = adata.var['gene_id'].astype(str)
    #
    adata.layers['raw'] = adata.X.copy()
    #
    return(adata)


def get_adata_Twigger():
    data_dir = ihbca_raw_data_dir + "/twigger_data/formatted/"
    adata = ad.read_mtx(data_dir+'twigger_counts.mtx').T
    metadata = pd.read_csv(data_dir+'twigger_phenodata.csv', index_col=0)
    genemap = pd.read_csv(data_dir+'twigger_features.csv', index_col=0)
    adata.obs = metadata
    adata.var = genemap
    adata.var['gene_id'] = adata.var['gene_id'].astype(str)
    #
    adata.layers['raw'] = adata.X.copy()
    #
    return(adata)


def get_adata_Kumar():
    data_dir = ihbca_raw_data_dir + "/kumar_data/formatted/"
    adata = ad.read_mtx(data_dir+'kumar_counts.mtx').T
    metadata = pd.read_csv(data_dir+'kumar_phenodata.csv', index_col=0)
    genemap = pd.read_csv(data_dir+'kumar_features.csv', index_col=0)
    adata.obs = metadata
    adata.var = genemap
    adata.var['gene_id'] = adata.var['gene_id'].astype(str)
    #
    adata.layers['raw'] = adata.X.copy()
    #
    return(adata)


def get_adata_Reed():
    data_dir = ihbca_raw_data_dir + "/reed_data/formatted/"
    adata = sc.read_h5ad(data_dir + 'HBCA_scVI_processing_date_2023-06-19.h5ad')
    adata.obs['batch'] = adata.obs['processing_date']
    adata.obs['cellID'] = adata.obs.index
    adata.obs['age'] = adata.obs.patient_age
    parous_dict = {'0': False, '1': True, '2': True, '3': True, '4': True, '5': True, 'unknown': None}
    adata.obs['parous'] = adata.obs.parity.map(parous_dict)
    risk_status_dict = {'Mammoplasty WT': 'AR', 
                        'Mastectomy WT': 'HR-Unk',
                        'Mammoplasty BRCA1': 'HR-BR1',
                        'Mammoplasty BRCA2': 'HR-BR2',
                        'Mastectomy unknown': 'HR-Unk',
                        'Contralateral BRCA1': 'HR-cBR1'}
    adata.obs['risk_status'] = adata.obs.tissue_condition.map(risk_status_dict)
    adata.obs['tissue_origin'] = 'frozen'
    FACS_dict = {'Supernatant unsorted': 'not_sorted',
                    'Supernatant live-sorted': 'live_sorted',
                    'Organoid unsorted': 'not_sorted',
                    'Organoid LP sorted': 'cell_type_sorted'}
    adata.obs['FACS_status'] = adata.obs.before.map(FACS_dict)
    adata.obs['sample_type'] = adata.obs.sample_type_coarse
    adata.var['gene_id'] = adata.var['gene_ids']
    adata.X = adata.layers['raw'].copy()
    #
    return(adata)


#filtering
def filter_ENSG(adata):
    # filter_geneID = [isinstance(geneID, str) for geneID in adata.var['gene_id'].values]
    adata.var['gene_id'] = adata.var['gene_id'].astype(str)
    filter_geneID = [geneID != 'nan' for geneID in adata.var['gene_id'].values]
    adata = adata[:,filter_geneID].copy()
    adata.var_names = adata.var['gene_id'].values
    adata.var_names_make_unique()
    return(adata)





##### Setup iHBCA data #####
ihbca_raw_data_dir = '/home/adr44/rds/hpc-work/hbca/data/integration/iHBCA_raw_data/'

#load data
adata_gray = get_adata_Gray()
adata_murrow = get_adata_Murrow()
adata_nee = get_adata_Nee()
adata_pal = get_adata_Pal()
adata_twigger = get_adata_Twigger()
adata_kumar = get_adata_Kumar()
adata_reed = get_adata_Reed()

#dataset merging and harmonisation]

adata_gray = filter_ENSG(adata_gray)
adata_murrow = filter_ENSG(adata_murrow)
adata_nee = filter_ENSG(adata_nee)
adata_pal = filter_ENSG(adata_pal)
adata_twigger = filter_ENSG(adata_twigger)
adata_kumar = filter_ENSG(adata_kumar)
adata_reed = filter_ENSG(adata_reed)




# Merge adata
adatas = {
    "Reed": adata_reed,
    "Kumar": adata_kumar,
    "Nee": adata_nee,
    "Twigger": adata_twigger,
    "Murrow": adata_murrow,
    "Pal": adata_pal,
    "Gray": adata_gray}

adata = ad.concat(adatas, label="dataset", join="inner")

# #remove other adata objects from memory
del adata_gray
del adata_murrow
del adata_nee
del adata_pal
del adata_twigger
del adata_kumar
del adata_reed

print('Shape of merged adata:', adata)
print('Shape of merged counts:', adata.X.shape)

#setup
adata.layers['raw'] = adata.X.copy()

sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)

adata.layers['lognorm'] = adata.X.copy()

#do a preintegration hvg/pca/umap
sc.pp.highly_variable_genes(adata, layer='lognorm', n_top_genes=5000)
sc.pp.pca(adata, n_comps=50, use_highly_variable=True, svd_solver='arpack')
sc.pp.neighbors(adata, use_rep='X_pca', n_neighbors=15)
sc.tl.umap(adata)

#fix issue with parous column
adata.obs['parous'] = adata.obs['parous'].astype(str)
adata.obs['parity'] = adata.obs['parity'].astype(str)
adata.obs['age'] = adata.obs['age'].astype(str)

print('Done merging and harmonisation of iHBCA data.')

#save preintegration data
adata.write_h5ad(ihbca_raw_data_dir + '/preintegration_HBCA_inner_ENSEMBL.h5ad')
# adata.read_h5ad(ihbca_raw_data_dir + '/preintegration_HBCA_inner_ENSEMBL.h5ad')

print('Done saving pre-integration iHBCA data.')
