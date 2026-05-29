#Code to add the newly merged iHBCA colData to the old merged iHBCA object

# full merging keeps failing likely due to the size of the object, 
# so try this more simple (but less elegant) approach


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

from scipy.io import mmwrite
import gzip







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
                        'Mastectomy BRCA1': 'HR-BR1',
                        'Mastectomy BRCA2': 'HR-BR2',
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
ihbca_raw_data_dir = '/home/adr44/rds/rds-mammary-TaVjTYER47U/ADR44/data/hbca/iHBCA_raw_data/'

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

#Add dataset name to the obs
adata_gray.obs['dataset'] = 'gray'
adata_murrow.obs['dataset'] = 'murrow'
adata_nee.obs['dataset'] = 'nee'
adata_pal.obs['dataset'] = 'pal'
adata_twigger.obs['dataset'] = 'twigger'
adata_kumar.obs['dataset'] = 'kumar'
adata_reed.obs['dataset'] = 'reed'

#Merge the colData only
colDatas = {'gray': adata_gray.obs,
            'murrow': adata_murrow.obs,
            'nee': adata_nee.obs,
            'pal': adata_pal.obs,
            'twigger': adata_twigger.obs,
            'kumar': adata_kumar.obs,
            'reed': adata_reed.obs}

merged_colData = pd.concat(colDatas.values(), axis=0, keys=colDatas.keys(), names=['dataset'])
#set the index to the cellID but keep as a column
merged_colData = merged_colData.set_index('cellID', drop=False)

#save the merged colData
os.makedirs(ihbca_raw_data_dir + '/merged', exist_ok=True)
merged_colData.to_csv(ihbca_raw_data_dir + '/merged/merged_colData.csv')
merged_colData = pd.read_csv(ihbca_raw_data_dir + '/merged/merged_colData.csv', index_col=0, dtype=str)
#rename first ((non-index) column to 'cellID'
merged_colData.columns = ['cellID'] + list(merged_colData.columns[1:])

#Load the old merged iHBCA object
adata_iHBCA = sc.read_h5ad('/home/adr44/rds/hpc-work/hbca/DATA_FINAL/cellxgene/integrated_HBCA_cellxgene.h5ad')
adata_iHBCA.X = adata_iHBCA.raw.X

#remove adata_iHBCA.raw to save on RAM
adata_iHBCA.raw = None

#check compatibility
print('Old merged iHBCA object:')
print(adata_iHBCA.obs.shape)
print(adata_iHBCA.obs.head())
print('New merged colData:')
print(merged_colData.shape)
print(merged_colData.head())

#The new object has more cells (and metadata) than the old object
new_cells = set(merged_colData.cellID)
old_cells = set(adata_iHBCA.obs.index)
old_cells_not_in_new = old_cells - new_cells
print('Old cells not in new:', len(old_cells_not_in_new))
new_cells_not_in_old = new_cells - old_cells
print('New cells not in old:', len(new_cells_not_in_old))

#explore these cells - they are the misslabelled Pal cells (and those that due to misslabelling were not included in the old object)
df_old_cells_not_in_new = adata_iHBCA.obs.loc[list(old_cells_not_in_new)]
df_new_cells_not_in_old = merged_colData.loc[list(new_cells_not_in_old)]
print(df_old_cells_not_in_new.head())
print(df_new_cells_not_in_old.head())

#Fix cell inclusion in the old object to the new one
old_cells_to_keep = old_cells - old_cells_not_in_new
old_cells_to_keep_str = [str(cell) for cell in old_cells_to_keep] 
adata_pal_new = adata_pal[adata_pal.obs.index.isin(list(new_cells_not_in_old)), :]
adata_pal_new.obs['donor_id'] = adata_pal_new.obs['patientID']

adata_pal_new.layers = None

#remove other adata objects from memory (except pal)
del adata_gray
del adata_murrow
del adata_nee
# del adata_pal
del adata_twigger
del adata_kumar
del adata_reed

#Add the new cells to the old object
adata_iHBCA = adata_iHBCA.concatenate(adata_pal_new, join='inner')

#Reorder to match merged_colData
adata_iHBCA.obs['cellID'] = adata_iHBCA.obs.index
adata_iHBCA.obs['cellID_2'] = adata_iHBCA.obs.cellID.str.split('-').str[:-1].str.join('-')
adata_iHBCA.obs.set_index('cellID_2', inplace=True)
adata_iHBCA.obs.index.name = 'cellID'

#subset to the cells in the new colData
adata_iHBCA = adata_iHBCA[merged_colData.cellID, :]

#Again #check compatibility
print('Updated merged iHBCA object:')
print(adata_iHBCA.obs.shape)
print(adata_iHBCA.obs.head())
print('New merged colData:')
print(merged_colData.shape)
print(merged_colData.head())

#check if the indices are the same
print('Check if the indices are the same:')
print(Counter(adata_iHBCA.obs.index == merged_colData.index))

#add the new colData to the old merged iHBCA object
adata_iHBCA.obs = merged_colData

#Fix the var colnames
adata_iHBCA.var.columns = ['highly_variable-old', 'means-old', 'dispersions-old', 'dispersions_norm-old', 'symbol', 'gene_id']
adata_iHBCA.var = adata_iHBCA.var[['symbol', 'gene_id']]

#Fix age and parous columns
adata_iHBCA.obs['parous'] = adata_iHBCA.obs['parous'].astype(str)
adata_iHBCA.obs['parity'] = adata_iHBCA.obs['parity'].astype(str)
adata_iHBCA.obs['age'] = adata_iHBCA.obs['age'].astype(str)
adata_iHBCA.obs['pred_spikein'] = adata_iHBCA.obs['pred_spikein'].astype(str)

#Save the new merged iHBCA object
adata_iHBCA.write(ihbca_raw_data_dir + '/merged/merged_HBCA_inner_ENSEMBL_simple_prenormalisation.h5ad')
# adata_iHBCA = sc.read_h5ad(ihbca_raw_data_dir + '/merged/merged_HBCA_inner_ENSEMBL_simple_prenormalisation.h5ad')

#Make logcounts and simple preintegration UMAP
adata_iHBCA.layers['raw'] = adata_iHBCA.X.copy()

sc.pp.normalize_per_cell(adata_iHBCA, counts_per_cell_after=1e4)
sc.pp.log1p(adata_iHBCA)

adata_iHBCA.layers['lognorm'] = adata_iHBCA.X.copy()

sc.pp.highly_variable_genes(adata_iHBCA, layer='lognorm', n_top_genes=5000)
sc.pp.pca(adata_iHBCA, n_comps=50, use_highly_variable=True, svd_solver='arpack')
sc.pp.neighbors(adata_iHBCA, use_rep='X_pca', n_neighbors=15)
sc.tl.umap(adata_iHBCA)

#Plot the preintegration UMAP
sc.settings.figdir = ihbca_raw_data_dir + '/merged/plots/'
os.makedirs(sc.settings.figdir, exist_ok=True)
sc.pl.umap(adata_iHBCA, color=['dataset', 'level1'], ncols=1, save='_preintegration.png')




#save the new merged iHBCA object
adata_iHBCA.write(ihbca_raw_data_dir + '/merged/preintegration_HBCA_inner_ENSEMBL_simple.h5ad')
# adata_iHBCA = sc.read_h5ad(ihbca_raw_data_dir + '/merged/preintegration_HBCA_inner_ENSEMBL_simple.h5ad')

#Save just the raw counts matrix (.mtx.gz)
# scipy.sparse.save_npz(ihbca_raw_data_dir + '/merged/preintegration_HBCA_inner_ENSEMBL_simple_counts_matrix.npz', adata_iHBCA.layers['raw'])
# with gzip.open(ihbca_raw_data_dir + '/merged/preintegration_HBCA_inner_ENSEMBL_simple_counts_matrix.mtx.gz', 'wb') as f:
#     mmwrite(f, scipy.sparse.csr_matrix(adata_iHBCA.layers['raw']))
scipy.sparse.save_npz(
    ihbca_raw_data_dir + '/merged/preintegration_HBCA_inner_ENSEMBL_simple_counts_matrix.npz',
    scipy.sparse.csr_matrix(adata_iHBCA.layers['raw']),
    compressed=True
)

# ### FIX THE RISK STATUS LABELS
# # adata_iHBCA.obs['risk_status2'] = adata_iHBCA.obs['risk_status']
# risk_status_dict = {'Mammoplasty WT': 'AR', 
#                     'Mastectomy WT': 'HR-Unk',
#                     'Mastectomy BRCA1': 'HR-BR1',
#                     'Mastectomy BRCA2': 'HR-BR2',
#                     'Mastectomy unknown': 'HR-Unk',
#                     'Contralateral BRCA1': 'HR-cBR1'}
# adata_iHBCA.obs['risk_status'][adata_iHBCA.obs.dataset == 'reed'] = adata_iHBCA.obs['tissue_condition'][adata_iHBCA.obs.dataset == 'reed'].map(risk_status_dict)

print('Saved the new merged iHBCA object.')

#Make a subsetted object and save (100k cells)
np.random.seed(0)
adata_sub = adata_iHBCA[np.random.choice(adata_iHBCA.obs.index, 100000, replace=False), :]
#save
adata_sub.write(ihbca_raw_data_dir + '/merged/preintegration_HBCA_inner_ENSEMBL_simple_sub100k.h5ad')

#Make a subsetted object and save (400k cells)
np.random.seed(0)
adata_sub = adata_iHBCA[np.random.choice(adata_iHBCA.obs.index, 400000, replace=False), :]
#save
adata_sub.write(ihbca_raw_data_dir + '/merged/preintegration_HBCA_inner_ENSEMBL_simple_sub400k.h5ad')


#subset 40k cells per dataset
np.random.seed(0)
cells_to_subset = []
for dataset in adata_iHBCA.obs['dataset'].unique():
    all_cells_dataset = adata_iHBCA.obs[adata_iHBCA.obs['dataset'] == dataset].index
    sub40k_cells_dataset = np.random.choice(all_cells_dataset, 40000, replace=False)
    cells_to_subset.append(sub40k_cells_dataset)
#subset
cells_to_subset = np.concatenate(cells_to_subset)
adata_sub = adata_iHBCA[adata_iHBCA.obs.index.isin(cells_to_subset), :]
#save
adata_sub.write(ihbca_raw_data_dir + '/merged/preintegration_HBCA_inner_ENSEMBL_simple_40kPerDataset.h5ad')






#Save a smaller version with just the lognorm layer
adata_iHBCA.layers = {}
adata_iHBCA.write(ihbca_raw_data_dir + '/merged/preintegration_HBCA_inner_ENSEMBL_simple_lognorm_only.h5ad')
adata_iHBCA.obsm['X_pca'].tofile(ihbca_raw_data_dir + '/merged/merged_iHBCA_pca50.csv', sep=',')




