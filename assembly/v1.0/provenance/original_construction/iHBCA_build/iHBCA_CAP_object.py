# Build the CAP object for iHBCA



#libraries
import scanpy as sc
import pandas as pd
import numpy as np
import scipy
import os



#load data
ihbca_raw_data_dir = '/home/adr44/rds/rds-mammary-TaVjTYER47U/ADR44/data/hbca/iHBCA_raw_data/'
cellxgene_dir = '/home/adr44/rds/rds-mammary-TaVjTYER47U/ADR44/data/hbca/DATA_FINAL/cellxgene/'

# adata = sc.read_h5ad(cellxgene_dir + 'integrated_HBCA_cellxgene.h5ad') #OLD use the updated object now
adata = sc.read_h5ad(ihbca_raw_data_dir + '/merged/preintegration_HBCA_inner_ENSEMBL_simple_lognorm_only.h5ad')
counts_matrix = scipy.sparse.load_npz(ihbca_raw_data_dir + '/merged/preintegration_HBCA_inner_ENSEMBL_simple_counts_matrix.npz')
adata.X = counts_matrix
CAP_obs = pd.read_csv(cellxgene_dir + 'integrated_HBCA_cellxgene_CAP_metadata_NEW.csv', 
                      dtype=str, #errors keep occuring otherwise
                      index_col=0)
adata.obs = CAP_obs.loc[adata.obs_names]

#add cell type annotations
ihbca_level1_5_obs = pd.read_csv('/home/adr44/rds/hpc-work/hbca/data/integration/ihbca_annotation/output/level1.5/adata/ihbca_level1.5_annotations.csv', index_col=0)
adata.obs['celltype_level1'] = ihbca_level1_5_obs.loc[adata.obs_names, 'level1_annotation']
adata.obs['celltype_level2'] = ihbca_level1_5_obs.loc[adata.obs_names, 'level1.5_annotation']

#Add the scVI_100 UMAP and X_scVI embeddings
scvi_embed = pd.read_csv('/home/adr44/rds/hpc-work/hbca/data/integration/scvi/output/iHBCA/global/n_dims_100/X_scVI.csv', index_col=0)
scvi_UMAP = pd.read_csv('/home/adr44/rds/hpc-work/hbca/data/integration/scvi/output/iHBCA/global/n_dims_100/UMAP.csv', index_col=0)
adata.obsm['X_scVI'] = scvi_embed.values
adata.obsm['X_umap'] = scvi_UMAP.values

#Var
adata.var = adata.var[['symbol', 'gene_id']]

#save temp CAP object
adata.write_h5ad(cellxgene_dir + 'iHBCA_CAP_object.h5ad')
# adata = sc.read_h5ad(cellxgene_dir + 'iHBCA_CAP_object.h5ad') #test loading


