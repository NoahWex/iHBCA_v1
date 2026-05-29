#Code to redo scVI integration of HBCA data


# conda env py-scanpy

#libraries
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
import anndata as ad
# import harmonypy as hm

import glob
import os
from collections import Counter
import matplotlib.pyplot as plt




##### Functions #####

def post_dimred_processing(adata, dimred='X_pca', run_clustering=True):
    # neighbourhood graph
    print("knn graph")
    sc.pp.neighbors(adata, use_rep=dimred, n_neighbors=15)
    # UMAP
    print("UMAP")
    sc.tl.umap(adata)
    #
    # Leiden clustering
    if run_clustering:
        list_leiden_res = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        print("Leiden clustering")
        for leiden_res in list_leiden_res:
            print(leiden_res)
            sc.tl.leiden(adata, resolution=leiden_res, key_added='leiden_'+str(leiden_res).replace(".", "_"))
    #
    print("Diffusion map")
    sc.tl.diffmap(adata)
    #
    return adata








##### Analysis #####

#load data
ihbca_raw_data_dir = '/home/adr44/rds/hpc-work/hbca/data/integration/iHBCA_raw_data/'
adata = sc.read_h5ad(ihbca_raw_data_dir + '/merged/preintegration_HBCA_inner_ENSEMBL_simple.h5ad')


#Harmony integration
batchID = 'batch'
sc.external.pp.harmony_integrate(adata, key=batchID, basis='X_pca')

#post processing
adata = post_dimred_processing(adata, dimred='X_pca_harmony', run_clustering=True)

#plot integration results
save_dir = '/home/adr44/rds/hpc-work/hbca/data/integration/harmony/output/'
os.makedirs(save_dir, exist_ok=True)
sc.settings.figdir = save_dir
sc.pl.umap(adata, color=['dataset', 'leiden_0_5', 'leiden_1_0', 'level1'], ncols=1, save='5000hvg_dataset_leiden.png')
sc.pl.umap(adata, color=['level1'], save='_5000hvg_level2.png', legend_loc='on data')
#save X_scVI, X_UMAP, and metadata
umap_df = pd.DataFrame(adata.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
umap_df.to_csv(save_dir + '/5000hvg_UMAP.csv')
pca_hm_df = pd.DataFrame(adata.obsm['X_pca_harmony'], columns=['PC' + str(i) for i in range(1, adata.obsm['X_pca_harmony'].shape[1] + 1)])
pca_hm_df.to_csv(save_dir + '/5000hvg_X_pca_harmony.csv')
adata.var.to_csv(save_dir + '/5000hvg_gene_data.csv')
adata.obs.to_csv(save_dir + '/5000hvg_cell_data.csv')



