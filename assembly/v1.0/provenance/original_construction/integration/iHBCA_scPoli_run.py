#Code to run scPoli integration for the iHBCA data


# conda env py-scvi

#libraries
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
import anndata as ad
from sklearn.metrics import classification_report
from sklearn.decomposition import KernelPCA
import mygene
mg = mygene.MyGeneInfo()

import glob
import os
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns
import sys

from scarches.models.scpoli import scPoli
import torch
device = torch.device("cuda")



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
    # print("Diffusion map")
    # sc.tl.diffmap(adata)
    #
    return adata





##### Analysis #####

#load data - iHBCA 
ihbca_raw_data_dir = '/home/adr44/rds/rds-mammary-TaVjTYER47U/ADR44/data/hbca/iHBCA_raw_data/'
adata = sc.read_h5ad(ihbca_raw_data_dir + '/merged/preintegration_HBCA_inner_ENSEMBL_simple.h5ad')

# #make random 100000 cell subset for testing things
# np.random.seed(0)
# adata_sub = adata[np.random.choice(adata.obs.index, 100000, replace=False), :]
# #save
# adata_sub.write(ihbca_raw_data_dir + '/merged/preintegration_HBCA_inner_ENSEMBL_simple_sub100k.h5ad')

# Run scPoli

# for latent_dim in [20, 50, 100, 200]:
latent_dim = int(sys.argv[1])

# specify adata object to use
adata_use = adata
# Set hidden layer size
if (latent_dim == 20) | (latent_dim == 50):
    hidden_layer_size = 128
else:
    hidden_layer_size = 512
# print params
print("latent_dim: " + str(latent_dim))
print("hidden_layer_size: " + str(hidden_layer_size)) 
# scPoli setup
early_stopping_kwargs = {
    "early_stopping_metric": "val_prototype_loss",
    "mode": "min",
    "threshold": 0,
    "patience": 20,
    "reduce_lr": True,
    "lr_patience": 13,
    "lr_factor": 0.1,
}
#make sure batch, dataset and risk_status are strings - otherwise you get an error
adata_use.obs['batch'] = adata_use.obs['batch'].astype(str)
adata_use.obs['risk_status'] = adata_use.obs['risk_status'].astype(str)
adata_use.obs['dataset'] = adata_use.obs['dataset'].astype(str)
# make sure no duplicate batch labels across datasets
adata_use.obs['batch'] = adata_use.obs['dataset'] + '_' + adata_use.obs['batch']
scpoli_model = scPoli(
    adata=adata_use,
    condition_keys=['batch', 'risk_status'],
    cell_type_keys=None,
    embedding_dims=[10, 10],
    latent_dim=latent_dim,
    hidden_layer_sizes=[hidden_layer_size],
    recon_loss='nb',
)
# Train scPoli
scpoli_model.train(
    n_epochs=50,
    pretraining_epochs=40,
    early_stopping_kwargs=early_stopping_kwargs,
    eta=5,
)
print("Training finished")
# Save model and adata
save_dir = '/home/adr44/rds/hpc-work/hbca/data/integration/scPoli/output/'
os.makedirs(save_dir + '/model/', exist_ok=True)
model_path = save_dir + '/model/iHBCA_allcells_ld=' + str(latent_dim) + '_model.pt'
#Add time to mode for uniqueness
time = pd.Timestamp.now().strftime("%Y%m%d%H%M%S")
scpoli_model.save(save_dir + '/model/iHBCA_allcells_ld=' + str(latent_dim) + '_model_t' + str(time) +'.pt') 
# scpoli_model = scPoli.load(model_path, adata=adata_use)
#post dimred processing
adata_use.obsm['X_scPoli'] = scpoli_model.get_latent(adata_use)
adata_use = post_dimred_processing(adata_use, dimred='X_scPoli', run_clustering=True)
#save X_scPoli, X_UMAP, and metadata
os.makedirs(save_dir + 'n_dims_' + str(latent_dim), exist_ok=True)
umap_df = pd.DataFrame(adata_use.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
umap_df.to_csv(save_dir + 'n_dims_' + str(latent_dim) + '/UMAP.csv')
scPoli_df = pd.DataFrame(adata_use.obsm['X_scPoli'], columns=['scPoli' + str(i) for i in range(latent_dim)])
scPoli_df.to_csv(save_dir + 'n_dims_' + str(latent_dim) + '/X_scVI.csv')
adata_use.var.to_csv(save_dir + 'n_dims_' + str(latent_dim) + '/gene_data.csv')
adata_use.obs.to_csv(save_dir + 'n_dims_' + str(latent_dim) + '/cell_data.csv')


