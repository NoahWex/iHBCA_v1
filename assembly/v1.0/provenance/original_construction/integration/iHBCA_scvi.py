#Code to redo scVI integration of HBCA data


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

import scvi
import torch
import sys
device = torch.device("cuda")


#set matrix precision for speed up
torch.set_float32_matmul_precision('medium')




##### Functions #####

#integration etc.

def run_scVI(adata, batchID='batch', n_dims=20, n_layers=2):
    # subset to HVGs
    train_adata = adata[:, adata.var['highly_variable']].copy()
    scvi.model.SCVI.setup_anndata(train_adata, layer='raw', batch_key=batchID)
    #
    arches_params = {
        "use_layer_norm": "both",
        "use_batch_norm": "none",
        "encode_covariates": True,
        "dropout_rate": 0.2,
        "n_layers": n_layers,
    }
    #
    vae = scvi.model.SCVI(train_adata, n_latent=n_dims, gene_likelihood="nb", **arches_params)
    # vae.train() # 
    vae.train(early_stopping=True,
        train_size=0.9,
        early_stopping_patience=45,
        max_epochs=400, 
        batch_size=1024, 
        limit_train_batches=20)
    #
    adata.obsm['X_scVI'] = vae.get_latent_representation()
    #
    return {'adata': adata, 'vae': vae}


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


##### Load iHBCA data #####

ihbca_raw_data_dir = '/home/adr44/rds/rds-mammary-TaVjTYER47U/ADR44/data/hbca/iHBCA_raw_data/'
adata = sc.read_h5ad(ihbca_raw_data_dir + '/merged/preintegration_HBCA_inner_ENSEMBL_simple.h5ad')



#### Run scVI integration of iHBCA data ####


#input variables for scVI
batchID = "batch"
n_layers = 2 #seems best performing in testing


#setup

# batch correction needs at least 3 cells per batch:
print("Remove batches with less than 3 cells")
dict_donorID_count = Counter(adata.obs[batchID])
newDict = dict()
# identify donors with less than 3 cells.
for (key, value) in dict_donorID_count.items():
    if value < 3:
        newDict[key] = value

filter_sub = [x not in newDict.keys() for x in adata.obs[batchID]]
adata = adata[filter_sub,:]

print("Rerun HVG detection")
sc.pp.highly_variable_genes(adata, layer='lognorm', n_top_genes=5000)


#run global scVI with a few different (n_dims) parametrisations. Testing showed best results with n_layers=2

save_dir = '/home/adr44/rds/hpc-work/hbca/data/integration/scvi/output/iHBCA/global/'
os.makedirs(save_dir, exist_ok=True)

# for n_dims in [20, 50, 100, 200]:

# params from bash input
n_dims = int(sys.argv[1])

# Run scVI
print("n_dims: " + str(n_dims))
os.makedirs(save_dir + 'n_dims_' + str(n_dims), exist_ok=True)
sc.settings.figdir = save_dir + 'n_dims_' + str(n_dims) + '/'
outs = run_scVI(adata, batchID=batchID, n_dims=n_dims, n_layers=n_layers)
adata_scVI = outs['adata']
vae = outs['vae']
#save model
vae.save(save_dir + 'n_dims_' + str(n_dims) + '/vae', overwrite=True)
# vae = scvi.model.SCVI.load(save_dir + 'n_dims_' + str(n_dims) + '/vae', adata=adata_scVI)
#post processing
adata_scVI = post_dimred_processing(adata_scVI, dimred='X_scVI', run_clustering=True)
#plot training loss
plt.clf()
plt.plot(vae.history['reconstruction_loss_train']['reconstruction_loss_train'], label='train')
plt.plot(vae.history['reconstruction_loss_validation']['reconstruction_loss_validation'], label='validation')
plt.legend()
plt.savefig(save_dir + '/n_dims_' + str(n_dims) + '/training_loss.png')
#plot integration results
sc.pl.umap(adata_scVI, color=[batchID, 'dataset', 'leiden_0_5', 'leiden_1_0', 'level2'], ncols=1, save='_n_dims_' + str(n_dims) + '.png')
sc.pl.umap(adata_scVI, color=['level2'], save='_n_dims_' + str(n_dims) + '_level2.png', legend_loc='on data')
#save X_scVI, X_UMAP, and metadata
umap_df = pd.DataFrame(adata_scVI.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
umap_df.to_csv(save_dir + 'n_dims_' + str(n_dims) + '/UMAP.csv')
scVI_df = pd.DataFrame(adata_scVI.obsm['X_scVI'], columns=['scVI_' + str(i) for i in range(n_dims)])
scVI_df.to_csv(save_dir + 'n_dims_' + str(n_dims) + '/X_scVI.csv')
adata_scVI.var.to_csv(save_dir + 'n_dims_' + str(n_dims) + '/gene_data.csv')
adata_scVI.obs.to_csv(save_dir + 'n_dims_' + str(n_dims) + '/cell_data.csv')









