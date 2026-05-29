# Aim to compare the cellTypist annotations with leiden clustering results of two embeddings.



#libraries
import scanpy as sc
import celltypist
from celltypist import models

import pandas as pd
import numpy as np
import time
import os



save_dir = '/home/adr44/rds/hpc-work/hbca/data/integration/integration_comparison2/output/leiden/'



#Load data
ihbca_raw_data_dir = '/home/adr44/rds/rds-mammary-TaVjTYER47U/ADR44/data/hbca/iHBCA_raw_data/'
print('Loading full dataset')
adata = sc.read_h5ad(ihbca_raw_data_dir + '/merged/preintegration_HBCA_inner_ENSEMBL_simple_lognorm_only.h5ad')

# pca already computed in the original data (X_pca)

#Run leiden clustering on both embeddings and PCA for a comparison
t1 = time.time()
for emb in ['pca']: #'scvi20', 'scPoli50', 
    print('Computing leiden clustering for embedding: ' + emb)
    sc.pp.neighbors(adata, use_rep='X_' + emb, n_neighbors=15)
    for resolution in [0.1, 0.5, 1.0]:
        print('  resolution: ' + str(resolution))
        #run only if the clustering result does not already exist
        if os.path.exists(save_dir + 'leiden_clustering_' + emb + '_res' + str(resolution) + '.csv'):
            print('  Clustering result already exists, skipping computation')
        else:
            sc.tl.leiden(adata, resolution=resolution, key_added='leiden_' + emb + '_res' + str(resolution)) #, flavor='igraph')
            temp = adata.obs['leiden_' + emb + '_res' + str(resolution)]
            temp.to_csv(save_dir + 'leiden_clustering_' + emb + '_res' + str(resolution) + '.csv')

print('Time taken: ' + str(time.time() - t1) + ' seconds')

