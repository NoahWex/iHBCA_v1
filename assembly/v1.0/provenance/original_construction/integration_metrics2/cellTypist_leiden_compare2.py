# Aim to compare the cellTypist annotations with leiden clustering results of two embeddings.



#libraries
import scanpy as sc
import celltypist
from celltypist import models

import pandas as pd
import numpy as np
import time



save_dir = '/home/adr44/rds/hpc-work/hbca/data/integration/integration_comparison2/output/leiden/'



#Load data
ihbca_raw_data_dir = '/home/adr44/rds/rds-mammary-TaVjTYER47U/ADR44/data/hbca/iHBCA_raw_data/'
print('Loading full dataset')
adata = sc.read_h5ad(ihbca_raw_data_dir + '/merged/preintegration_HBCA_inner_ENSEMBL_simple_lognorm_only.h5ad')

#Add scPoli 50 dim embedding as the best performing integration
scpoli_emb = pd.read_csv('/home/adr44/rds/hpc-work/hbca/data/integration/scPoli/output/n_dims_50/X_scVI.csv', index_col=0)
adata.obsm['X_scPoli50'] = scpoli_emb.values


#Run leiden clustering on both embeddings and PCA for a comparison
t1 = time.time()
for emb in ['scPoli50']:
    print('Computing leiden clustering for embedding: ' + emb)
    sc.pp.neighbors(adata, use_rep='X_' + emb, n_neighbors=15)
    for resolution in [0.1, 0.5, 1.0]:
        print('  resolution: ' + str(resolution))
        sc.tl.leiden(adata, resolution=resolution, key_added='leiden_' + emb + '_res' + str(resolution)) #, flavor='igraph')

print('Time taken: ' + str(time.time() - t1) + ' seconds')

#Save the leiden clustering results for comparison
adata.obs.to_csv(save_dir + 'leiden_clustering_scPoli50.csv')

