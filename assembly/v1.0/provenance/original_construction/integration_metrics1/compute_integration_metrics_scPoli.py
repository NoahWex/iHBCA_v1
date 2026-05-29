#Code to summarise and quantify the results of different iHBCA integrations
# Here we look at scPoli integration


# conda env py-scanpy

#libraries
import scanpy as sc
import scib
import pandas as pd

import os
import sys




##### Load iHBCA data #####

subset = str(sys.argv[1])
# subset = 'sub100k'

if subset == 'full':
    ihbca_raw_data_dir = '/home/adr44/rds/rds-mammary-TaVjTYER47U/ADR44/data/hbca/iHBCA_raw_data/'
    #load the subsampled data for testing
    print('Loading full dataset')
    adata = sc.read_h5ad(ihbca_raw_data_dir + '/merged/preintegration_HBCA_inner_ENSEMBL_simple.h5ad')
else:
    ihbca_raw_data_dir = '/home/adr44/rds/rds-mammary-TaVjTYER47U/ADR44/data/hbca/iHBCA_raw_data/'
    #load the subsampled data for testing
    print('Loading subset: ' + subset)
    adata = sc.read_h5ad(ihbca_raw_data_dir + '/merged/preintegration_HBCA_inner_ENSEMBL_simple_' + subset + '.h5ad') #sub100k

save_dir = '/home/adr44/rds/hpc-work/hbca/data/integration/integration_comparison/output/' + subset + '/'




##### Harmonise level1 labels to get a general cell type label nomenclaure #####
#note: I have made both lymphoid and myeloid immune cells as 'Immune'  
# and those labelled as all endothelial cells as 'Vascular endothelial' as this is the majority
preharmonise_level1_labels = ['Fibroblast', 'Basal', 'Luminal progenitor',
                  'Endothelial [vascular]', 'Luminal hormone sensing',
                  'Vascular mural', 'Lymphoid', 'Vascular Endothelial',
                  'Luminal Hormone Responsive', 'Luminal Secretory',
                  'Fibroblasts', 'Perivascular', 'Luminal Progenitor',
                  'Stroma', 'LC2', 'Immune', 'Endothelial', 'Myeloid',
                  'Luminal2', 'Luminal1', 'Mature Luminal', 'Pericytes',
                  'HRpos_Luminal', 'LC1', 'Secretory_Luminal', 'AV', 'LP',
                  'BA', 'Endothelial [lymphatic]', 'HS', 'HR',
                  'Lymphatic Endothelial', 'VasLymph', 'Epithelial',
                  'Vascular_Endothelial', 'IM', 'VA', 'EN', 'FB',
                  'Vascular_Accessory', 'Lymphocyte', 'Lymphatic',
                  'Macrophage', 'Lymphatic_Endothelial']
harmonise_level1_labels = ['Fibroblast', 'Basal-myoepithelial', 'Luminal adaptive secretory precursor',
                  'Vascular endothelial', 'Luminal hormone sensing',
                  'Perivascular', 'Immune', 'Vascular endothelial',
                  'Luminal hormone sensing', 'Luminal adaptive secretory precursor',
                  'Fibroblast', 'Perivascular', 'Luminal adaptive secretory precursor',
                  'Unlabelled', 'Lactocyte', 'Immune', 'Vascular endothelial', 'Immune',
                  'Luminal hormone sensing', 'Luminal adaptive secretory precursor', 'Luminal hormone sensing', 'Perivascular',
                  'Luminal hormone sensing', 'Lactocyte', 'Luminal adaptive secretory precursor', 'Luminal adaptive secretory precursor', 'Luminal adaptive secretory precursor',
                  'Basal-myoepithelial', 'Lymphatic endothelial', 'Luminal hormone sensing', 'Luminal hormone sensing',
                  'Lymphatic endothelial', 'Lymphatic endothelial', 'Unlabelled',
                  'Vascular endothelial', 'Immune', 'Perivascular', 'Vascular endothelial', 'Fibroblast',
                  'Perivascular', 'Immune', 'Lymphatic endothelial',
                  'Immune', 'Lymphatic endothelial']
harmonise_level1_labels_dict = dict(zip(preharmonise_level1_labels, harmonise_level1_labels))

adata.obs['harmonise_level1'] = adata.obs['level1'].map(harmonise_level1_labels_dict)







####### Compute the scPoli integration metrics #####

scpoli_dir = '/home/adr44/rds/hpc-work/hbca/data/integration/scPoli/output/'

latent_dim = int(sys.argv[2])
print(latent_dim)

# #remove the counts matrixes to save RAM
# adata.layers = {}
# adata.X = None

if (latent_dim == 20) | (latent_dim == 50):
    hidden_layer_size = 128
else:
    hidden_layer_size = 512

print("latent_dim: " + str(latent_dim))
print("hidden_layer_size: " + str(hidden_layer_size))
#add scpoli data
adata_integration = adata.copy()
scpoli_cell_data = pd.read_csv(scpoli_dir + '/n_dims_' + str(latent_dim) + '/cell_data.csv', index_col=0, dtype=str)
cells_to_keep = scpoli_cell_data.index.intersection(adata_integration.obs_names)
cells_to_keep_bool = scpoli_cell_data.index.isin(cells_to_keep)
adata_integration.obs = scpoli_cell_data.loc[adata_integration.obs_names,]

#reorder adata and adata_integration to match scPoli order
adata_integration = adata_integration[cells_to_keep,]
adata = adata[cells_to_keep,]

#add scpoli embedding
scpoli_emb = pd.read_csv(scpoli_dir + '/n_dims_' + str(latent_dim) + '/X_scVI.csv', index_col=0)
adata_integration.obsm['X_emb'] = scpoli_emb.iloc[cells_to_keep_bool,].values

#Make scib metadata into categorical for the metrics
adata_integration.obs['harmonise_level1'] = adata_integration.obs['level1'].map(harmonise_level1_labels_dict)
adata_integration.obs['harmonise_level1'] = adata_integration.obs['harmonise_level1'].astype('category')
adata_integration.obs['leiden_0_5'] = adata_integration.obs['leiden_0_5'].astype('category')
adata_integration.obs['batch'] = adata_integration.obs['dataset'].astype(str) + '_' + adata_integration.obs['batch'].astype(str)
adata_integration.obs['batch'] = adata_integration.obs['batch'].astype('category')

#similarly for the original adata
adata.obs['harmonise_level1'] = adata_integration.obs['harmonise_level1']
adata.obs['leiden_0_5'] = adata_integration.obs['leiden_0_5']
adata.obs['batch'] = adata_integration.obs['batch']

# plot some UMAPs to check data merging
sc.settings.figdir = save_dir + 'test_umaps/'
sc.pl.umap(adata, color=['harmonise_level1', 'leiden_0_5', 'batch'], ncols=1, save='_original.png')
sc.pl.umap(adata_integration, color=['harmonise_level1', 'leiden_0_5', 'batch'], ncols=1, save='_scPoli.png')
sc.pl.embedding(adata_integration, basis='X_emb', color=['harmonise_level1'], ncols=1, save='_scPoli_emb.png')

#Create KNN graphs on the subsetted objects
sc.pp.neighbors(adata, n_neighbors=15, use_rep='X_pca')
sc.pp.neighbors(adata_integration, n_neighbors=15, use_rep='X_emb')

#Integration metrics
integration_metrics = scib.metrics.metrics(adata=adata, 
                    adata_int=adata_integration, 
                    batch_key='batch', 
                    label_key='harmonise_level1', 
                    embed='X_emb',
                    cluster_key='leiden_compare',
                    ari_=True, nmi_=True, nmi_method='arithmetic', nmi_dir=None, silhouette_=True, si_metric='euclidean', 
                    pcr_=True, cell_cycle_=False, hvg_score_=False, 
                    isolated_labels_=True, isolated_labels_f1_=True, isolated_labels_asw_=True, n_isolated=None, 
                    graph_conn_=True, trajectory_=False, kBET_=True, lisi_graph_=True,
                    organism='human', 
                    n_cores=8,
                    type_='emb')
#Save the metrics
os.makedirs(save_dir + 'scPoli_metrics/', exist_ok=True)
integration_metrics.to_csv(save_dir + 'scPoli_metrics/' + 'n_dims_' + str(latent_dim) + '_integration_metrics.csv')
# integration_metrics = pd.read_csv(save_dir + 'scPoli_metrics/' + 'n_dims_' + str(latent_dim) + '_integration_metrics.csv')


