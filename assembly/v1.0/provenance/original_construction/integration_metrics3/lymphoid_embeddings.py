#Compare the level1.5 annotations on the scVI and scPoli data of lymphoid cells

#libraries
import scanpy as sc
import pandas as pd
import numpy as np
import os

import seaborn as sns
import matplotlib.pyplot as plt

#Save dirs and data
save_dir = '/home/adr44/rds/hpc-work/hbca/data/integration/integration_comparison3/output/lymphoid_embeddings/'
level1_dir = '/home/adr44/rds/hpc-work/hbca/data/integration/ihbca_annotation/output/level1/'
level1_5_dir = '/home/adr44/rds/hpc-work/hbca/data/integration/ihbca_annotation/output/level1.5/'
ihbca_raw_data_dir = '/home/adr44/rds/rds-mammary-TaVjTYER47U/ADR44/data/hbca/iHBCA_raw_data/'

celltype = 'TL'
adata_lym = sc.read_h5ad(level1_dir + 'adata/ihbca_level1_' + celltype + '_cells.h5ad')
adata_all = sc.read_h5ad(ihbca_raw_data_dir + '/merged/preintegration_HBCA_inner_ENSEMBL_simple_lognorm_only.h5ad')

scpoli_emb = pd.read_csv('/home/adr44/rds/hpc-work/hbca/data/integration/scPoli/output/n_dims_50/X_scVI.csv', index_col=0)
scpoli_emb = scpoli_emb.loc[adata_all.obs_names.isin(adata_lym.obs_names), :]
adata_lym.obsm['X_scPoli50'] = scpoli_emb.values

scpoli_umap = pd.read_csv('/home/adr44/rds/hpc-work/hbca/data/integration/scPoli/output/n_dims_50/UMAP.csv', index_col=0)
scpoli_umap = scpoli_umap.loc[adata_all.obs_names.isin(adata_lym.obs_names), :]
adata_lym.obsm['X_umap_scPoli50'] = scpoli_umap.values

#Add level1.5 annotations
level1_5_annotations = pd.read_csv(level1_5_dir + 'adata/' + celltype + '_level1.5_annotation.csv', index_col=0)
adata_lym.obs['level1.5_annotation'] = level1_5_annotations.loc[adata_lym.obs_names, 'level1.5_annotation'].values

#Generate subUMAPs for both embeddings
sc.settings.figdir = save_dir
adata_lym.obsm['X_scVI100'] = adata_lym.obsm['X_scVI']
for emb in ['scVI100', 'scPoli50']:
    print('Generating UMAP for embedding: ' + emb)
    sc.pp.neighbors(adata_lym, use_rep='X_' + emb, n_neighbors=15)
    sc.tl.umap(adata_lym, random_state=42) # key_added='X_umap_' + emb + '_sub'
    adata_lym.obsm['X_umap_' + emb + '_sub'] = adata_lym.obsm['X_umap']

#Plot UMAPs coloured by level1.5 annotation
level1_5_colour_dictionary = {
    'Lactocyte': "#AC79F3",
    'Luminal adaptive secretory precursor': '#DDA0DD',
    'Luminal hormone sensing': '#EE8298',
    'Basal-myoepithelial': '#E6554A',
    'Fibroblast': '#994E2E',
    'Perivascular': '#F89440',
    'Vascular endothelial': '#F1C232',
    'Lymphatic endothelial': "#E4ED71",
    'CD4T': '#99A800',
    'CD8T': '#64C6A6',
    'NK/NKT': "#66B505",
    'B-cell': "#188048",
    'Plasma cell': "#2379558E",
    'Macrophage': "#45D5E2",
    'Macrophage lipid-associated': "#72C4E5",
    'Monocyte': "#4586EE",
    'Conventional dendritic cell': "#092C8C",
    'Plasmacytoid dendritic cell': "#77A7FA",
    'Mast cell': "#0048FF",
    'Neutrophil': "#370587",  
    'Single nuclei': '#EAEAEA',  # Light grey
    'Doublet': '#EAEAEA'  # Light grey
}

for emb in ['scVI100', 'scPoli50']:
    #Plot UMAP coloured by level1.5 annotation
    sc.pl.embedding(adata_lym, basis='X_umap_' + emb + '_sub', 
                    color='level1.5_annotation', palette=level1_5_colour_dictionary,
                    title='Lymphoid cells UMAP (' + emb + ')', 
               save='_lymphoid_level1.5_' + emb + '.png', show=False)
    sc.pl.embedding(adata_lym, basis='X_umap_' + emb + '_sub', 
                    color='level1.5_annotation', palette=level1_5_colour_dictionary,
                    title='Lymphoid cells UMAP (' + emb + ')', 
               save='_lymphoid_level1.5_' + emb + '.pdf', show=False)


