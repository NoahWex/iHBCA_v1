#Code to generate some level1 annotation confusion matrices between the different dataset annotations

#libraries
import scanpy as sc
import pandas as pd
import os

import seaborn as sns
import matplotlib.pyplot as plt



#directories and data
save_dir = '/home/adr44/rds/hpc-work/hbca/data/integration/integration_comparison3/output/level1_confusion_mat/'

# #temp 
# adata = sc.read_h5ad('/home/adr44/rds/rds-mammary-TaVjTYER47U/ADR44/data/hbca/DATA_FINAL/cellxgene/HBCA_cellxgene_global.h5ad')
# adata.obs.to_csv('/home/adr44/rds/rds-mammary-TaVjTYER47U/ADR44/data/hbca/DATA_FINAL/cellxgene/HBCA_cellxgene_global_obs.csv')

HBCA_obs = pd.read_csv('/home/adr44/rds/rds-mammary-TaVjTYER47U/ADR44/data/hbca/DATA_FINAL/cellxgene/HBCA_cellxgene_global_obs.csv', index_col=0)
ihbca_level1_obs = pd.read_csv('/home/adr44/rds/hpc-work/hbca/data/integration/ihbca_annotation/output/level1/adata/ihbca_level1_annotations.csv', index_col=0)

#fix iHBCA level1 annotation for the HBCA ('reed') dataset. The old annotations are currently in ihbca_level1_obs
#reorder HBCA_obs to match ihbca_level1_obs
HBCA_obs = HBCA_obs.loc[ihbca_level1_obs.index[ihbca_level1_obs['dataset'] == 'reed']]
HBCA_obs.level1[HBCA_obs.level1 == 'Luminal adaptive secretory precurser'] = 'Luminal adaptive secretory precursor'
ihbca_level1_obs.loc[HBCA_obs.index, 'level1'] = HBCA_obs.loc[:, 'level1']

#Perdataset order of labels
dataset_level1_orders = {
    'ihbca': ['Lactocyte', 'Luminal adaptive secretory precursor', 'Luminal hormone sensing', 'Basal-myoepithelial',
             'Fibroblast', 'Perivascular', 'Vascular endothelial', 'Lymphatic endothelial',
             'T-lymphocyte', 'B-lymphocyte',
             'Myeloid'],
    'reed': ['Luminal adaptive secretory precursor', 'Luminal hormone sensing', 'Basal-myoepithelial',
             'DDC1', 'DDC2',
             'Fibroblast', 'Perivascular', 'Vascular endothelial', 'Lymphatic endothelial',
             'Lymphoid',
             'Myeloid',
             'Doublet', 'stripped_nuclei'],
    'kumar': ['Luminal Secretory', 'Luminal Hormone Responsive', 'Basal',
              'Fibroblast', 'Perivascular', 'Vascular Endothelial', 'Lymphatic Endothelial',
              'Lymphoid',
              'Myeloid'],
    'nee': ['Luminal1', 'Luminal2', 'Basal',
            'Fibroblasts', 'Pericytes', 'Endothelial', 'Lymphatic',
            'Immune'],
    'pal': ['Luminal Progenitor', 'Mature Luminal', 'Basal',
            'Epithelial',
            'Fibroblast', 'Stroma'],
    'twigger': ['LC1', 'LC2', 'LP', 'HR', 'BA',
                'FB', 'VA', 'EN',
                'IM'],
    'murrow': ['Secretory_Luminal', 'HRpos_Luminal', 'Basal', 
               'Fibroblast', 'Vascular_Accessory',  'Vascular_Endothelial', 'Lymphatic_Endothelial',
               'Lymphocyte',
               'Macrophage'],
    'gray': ['AV', 'HS', 'BA', 
             'Fibroblast', 'VasLymph', 
             'Immune']}

# For each dataset create confusion matrix (level1 vs level1_annotation)
os.makedirs(save_dir + 'csv/', exist_ok=True)
os.makedirs(save_dir + 'heatmap/', exist_ok=True)
for dataset in ['reed', 'kumar', 'nee', 'pal', 'twigger', 'murrow', 'gray']:
    print('Processing dataset: ' + dataset)
    dataset_level1_obs = ihbca_level1_obs[ihbca_level1_obs['dataset'] == dataset]
    confusion_mat = pd.crosstab(dataset_level1_obs['level1_annotation'], dataset_level1_obs['level1'], rownames=['iHBCA Level 1 Annotation'], colnames=[dataset + ' Level 1 Annotation'])
    #reorder confusion matrix
    confusion_mat = confusion_mat.reindex(index=dataset_level1_orders['ihbca'], columns=dataset_level1_orders[dataset], fill_value=0)
    confusion_mat.to_csv(save_dir + 'csv/confusion_matrix_' + dataset + '_level1.csv')
    #normalize confusion matrix
    confusion_mat_norm = confusion_mat.astype(int)
    confusion_mat_norm = confusion_mat_norm.div(confusion_mat_norm.sum(axis=1), axis=0)
    #plot heatmap (normalized)
    plt.figure(figsize=(12,10))
    sns.heatmap(confusion_mat_norm, annot=True, cmap='Blues') 
    plt.title('Confusion Matrix: iHBCA Level 1 vs ' + dataset + ' Level 1')
    plt.tight_layout()
    plt.savefig(save_dir + 'heatmap/confusion_matrix_' + dataset + '_level1_heatmap.png')
    plt.savefig(save_dir + 'heatmap/confusion_matrix_' + dataset + '_level1_heatmap.pdf')
    plt.close()
    #plot heatmap (unnormalized)
    plt.figure(figsize=(12,10))
    sns.heatmap(confusion_mat, annot=True, cmap='Blues') 
    plt.title('Confusion Matrix: iHBCA Level 1 vs ' + dataset + ' Level 1 (raw counts)')
    plt.tight_layout()
    plt.savefig(save_dir + 'heatmap/confusion_matrix_' + dataset + '_level1_heatmap_raw.png')
    plt.savefig(save_dir + 'heatmap/confusion_matrix_' + dataset + '_level1_heatmap_raw.pdf')


