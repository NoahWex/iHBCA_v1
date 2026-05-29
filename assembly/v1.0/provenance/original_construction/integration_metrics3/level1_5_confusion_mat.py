#Code to generate some level1 annotation confusion matrices between the different dataset annotations

#libraries
import scanpy as sc
import pandas as pd
import os

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm



#directories and data
save_dir = '/home/adr44/rds/hpc-work/hbca/data/integration/integration_comparison3/output/level1_5_confusion_mat/'

# #temp 
# adata = sc.read_h5ad('/home/adr44/rds/rds-mammary-TaVjTYER47U/ADR44/data/hbca/DATA_FINAL/cellxgene/HBCA_cellxgene_global.h5ad')
# adata.obs.to_csv('/home/adr44/rds/rds-mammary-TaVjTYER47U/ADR44/data/hbca/DATA_FINAL/cellxgene/HBCA_cellxgene_global_obs.csv')

HBCA_obs = pd.read_csv('/home/adr44/rds/rds-mammary-TaVjTYER47U/ADR44/data/hbca/DATA_FINAL/cellxgene/HBCA_cellxgene_global_obs.csv', index_col=0)
ihbca_level1_5_obs = pd.read_csv('/home/adr44/rds/hpc-work/hbca/data/integration/ihbca_annotation/output/level1.5/adata/ihbca_level1.5_annotations.csv', index_col=0)

#fix iHBCA level1 annotation for the HBCA ('reed') dataset. The old annotations are currently in ihbca_level1_obs
#reorder HBCA_obs to match ihbca_level1_5_obs
HBCA_obs = HBCA_obs.loc[ihbca_level1_5_obs.index[ihbca_level1_5_obs['dataset'] == 'reed']]
HBCA_obs.level1[HBCA_obs.level1 == 'Luminal adaptive secretory precurser'] = 'Luminal adaptive secretory precursor'
ihbca_level1_5_obs.loc[HBCA_obs.index, 'level1'] = HBCA_obs.loc[:, 'level1']
ihbca_level1_5_obs.loc[HBCA_obs.index, 'level2'] = HBCA_obs.loc[:, 'level2']

#Perdataset order of labels
dataset_level1_5_orders = {
    'ihbca': ['Lactocyte', 'Luminal adaptive secretory precursor', 'Luminal hormone sensing', 'Basal-myoepithelial',
                'Fibroblast', 'Perivascular', 'Vascular endothelial', 'Lymphatic endothelial',
                'CD4T', 'CD8T', 'NK/NKT', 'B-cell', 'Plasma cell', 
                'Macrophage', 'Macrophage lipid-associated', 'Monocyte', 'Conventional dendritic cell', 'Plasmacytoid dendritic cell',
                'Mast cell', 'Neutrophil'],
    'reed': ['Luminal adaptive secretory precursor', 'Luminal hormone sensing', 'Basal-myoepithelial',
             'DDC1', 'DDC2',
             'Fibroblast', 'Perivascular', 'Vascular endothelial', 'Lymphatic endothelial',
             'Lymphoid',
             'Myeloid'],
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
               'Fibroblast',  'Vascular_Accessory', 'Vascular_Endothelial', 'Lymphatic_Endothelial',
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
    dataset_level1_5_obs = ihbca_level1_5_obs[ihbca_level1_5_obs['dataset'] == dataset]
    confusion_mat = pd.crosstab(dataset_level1_5_obs['level1.5_annotation'], dataset_level1_5_obs['level1'], rownames=['iHBCA Level 1.5 Annotation'], colnames=[dataset + ' Level 1 Annotation'])
    #reorder confusion matrix
    confusion_mat = confusion_mat.reindex(index=dataset_level1_5_orders['ihbca'], columns=dataset_level1_5_orders[dataset], fill_value=0)
    confusion_mat.to_csv(save_dir + 'csv/confusion_matrix_' + dataset + '_level1.csv')
    #normalize confusion matrix
    confusion_mat_norm = confusion_mat.astype(int)
    confusion_mat_norm = confusion_mat_norm.div(confusion_mat_norm.sum(axis=1), axis=0)
    #plot heatmap (normalized)
    plt.figure(figsize=(12,10))
    sns.heatmap(confusion_mat_norm, annot=True, cmap='Blues') #, cmap='viridis'
    plt.title('Confusion Matrix: iHBCA Level 1 vs ' + dataset + ' Level 1.5')
    plt.tight_layout()
    plt.savefig(save_dir + 'heatmap/confusion_matrix_' + dataset + '_level1.5_heatmap.png')
    plt.savefig(save_dir + 'heatmap/confusion_matrix_' + dataset + '_level1.5_heatmap.pdf')
    plt.close()
    #plot heatmap (unnormalized)
    plt.figure(figsize=(12,10))
    sns.heatmap(confusion_mat, annot=True, cmap='Blues') #, cmap='viridis'
    plt.title('Confusion Matrix: iHBCA Level 1 vs ' + dataset + ' Level 1.5 (raw counts)')
    plt.tight_layout()
    plt.savefig(save_dir + 'heatmap/confusion_matrix_' + dataset + '_level1.5_heatmap_raw.png')
    plt.savefig(save_dir + 'heatmap/confusion_matrix_' + dataset + '_level1.5_heatmap_raw.pdf')

#Do the same but for level2 annotations from each dataset (use level1 in case of NA in level2)
ihbca_level1_5_obs['level2_filled'] = ihbca_level1_5_obs['level2'].copy()
ihbca_level1_5_obs['level2_filled'][ihbca_level1_5_obs['level2_filled'].isna()] = ihbca_level1_5_obs['level1'][ihbca_level1_5_obs['level2_filled'].isna()]

order_of_celltypes_list = {'reed': ["LASP1", "LASP2", "LASP3", "LASP4", "LASP5",
                                           "LHS1", "LHS2", "LHS3",
                                           "BMYO1", "BMYO2", 
                                           "FB1", "FB2", "FB3", "FB4",
                                           "PV1", "PV2", "PV3", "PV4", "PV5",
                                           "VEV", "VEC", "VEA","VEAT", 
                                           "LE1", "LE2",
                                           'CD4_naive', "CD4_Th", 
                                           'CD8_Tem', 'CD8_Trm',
                                           'CD8_Tc1',
                                           'NKT', 'NK',
                                           'ILC', 
                                           'B_naive', 'B_mem_switched', 'B_mem_unswitched', 
                                           'Plasma_cell', 
                                           'Macro', 'Macro-lipo', 'DC'],
                                'kumar': ['Lumsec-major','Lumsec-basal', 'Lumsec-myo', 
                                            'Lumsec-KIT', 'Lumsec-lac', 'Lumsec-HLA', 
                                            'Lumsec-prol',
                                            'LummHR-major', 'LummHR-active', 'LummHR-SCGB',
                                            'basal', 
                                            'pericytes', 'vsmc',
                                            'Fibro-major', 'Fibro-matrix', 'Fibro-prematrix',
                                            'Fibro-SFRP4',
                                            'Vas-arterial', 'Vas-capillary', 'Vas-venous',
                                            'Lymph-major', 'Lymph-immune', 'Lymph-valve1',
                                            'Lymph-valve2', 
                                            'CD4-activated', 'CD4-naive', 'CD4-Tem', 
                                            'CD4-Th', 'CD4-Th-like', 'CD4-Treg', 
                                            'CD8-activated', 'CD8-Tem', 'CD8-Trm',
                                            'T_prol', 'GD',
                                            'NK', 'NKT', 'NK-ILCs',
                                            'b_naive', 'bmem_switched', 'bmem_unswitched', 
                                            'plasma_IgA', 'plasma_IgG', 
                                            'Macro-IFN', 'Macro-lipo', 'Macro-m1', 
                                            'Macro-m1-CCL', 'Macro-m2', 'Macro-m2-CXCL', 
                                            'Mono-classical', 'Mono-non-classical',
                                            'mDC', 'pDC', 'cDC1', 'cDC2',  
                                            'Mast', 
                                            'Neutrophil',
                                            'mye-prol'],
                                'nee': ['Luminal1-ALDH1A3', 'Luminal1-LTF',
                                          'Luminal2-AREG', 'Luminal2-MUCL1',
                                          'Basal', 'Basal-Myoepithelial',
                                          'Pericytes',
                                          'Fibroblasts',
                                          'Endothelial',
                                          'Lymphatic',
                                          'Immune'],
                                'gray': ['AP', 'BL', 
                                           'HSa', 'HSb', 
                                           'BAa', 'BAb',
                                           'F1', 'F2', 'F3', 
                                           'VL3_PE', 'VL2_VE', 'VL1_LE',
                                           'I3_Tcell', 'I2_NK',
                                           'I4_Bcell', 'I5_PlasmaCell', 
                                           'I1_Myeloid'],
                                'murrow': ['Secretory_Luminal', 'HRpos_Luminal', 'Basal', 
                                             'Fibroblast',
                                             'Vascular_Accessory',
                                             'Vascular_Endothelial', 'Lymphatic_Endothelial',
                                             'Lymphocyte', 
                                             'Macrophage'],
                                'twigger': ['LC1', 'LC2', 'LP', 'HR', 'BA', 
                                              'FB',
                                              'VA',
                                              'EN',
                                              'IM'],
                                'pal': ['Epithelial',
                                          'Luminal Progenitor', 'Mature Luminal', 'Basal', 
                                          'Fibroblast',
                                          'Stroma']}


for dataset in ['reed', 'kumar', 'nee', 'pal', 'twigger', 'murrow', 'gray']:
        print('Processing dataset: ' + dataset)
        dataset_level1_5_obs = ihbca_level1_5_obs[ihbca_level1_5_obs['dataset'] == dataset]
        confusion_mat = pd.crosstab(dataset_level1_5_obs['level1.5_annotation'], dataset_level1_5_obs['level2_filled'], rownames=['iHBCA Level 1.5 Annotation'], colnames=[dataset + ' Level 2 Annotation'])
        #reorder confusion matrix
        confusion_mat = confusion_mat.reindex(index=dataset_level1_5_orders['ihbca'], columns=order_of_celltypes_list[dataset], fill_value=0)
        confusion_mat.to_csv(save_dir + 'csv/confusion_matrix_' + dataset + '_level2.csv')
        #normalize confusion matrix
        confusion_mat_norm = confusion_mat.astype(int)
        confusion_mat_norm = confusion_mat_norm.div(confusion_mat_norm.sum(axis=1), axis=0)
        #plot heatmap (normalized)
        plt.figure(figsize=(30,14))
        sns.heatmap(confusion_mat_norm, annot=True, cmap='Blues') #, cmap='viridis'
        plt.title('Confusion Matrix: iHBCA Level 1.5 vs ' + dataset + ' Level 2')
        plt.tight_layout()
        plt.savefig(save_dir + 'heatmap/confusion_matrix_' + dataset + '_level2_heatmap.png')
        plt.savefig(save_dir + 'heatmap/confusion_matrix_' + dataset + '_level2_heatmap.pdf')
        plt.close()
        #plot heatmap (unnormalized)
        plt.figure(figsize=(30,14))
        sns.heatmap(confusion_mat, annot=True, cmap='Blues',
                    norm=LogNorm(
                                vmin=1,
                                vmax=confusion_mat.values.max()
                        ))
        plt.title('Confusion Matrix: iHBCA Level 1.5 vs ' + dataset + ' Level 2 (raw counts)')
        plt.tight_layout()
        plt.savefig(save_dir + 'heatmap/confusion_matrix_' + dataset + '_level2_heatmap_raw.png')
        plt.savefig(save_dir + 'heatmap/confusion_matrix_' + dataset + '_level2_heatmap_raw.pdf')
