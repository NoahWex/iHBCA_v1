#Code to summarise and quantify the results of different iHBCA integrations
# Here we look at no integration as a control


# conda env py-scanpy

#libraries
import scanpy as sc
import scib
import pandas as pd

import os
import sys



##### Load iHBCA data #####

subset = str(sys.argv[1])
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






##### Compute the no integration metrics #####

#Create KNN graphs on the subsetted objects
sc.pp.neighbors(adata, n_neighbors=15, use_rep='X_pca')

# Create adata_integration object fro comparison
adata_integration = adata.copy()

#Make scib metadata into categorical for the metrics
adata_integration.obs['harmonise_level1'] = adata_integration.obs['level1'].map(harmonise_level1_labels_dict)
adata_integration.obs['harmonise_level1'] = adata_integration.obs['harmonise_level1'].astype('category')
adata_integration.obs['leiden_0_5'] = adata_integration.obs['leiden_0_5'].astype('category')
adata_integration.obs['batch'] = adata_integration.obs['dataset'].astype(str) + '_' + adata_integration.obs['batch'].astype(str)
adata_integration.obs['batch'] = adata_integration.obs['batch'].astype('category')

adata.obs['harmonise_level1'] = adata_integration.obs['harmonise_level1']
adata.obs['leiden_0_5'] = adata_integration.obs['leiden_0_5']
adata.obs['batch'] = adata_integration.obs['batch']

#Integration metrics
# Cell cycle fails due to missing genes so do not perform
# kBET fails due to missing packages which are not easily installed to this environment
# hvg_score fails for some reason...? https://github.com/theislab/scib/issues/125
# trajectory fails (likely as dpt needs to be run first) but I don't think it is relevant here
# lisi_graph fails due to /lib64/libc.so.6: version `GLIBC_2.34' not found (tied to OS so cannot fix)
    # FIXED USING:  https://github.com/theislab/scib/issues/375
integration_metrics = scib.metrics.metrics(adata=adata, 
                        adata_int=adata_integration, 
                        batch_key='batch', 
                        label_key='harmonise_level1', 
                        embed='X_pca',
                        cluster_key='leiden_compare',
                        ari_=True, nmi_=True, nmi_method='arithmetic', nmi_dir=None, silhouette_=True, si_metric='euclidean', 
                        pcr_=True, cell_cycle_=False, hvg_score_=False, 
                        isolated_labels_=True, isolated_labels_f1_=True, isolated_labels_asw_=True, n_isolated=None, 
                        graph_conn_=True, trajectory_=False, kBET_=True, lisi_graph_=True,
                        organism='human', 
                        n_cores=8,
                        type_='emb')

#tester
integration_metrics = scib.metrics.metrics(adata=adata, 
                        adata_int=adata_integration, 
                        batch_key='batch', 
                        label_key='harmonise_level1', 
                        embed='X_pca',
                        cluster_key='leiden_compare',
                        ari_=False, nmi_=False, nmi_method='arithmetic', nmi_dir=None, silhouette_=False, si_metric='euclidean', 
                        pcr_=False, cell_cycle_=False, hvg_score_=False, 
                        isolated_labels_=False, isolated_labels_f1_=False, isolated_labels_asw_=False, n_isolated=None, 
                        graph_conn_=False, trajectory_=False, kBET_=True, lisi_graph_=False,
                        organism='human', 
                        n_cores=8,
                        type_='emb')


#Save the metrics
os.makedirs(save_dir + 'unintegrated_metrics/', exist_ok=True)
integration_metrics.to_csv(save_dir + 'unintegrated_metrics/' + 'integration_metrics.csv')




