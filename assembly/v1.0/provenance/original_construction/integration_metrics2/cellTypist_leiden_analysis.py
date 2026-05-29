# Load Leiden clustering results and compare with cellTypist annotations

#libraries
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.spatial.distance import jensenshannon


import os

#directories
save_dir = '/home/adr44/rds/hpc-work/hbca/data/integration/integration_comparison2/output/leiden_analysis/'
leiden_dir = '/home/adr44/rds/hpc-work/hbca/data/integration/integration_comparison2/output/leiden/'

#Load adata
print('Loading data...')
ihbca_raw_data_dir = '/home/adr44/rds/rds-mammary-TaVjTYER47U/ADR44/data/hbca/iHBCA_raw_data/'
adata = sc.read_h5ad(ihbca_raw_data_dir + '/merged/preintegration_HBCA_inner_ENSEMBL_simple_lognorm_only.h5ad')

# #set up (old leiden results saved in a different format)
# pca = pd.read_csv('/home/adr44/rds/hpc-work/hbca/data/integration/integration_comparison2/output/leiden/leiden_clustering_pca.csv', index_col=0)
# scPoli = pd.read_csv('/home/adr44/rds/hpc-work/hbca/data/integration/integration_comparison2/output/leiden/leiden_clustering_scPoli50.csv', index_col=0)
# for res in [0.1, 0.5, 1.0]:
#     temp1 = pca['leiden_pca_res' + str(res)]
#     temp1.to_csv(leiden_dir + 'leiden_clustering_pca_res' + str(res) + '.csv')
#     temp2 = scPoli['leiden_scPoli50_res' + str(res)]
#     temp2.to_csv(leiden_dir + 'leiden_clustering_scPoli50_res' + str(res) + '.csv')

#load leiden clustering results
for emb in ['pca', 'scVI100', 'scPoli50']:
    for resolution in [0.1, 0.5, 1.0]:
        leiden_file = leiden_dir + 'leiden_clustering_' + emb + '_res' + str(resolution) + '.csv'
        if os.path.exists(leiden_file):
            print('Loading leiden clustering for embedding: ' + emb + ' at resolution: ' + str(resolution))
            leiden_labels = pd.read_csv(leiden_file, index_col=0)
            adata.obs['leiden_' + emb + '_res' + str(resolution)] = leiden_labels.values
            adata.obs['leiden_' + emb + '_res' + str(resolution)] = adata.obs['leiden_' + emb + '_res' + str(resolution)].astype('category')
        else:
            print('Leiden clustering file not found for embedding: ' + emb + ' at resolution: ' + str(resolution))

#Load cellTypist annotations
celltypist_reed = pd.read_csv('/home/adr44/rds/hpc-work/hbca/data/integration/integration_comparison2/output/celltypist/celltypist_reed_predicted_labels.csv', index_col=0)
celltypist_kumar = pd.read_csv('/home/adr44/rds/hpc-work/hbca/data/integration/integration_comparison2/output/celltypist/celltypist_kumar_predicted_labels.csv', index_col=0)
adata.obs['cellTypist_annotation_reed'] = celltypist_reed.predicted_labels.values
adata.obs['cellTypist_annotation_kumar'] = celltypist_kumar.predicted_labels.values

#Load harmonised level1 annotations
CAP_df = pd.read_csv('/home/adr44/rds/rds-mammary-TaVjTYER47U/ADR44/data/hbca/DATA_FINAL/cellxgene/integrated_HBCA_cellxgene_CAP_metadata.csv', index_col=0)


#load the scVI100 and scPoli50 UMAP embeddings for plotting
scvi_umap = pd.read_csv('/home/adr44/rds/hpc-work/hbca/data/integration/scvi/output/iHBCA/global/n_dims_100/UMAP.csv', index_col=0)
adata.obsm['X_umap_scVI100'] = scvi_umap.values
scpoli_umap = pd.read_csv('/home/adr44/rds/hpc-work/hbca/data/integration/scPoli/output/n_dims_50/UMAP.csv', index_col=0)
adata.obsm['X_umap_scPoli50'] = scpoli_umap.values


# Basic plots
sc.settings.figdir = save_dir + '/umap'
os.makedirs(sc.settings.figdir, exist_ok=True)
colour_list = ['leiden_pca_res0.1', 'leiden_pca_res0.5', 'leiden_pca_res1.0',
               'leiden_scVI100_res0.1', 'leiden_scVI100_res0.5', 'leiden_scVI100_res1.0',
               'leiden_scPoli50_res0.1', 'leiden_scPoli50_res0.5', 'leiden_scPoli50_res1.0']
#limit to those that exist in the data
colour_list = [x for x in colour_list if x in adata.obs.columns]
sc.pl.embedding(adata, basis='X_umap', color=colour_list, legend_loc='on data', save='_leiden_pca.png')
sc.pl.embedding(adata, basis='X_umap_scVI100', color=colour_list, legend_loc='on data', save='_leiden_scVI100.png')
sc.pl.embedding(adata, basis='X_umap_scPoli50', color=colour_list, legend_loc='on data', save='_leiden_scPoli50.png')

sc.pl.embedding(adata, basis='X_umap', color=['cellTypist_annotation_reed', 'cellTypist_annotation_kumar'], legend_loc='on data', save='_cellTypist.png')
sc.pl.embedding(adata, basis='X_umap_scVI100', color=['cellTypist_annotation_reed', 'cellTypist_annotation_kumar'], legend_loc='on data', save='_cellTypist_scVI100.png')
sc.pl.embedding(adata, basis='X_umap_scPoli50', color=['cellTypist_annotation_reed', 'cellTypist_annotation_kumar'], legend_loc='on data', save='_cellTypist_scPoli50.png')


# Make confusion map of leiden clusters vs cellTypist annotations
os.makedirs(save_dir + 'heatmap/', exist_ok=True)
for emb in ['pca', 'scVI100', 'scPoli50']:
    for resolution in [0.1, 0.5, 1.0]:
        key = 'leiden_' + emb + '_res' + str(resolution)
        if key in adata.obs.columns:
            print('Making heatmap for embedding: ' + emb + ' at resolution: ' + str(resolution))
            ct_vs_leiden = pd.crosstab(adata.obs['cellTypist_annotation_reed'], adata.obs[key], normalize='index')
            plt.figure(figsize=(10,8))
            sns.heatmap(ct_vs_leiden, cmap='Reds')
            plt.title('CellTypist vs Leiden Clusters (' + emb + ' res ' + str(resolution) + ')')
            plt.savefig(save_dir + 'heatmap/heatmap_cellTypist_reed_' + emb + '_res' + str(resolution) + '.png')
            plt.close()
            ct_vs_leiden = pd.crosstab(adata.obs['cellTypist_annotation_kumar'], adata.obs[key], normalize='index')
            plt.figure(figsize=(10,8))
            sns.heatmap(ct_vs_leiden, cmap='Reds')
            plt.title('CellTypist vs Leiden Clusters (' + emb + ' res ' + str(resolution) + ')')
            plt.savefig(save_dir + 'heatmap/heatmap_cellTypist_kumar_' + emb + '_res' + str(resolution) + '.png')
            plt.close()



# Generate a function that calculates the absolute variation from expected distributions by a given metadata label and a clustering
def calculate_variation_from_expected(adata, cluster_key, label_key, ignored_labels=[]):
    contingency_table = pd.crosstab(adata.obs[cluster_key], adata.obs[label_key])
    contingency_table = contingency_table.loc[~contingency_table.index.isin(ignored_labels)]
    contingency_table_norm = contingency_table.div(contingency_table.sum(axis=1), axis=0)
    expected_distribution = contingency_table_norm.sum(axis=0) / contingency_table_norm.sum().sum()
    variation = contingency_table_norm.sub(expected_distribution, axis=1).abs().sum().sum() / 2
    return variation

# Generate the Jensen–Shannon distance between two distributions
def calculate_jsd_between_distributions(adata, cluster_key, label_key, ignored_labels=[]):
    contingency_table = pd.crosstab(adata.obs[label_key], adata.obs[cluster_key])
    contingency_table = contingency_table.loc[~contingency_table.index.isin(ignored_labels)]
    contingency_table_norm = contingency_table.div(contingency_table.sum(axis=1), axis=0)
    expected_distribution = contingency_table_norm.sum(axis=0) / contingency_table_norm.sum().sum()
    jsd_values = []
    for idx, row in contingency_table_norm.iterrows():
        jsd = jensenshannon(row.values, expected_distribution.values)
        jsd_values.append(jsd)
    mean_jsd = np.mean(jsd_values)
    return mean_jsd

# Calculate variation from expected for each leiden clustering vs cellTypist annotations
variation_results = []
ignored_labels = ['stripped_nuclei']
for emb in ['pca', 'scVI100', 'scPoli50']:
    for resolution in [0.1, 0.5, 1.0]:
        key = 'leiden_' + emb + '_res' + str(resolution)
        if key in adata.obs.columns:
            # Celltypist variation (biological)
            var_reed = calculate_variation_from_expected(adata, key, 'cellTypist_annotation_reed', ignored_labels=ignored_labels)
            jsd_reed = calculate_jsd_between_distributions(adata, key, 'cellTypist_annotation_reed', ignored_labels=ignored_labels)
            var_kumar = calculate_variation_from_expected(adata, key, 'cellTypist_annotation_kumar', ignored_labels=ignored_labels)
            jsd_kumar = calculate_jsd_between_distributions(adata, key, 'cellTypist_annotation_kumar', ignored_labels=ignored_labels)
            #Dataset variation (technical)
            var_dataset = calculate_variation_from_expected(adata, key, 'dataset', ignored_labels=ignored_labels)
            jsd_dataset = calculate_jsd_between_distributions(adata, key, 'dataset', ignored_labels=ignored_labels)
            #Batch variation (technical)
            var_batch = calculate_variation_from_expected(adata, key, 'batch', ignored_labels=ignored_labels)
            jsd_batch = calculate_jsd_between_distributions(adata, key, 'batch', ignored_labels=ignored_labels)
            variation_results.append({'embedding': emb, 'resolution': resolution,
                                      'variation_reed': var_reed,
                                      'variation_kumar': var_kumar,
                                      'jsd_reed': jsd_reed,
                                      'jsd_kumar': jsd_kumar,
                                      'variation_dataset': var_dataset,
                                      'jsd_dataset': jsd_dataset,
                                      'variation_batch': var_batch,
                                      'jsd_batch': jsd_batch})

variation_results = pd.DataFrame(variation_results)

#Normalize variation and JSD scores to [0,1] range for easier comparison
variation_results_normalized = variation_results.copy()
for col in ['variation_reed', 'variation_kumar', 'variation_dataset', 'variation_batch',
            'jsd_reed', 'jsd_kumar', 'jsd_dataset', 'jsd_batch']:
    for res in [0.1, 0.5, 1.0]:
        mask = variation_results_normalized['resolution'] == res
        # minimum value set to 0 which is the lowest possible variation
        min_val = 0 # variation_results_normalized.loc[mask, col].min() 
        max_val = variation_results_normalized.loc[mask, col].max()
        variation_results_normalized.loc[mask, col] = (variation_results_normalized.loc[mask, col] - min_val) / (max_val - min_val)


#Make a bar plot of the variation results
os.makedirs(save_dir + 'variation_plots/', exist_ok=True)
melted_var = variation_results.melt(id_vars=['embedding', 'resolution'])
melted_var_norm = variation_results_normalized.melt(id_vars=['embedding', 'resolution'])

for res in [0.1, 0.5, 1.0]:
    plt.figure(figsize=(12,6))
    sns.barplot(data=melted_var[melted_var['resolution'] == res], x='embedding', y='value', hue='variable')
    plt.title('Variation from Expected Distributions (Resolution ' + str(res) + ')')
    plt.ylabel('Variation Score')
    plt.savefig(save_dir + 'variation_plots/variation_from_expected_res' + str(res) + '.png')
    plt.close()
    #
    plt.figure(figsize=(12,6))
    sns.barplot(data=melted_var_norm[melted_var_norm['resolution'] == res], x='embedding', y='value', hue='variable')
    plt.title('Normalized Variation from Expected Distributions (Resolution ' + str(res) + ')')
    plt.ylabel('Normalized Variation Score')
    plt.savefig(save_dir + 'variation_plots/variation_from_expected_normalized_res' + str(res) + '.png')
    plt.close()
