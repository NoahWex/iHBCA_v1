# Do principle compononent regression for LASP populations after donor pseudobulking.
# Consider the effect of Seq_data (batch), dataset, tissue_origin, and FACS_status as covariates.

#libraries
import scib

import scanpy as sc
import pandas as pd
import numpy as np
import os

import seaborn as sns
import matplotlib.pyplot as plt

#directories and data
save_dir = '/home/adr44/rds/hpc-work/hbca/data/integration/integration_comparison3/output/PCR_pseudobulk/'
ihbca_raw_data_dir = '/home/adr44/rds/rds-mammary-TaVjTYER47U/ADR44/data/hbca/iHBCA_raw_data/'

adata = sc.read_h5ad(ihbca_raw_data_dir + '/merged/preintegration_HBCA_inner_ENSEMBL_simple_lognorm_only.h5ad')
adata.obs['sequencing_date'] = adata.obs['batch']

#Make subsets for LASP and Lymphoid populations and pseudobulk by donor
adata_lasp = adata[adata.obs['level1'].isin(['Luminal adaptive secretory precursor', 
                                             'Luminal secretory', 'LP', 'AV', 
                                             'Secretory_Luminal', 'Luminal1']), :].copy()

def pseudobulk(adata, col):
    #Pseudobulk by the given column (e.g. donor)
    pb_data = []
    pb_obs = []
    for key in adata.obs[col].unique():
        temp = adata[adata.obs[col] == key, :].copy()
        #Sum the counts
        summed_counts = np.sum(temp.X, axis=0)
        pb_data.append(summed_counts)
        pb_obs.append({col: key,
                       'n_cells': temp.n_obs,
                       'dataset': temp.obs['dataset'].unique()[0],
                       'tissue_origin': temp.obs['tissue_origin'].unique()[0],
                       'FACS_status': temp.obs['FACS_status'].unique()[0],
                       'sequencing_date': temp.obs['sequencing_date'].unique()[0]})
    pb_data = np.vstack(pb_data)
    pb_obs = pd.DataFrame(pb_obs).set_index(col)
    
    pb_adata = sc.AnnData(X=pb_data, obs=pb_obs, var=adata.var.copy())
    return pb_adata

sc.pp.highly_variable_genes(adata_lasp, n_top_genes=2000)
adata_lasp_pb = pseudobulk(adata_lasp, 'patientID')

# Perform PCR using scIB
pc_regression_results1 = []
pc_regression_results2 = []
for cov in ['dataset', 'tissue_origin', 'FACS_status', 'sequencing_date']:
    print('Performing PCR for covariate: ' + cov)
    result = scib.me.pcr(adata_lasp_pb, covariate=cov, n_comps=20)
    result2 = scib.me.pcr(adata_lasp, covariate=cov, n_comps=20)
    pc_regression_results1.append({'covariate': cov, 'PCR_score': result})
    pc_regression_results2.append({'covariate': cov, 'PCR_score': result2})

# Make two simple (all grey) barplots to visualise the PCR results
os.makedirs(save_dir, exist_ok=True)
pc_regression_df1 = pd.DataFrame(pc_regression_results1)
plt.figure(figsize=(10,12))
sns.barplot(data=pc_regression_df1, x='covariate', y='PCR_score', color='grey')
plt.title('PCR scores for LASP populations (pseudobulked by donor)')
plt.ylabel('PCR score')
plt.xticks(rotation=45)
plt.ylim(0, 1)
plt.tight_layout()
plt.savefig(save_dir + 'PCR_scores_LASP_pseudobulk.png', dpi=300)
plt.savefig(save_dir + 'PCR_scores_LASP_pseudobulk.pdf', dpi=300)

pc_regression_df2 = pd.DataFrame(pc_regression_results2)
plt.figure(figsize=(10,12))
sns.barplot(data=pc_regression_df2, x='covariate', y='PCR_score', color='grey')
plt.title('PCR scores for LASP populations (single cell)')
plt.ylabel('PCR score')
plt.xticks(rotation=45)
plt.ylim(0, 1)
plt.tight_layout()
plt.savefig(save_dir + 'PCR_scores_LASP_singlecell.png', dpi=300)
plt.savefig(save_dir + 'PCR_scores_LASP_singlecell.pdf', dpi=300)





## Tried initially to do manually - but scIB has a function to do this 
# # Perform PCA regression
# # Expression matrix: samples x genes
# X = np.asarray(adata_lasp_pb.X)

# # Center
# scaler = StandardScaler(with_std=False)
# X_centered = scaler.fit_transform(X)

# n_pcs = 20
# pca = PCA(n_components=n_pcs)
# Z = pca.fit_transform(X_centered)   # PC scores (n_samples x n_pcs)

# V = pca.components_.T           

# # Design matrix
# covariates = adata_lasp_pb.obs[
#     ['dataset', 'tissue_origin', 'FACS_status', 'sequencing_date']
# ]

# # One-hot encode
# Y = pd.get_dummies(covariates, drop_first=True)

# Y = Y.values  # n_samples x n_covariates

# #select only top 5 PCs for regression
# k=5
# Z_k = Z[:, :k]
# V_k = V[:, :k]
# lambda_i = pca.explained_variance_ratio_[:k]

# pc_regression_scores = {}

# for j, cov_name in enumerate(pd.get_dummies(covariates, drop_first=True).columns):
#     r2_per_pc = np.array([
#         r2_score(Y[:, j], Z_k[:, i])
#         for i in range(k)
#     ])
    
#     pc_regression_scores[cov_name] = np.sum(lambda_i * r2_per_pc)
