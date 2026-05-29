#Code to summarise and plot the metrics for iHBCA integrations


# conda env py-scvi 
# (py-scanpy env has a bug with seaborn version affecting the heatmaps)

#libraries
import scanpy as sc
import scib
import pandas as pd

import os
import sys
import matplotlib.pyplot as plt
import seaborn as sns




##### Load integration metrics #####

subset = str(sys.argv[1])
# subset = 'sub100k'
# subset = '40kPerDataset'


## scPoli (standard)
start = 0
for latent_dim in [20, 50, 100, 200]:
    if (latent_dim == 20) | (latent_dim == 50):
        hidden_layer_size = 128
    else:
        hidden_layer_size = 512
    scpoli_metrics = pd.read_csv('/home/adr44/rds/hpc-work/hbca/data/integration/integration_comparison/output/' + subset + '/scPoli_metrics/n_dims_' + str(latent_dim) + '_integration_metrics.csv')
    scpoli_metrics.columns = ['Integration_metric', 'scPoli_hid_' + str(hidden_layer_size) + '_lat_' + str(latent_dim)]
    if start == 0:
        integration_metrics = scpoli_metrics
        start = 1
    else:
        integration_metrics = pd.merge(integration_metrics, scpoli_metrics, on='Integration_metric', how='left')



## scPoli (tests)
for test_type in ['output_moreHidden', 'output_moreDropout', 'output_matchHidden', 'output_intDataset', 'output_intDataset2', 'output_celltype', 'output_intDataset_celltype']: #
    for latent_dim in [50]:
        if (latent_dim == 20) | (latent_dim == 50):
            hidden_layer_size = 128
        else:
            hidden_layer_size = 512
        if test_type == 'output_moreHidden':
            hidden_layer_size = '2048_' + str(hidden_layer_size)
        scpoli_metrics = pd.read_csv('/home/adr44/rds/hpc-work/hbca/data/integration/integration_comparison/' + test_type + '/' + subset + '/scPoli_metrics/n_dims_' + str(latent_dim) + '_integration_metrics.csv')  
        scpoli_metrics.columns = ['Integration_metric', 'scPoli_hid_' + str(hidden_layer_size) + '_lat_' + str(latent_dim) + '_' + test_type]  
        integration_metrics = pd.merge(integration_metrics, scpoli_metrics, on='Integration_metric', how='left')



## Process all the metrics

#Save 
os.makedirs('/home/adr44/rds/hpc-work/hbca/data/integration/integration_comparison/output_tests/' + subset + '/', exist_ok=True)
integration_metrics.to_csv('/home/adr44/rds/hpc-work/hbca/data/integration/integration_comparison/output_tests/' + subset + '/integration_metrics_scPoli.csv', index=False)
# integration_metrics = pd.read_csv('/home/adr44/rds/hpc-work/hbca/data/integration/integration_comparison/output_tests/' + subset + '/integration_metrics_scPoli.csv')

#drop NAs
integration_metrics = integration_metrics.dropna()
integration_metrics.set_index('Integration_metric', inplace=True)

#add the averages
integration_metrics.loc['Average',:] = integration_metrics.mean(axis=0)
integration_metrics_ord = integration_metrics.sort_values(by='Average', axis=1, ascending=False)

#melt
integration_metrics_melt = pd.melt(integration_metrics, var_name='Integration_method', value_name='Value') #, id_vars=['Integration_metric']



##### Plot integration metrics #####

#Make a heatmap with metrics on the x axis and integration methods on the y axis
plt.clf()
plt.figure(figsize=(20, 20))
sns.heatmap(integration_metrics.astype(float).T, annot=True, vmin=0, vmax=1) 
plt.tight_layout()
plt.savefig('/home/adr44/rds/hpc-work/hbca/data/integration/integration_comparison/output_tests/' + subset + '/integration_metrics_heatmap.png')

plt.clf()
plt.figure(figsize=(20, 20))
sns.heatmap(integration_metrics_ord.astype(float).T, annot=True, vmin=0, vmax=1) 
plt.tight_layout()
plt.savefig('/home/adr44/rds/hpc-work/hbca/data/integration/integration_comparison/output_tests/' + subset + '/integration_metrics_heatmap_ord.png')

#scaled heatmaps per row
integration_metrics_scaled = integration_metrics.copy().to_numpy()
integration_metrics_scaled = (integration_metrics_scaled) / (integration_metrics_scaled.max(axis=1))[:, None] # - integration_metrics_scaled.min(axis=1)[:, None]
integration_metrics_scaled_df = pd.DataFrame(integration_metrics_scaled, columns=integration_metrics.columns, index=integration_metrics.index)
integration_metrics_scaled_df.loc['Average_2',:] = integration_metrics_scaled_df.mean(axis=0)
integration_metrics_scaled_df_ord = integration_metrics_scaled_df.sort_values(by='Average', axis=1, ascending=False)

plt.clf()
plt.figure(figsize=(20, 20))
sns.heatmap(integration_metrics_scaled_df.astype(float).T, annot=True, vmin=0, vmax=1)
plt.tight_layout()
plt.savefig('/home/adr44/rds/hpc-work/hbca/data/integration/integration_comparison/output_tests/' + subset + '/integration_metrics_heatmap_scaled.png')

plt.clf()
plt.figure(figsize=(20, 20))
sns.heatmap(integration_metrics_scaled_df_ord.astype(float).T, annot=True, vmin=0, vmax=1)
plt.tight_layout()
plt.savefig('/home/adr44/rds/hpc-work/hbca/data/integration/integration_comparison/output_tests/' + subset + '/integration_metrics_heatmap_scaled_ord.png')
