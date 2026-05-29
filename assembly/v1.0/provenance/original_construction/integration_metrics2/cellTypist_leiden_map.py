#Generate cellTypist labels for Kumar and Reed labels.


#libraries
import scanpy as sc
import celltypist
from celltypist import models
import pandas as pd
import numpy as np

import os

#directories
ihbca_raw_data_dir = '/home/adr44/rds/rds-mammary-TaVjTYER47U/ADR44/data/hbca/iHBCA_raw_data/'
save_dir = '/home/adr44/rds/hpc-work/hbca/data/integration/integration_comparison2/output/celltypist/'
os.makedirs(save_dir, exist_ok=True)


#Load adata
print('Loading full dataset')
adata = sc.read_h5ad(ihbca_raw_data_dir + '/merged/preintegration_HBCA_inner_ENSEMBL_simple_lognorm_only.h5ad')

#Load custom cellTypist models
model_reed = models.Model.load('/home/adr44/rds/hpc-work/hbca/DATA_FINAL/celltypist/model_from_reed.pkl')
model_kumar = models.Model.load('/home/adr44/rds/hpc-work/hbca/DATA_FINAL/celltypist/model_from_kumar.pkl')

# make predictions
print('Predicting labels...')
pred_reed = celltypist.annotate(adata, model=model_reed, majority_voting=False)
pred_kumar = celltypist.annotate(adata, model=model_kumar, majority_voting=False)

# #Only works if majority_voting=True
# adata.obs['celltypist_reed'] = pred_reed.predicted_labels.values
# adata.obs['celltypist_kumar'] = pred_kumar.predicted_labels.values

#Save the cellTypist annotations
pred_reed.predicted_labels.to_csv(
    save_dir + "celltypist_reed_predicted_labels.csv"
)

pred_reed.decision_matrix.to_csv(
    save_dir + "celltypist_reed_decision_matrix.csv"
)

pred_reed.probability_matrix.to_csv(
    save_dir + "celltypist_reed_probability_matrix.csv"
)

pred_kumar.predicted_labels.to_csv(
    save_dir + "celltypist_kumar_predicted_labels.csv"
)

pred_kumar.decision_matrix.to_csv(
    save_dir + "celltypist_kumar_decision_matrix.csv"
)

pred_kumar.probability_matrix.to_csv(
    save_dir + "celltypist_kumar_probability_matrix.csv"
)




