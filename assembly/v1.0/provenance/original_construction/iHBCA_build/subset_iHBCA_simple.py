#Code to subset the iHBCA data and save

# full merging keeps failing likely due to the size of the object, 
# so try this more simple (but less elegant) approach


# conda env py-scvi

#libraries
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
import anndata as ad

import glob
import os
from collections import Counter
import matplotlib.pyplot as plt


# Load the full iHBCA data
ihbca_raw_data_dir = '/home/adr44/rds/rds-mammary-TaVjTYER47U/ADR44/data/hbca/iHBCA_raw_data/'
adata_iHBCA = sc.read_h5ad(ihbca_raw_data_dir + '/merged/preintegration_HBCA_inner_ENSEMBL_simple.h5ad')


#Make a subsetted object and save (400k cells)
np.random.seed(0)
adata_sub = adata_iHBCA[np.random.choice(adata_iHBCA.obs.index, 400000, replace=False), :]
#save
adata_sub.write(ihbca_raw_data_dir + '/merged/preintegration_HBCA_inner_ENSEMBL_simple_sub400k.h5ad')



#Make a subsetted object and save (100k cells)
np.random.seed(0)
adata_sub = adata_iHBCA[np.random.choice(adata_iHBCA.obs.index, 100000, replace=False), :]
#save
adata_sub.write(ihbca_raw_data_dir + '/merged/preintegration_HBCA_inner_ENSEMBL_simple_sub100k.h5ad')


