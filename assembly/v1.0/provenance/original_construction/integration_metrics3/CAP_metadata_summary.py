# Make some basic summary plots and stats for the CAP metadata file

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os


#data 
cellxgene_dir = '/home/adr44/rds/rds-mammary-TaVjTYER47U/ADR44/data/hbca/DATA_FINAL/cellxgene/'
metadata = pd.read_csv(cellxgene_dir + 'integrated_HBCA_cellxgene_CAP_metadata_NEW.csv', 
                       dtype=str, #errors keep occuring otherwise
                       index_col=0)


#Summar of the percentage of NA or unknown values per column
na_percentages = (metadata.isnull() | metadata.isna() | (metadata == 'unknown') | (metadata == 'Unknown')).sum() / len(metadata) * 100
print(na_percentages)

#Make some basic bar plots for Tier 1 and Tier 2 metadata columns
tier2_cols = ['age_exact', 'age_binary', 'menarche_age', 'parity_exact', 'parity_binary', 'gravidity', 
              'breast_feeding', 'time_since_last_pregnancy', 'age_at_first_pregnancy', 'menopause_status',
              'hormone_replacement_therapy', 'contraceptive_history', 'body_mass_index', 'smoking_history',
              'surgery_type', 'BRCA1_status', 'BRCA2_status', 'other_germline_mutation', 'sample_type']

metadata1 = metadata.loc[:, metadata.columns.isin(tier2_cols) == False]
metadata2 = metadata.loc[:, metadata.columns.isin(tier2_cols) == True]

save_dir = '/home/adr44/rds/hpc-work/hbca/data/integration/integration_comparison3/output/CAP_metadata_summary/'
os.makedirs(save_dir, exist_ok=True)

for df, tier in zip([metadata1, metadata2], ['Tier1', 'Tier2']):
    na_percentages = (df.isnull() | df.isna() | (df == 'unknown') | (df == 'Unknown')).sum() / len(df) * 100
    df2 = 100 - na_percentages
    plt.clf()
    plt.figure(figsize=(20, 15))
    sns.barplot(x=df2.index, y=df2.values, color='grey')
    plt.xticks(rotation=90)
    plt.ylabel('Percentage of filled values (%)')
    plt.ylim(0, 100)
    plt.title(f'Percentage of filled values per column - {tier} metadata')
    plt.tight_layout()
    plt.savefig(save_dir + f'CAP_metadata_{tier}_filled_values_percentage_barplot.png', dpi=300)
    df2.to_csv(save_dir + f'CAP_metadata_{tier}_filled_values_percentage.csv')
    total = df2.mean()
    print(f'Overall percentage of {tier} metadata filled values: {total:.2f}%')

# Overall percentage of Tier 1 and Tier 2 metadata filled values
# tier1_total = 96.64%
# tier2_total = 42.04%

