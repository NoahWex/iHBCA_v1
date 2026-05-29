# Update iHBCA metadata to match HCA guidelines

# conda env py-scanpy

#libraries
import scanpy as sc
import pandas as pd
import numpy as np
import os



#load data
ihbca_raw_data_dir = '/home/adr44/rds/rds-mammary-TaVjTYER47U/ADR44/data/hbca/iHBCA_raw_data/'
cellxgene_dir = '/home/adr44/rds/rds-mammary-TaVjTYER47U/ADR44/data/hbca/DATA_FINAL/cellxgene/'
# adata = sc.read_h5ad(cellxgene_dir + 'integrated_HBCA_cellxgene.h5ad') #OLD use the updated object now
adata = sc.read_h5ad(ihbca_raw_data_dir + '/merged/preintegration_HBCA_inner_ENSEMBL_simple_lognorm_only.h5ad')

# extra_metadata = pd.read_csv(cellxgene_dir + '../extra_metadata/metadata_iHBCA_extra_v2.csv', index_col=0)
kumar_extrametadata = pd.read_csv('/home/adr44/rds/rds-mammary-TaVjTYER47U/ADR44/data/hbca/iHBCA_raw_data/kumar_data/formatted/kumar_full_phenodata.csv')
nee_extrametadata = pd.read_csv('/home/adr44/rds/rds-mammary-TaVjTYER47U/ADR44/data/hbca/iHBCA_raw_data/nee_data/formatted/nee_full_phenodata.csv')
murrow_extrametadata = pd.read_csv('/home/adr44/rds/rds-mammary-TaVjTYER47U/ADR44/data/hbca/iHBCA_raw_data/murrow_data/formatted/murrow_full_phenodata.csv')

## Formatting

# #HCA Tier 1 metadata
# sample_id donor_id protocol_url institute sample_collection_site sample_collection_relative_time_point 
# library_id library_id_repository author_batch_notes organism_ontology_term_id manner_of_death sample_source 
# sex_ontology_term_id sample_collection_method tissue_type sampled_site_condition tissue_ontology_term_id 
# tissue_free_text sample_preservation_method suspension_type cell_enrichment cell_viability_percentage 
# cell_number_loaded sample_collection_year assay_ontology_term_id library_preparation_batch 
# library_sequencing_run sequenced_fragment sequencing_platform is_primary_data reference_genome 
# gene_annotation_version alignment_software intron_inclusion disease_ontology_term_id 
# self_reported_ethnicity_ontology_term_id development_stage_ontology_term_id

# #HBCA Tier 2 metadata
# age_exact age_binary menarche_age parity_exact parity_binary gravidity breast_feeding time_since_last_pregnancy 
# age_at_first_pregnancy menopause_status hormone_replacement_therapy contraceptive_history body_mass_index 
# smoking_history surgery_type BRCA1_status BRCA2_status other_germline_mutation sample_type


# #iHBCA v1 metadata
# donor_id dataset doi sample_type tissue_origin assay_ontology_term_id batch FACS_status ethnicity donor_age
# development_stage_ontology_term_id parity risk_status level2 cell_type_ontology_term_id original_celltype 
# map_celltype_map2gray map_celltype_map2murrow map_celltype_map2nee map_celltype_map2pal map_celltype_map2twigger
# map_celltype_map2kumar map_celltype_map2reed organism_ontology_term_id tissue_ontology_term_id
# disease_ontology_term_id suspension_type self_reported_ethnicity_ontology_term_id sex_ontology_term_id


#Some metadata formatting and adjustments

#Kumar
kumar_extrametadata['donor_age'] = kumar_extrametadata['age'].map({'O': 'old', 'Y': 'young'})
#kumar_extrametadata['cellID'] = kumar_extrametadata['Unnamed: 0']

#Nee
nee_extrametadata['celltype_enriched'] = nee_extrametadata['orig.ident'].str.split('_').apply(lambda x: x[-1] if len(x) >= 3 else 'na').replace({'BRCA': 'na'})
nee_extrametadata['celltype_enriched'] = nee_extrametadata['celltype_enriched'].map({'na': 'na', 'Epi': 'CL:0002327+', 'Stroma': 'CL:0002327-'})
nee_extrametadata['cellID'] = nee_extrametadata['Unnamed: 0']

#Murrow
murrow_extrametadata['celltype_enriched'] = murrow_extrametadata['Sort'].map({'Live_singlet': 'DAPI+', 
                                                                              'Epithelial': 'CL:0002327+',
                                                                              'Luminal': 'CL:0002326+',
                                                                              'Basal': 'CL:0002324+'})
murrow_extrametadata['assay_ontology_term_id'] = murrow_extrametadata['Batch'].map({'Batch_1': 'EFO:0009899', 'Batch_2': 'EFO:0009899', 'Batch_3': 'EFO:0009922', 'Batch_4': 'EFO:0009922'})
murrow_extrametadata['sequencing_platform'] = murrow_extrametadata['Batch'].map({'Batch_1': 'EFO:0008563', 'Batch_2': 'EFO:0008563', 'Batch_3': 'EFO:0008637', 'Batch_4': 'EFO:0008637'})
murrow_extrametadata['cellID'] = murrow_extrametadata['Unnamed: 0']


#iHBCA
donor_mapping = pd.read_csv('/home/adr44/rds/rds-mammary-TaVjTYER47U/ADR44/data/hbca/DATA_FINAL/extra_metadata/HBCA_donor_id_anon_map.csv', header=None)
adata.obs['dataset'] = adata.obs.dataset.str.lower().str.capitalize()
adata.obs['donor_id'] = adata.obs.patientID.copy().astype(str)
adata.obs['donor_id'][adata.obs['dataset']=='Reed'] = adata.obs['patientID'][adata.obs['dataset']=='Reed'].map(dict(zip(donor_mapping[0], donor_mapping[1])))
adata.obs['cellID'] = adata.obs.index
adata.obs['cellID_first_part'] = adata.obs['cellID'].str.rsplit('_', n=1).str[0]
adata.obs['cellID_second_part'] = adata.obs['cellID'].str.rsplit('_', n=1).str[1]
adata.obs['barcode'] = adata.obs['cellID_second_part']
adata.obs['barcode'][adata.obs['dataset']=='Reed'] = adata.obs['cellID'].str.split('-', n=1).str[0][adata.obs['dataset']=='Reed']
adata.obs['sampleID'] = adata.obs['cellID_first_part']
adata.obs['sampleID'][adata.obs['dataset']=='Reed'] = adata.obs['cellID'].str.split('-', n=1).str[1][adata.obs['dataset']=='Reed']
adata.obs['cellID_new'] = adata.obs['dataset'].astype(str) + '__' + adata.obs['sampleID'] + '__' + adata.obs['barcode']
adata.obs['cell_enrichment'] = 'na'
adata.obs['cell_enrichment'][adata.obs['dataset']=='Reed' & adata.obs['FACS_status']=='cell_type_sorted'] = 'CL:4033057+' #LASP cells
adata.obs['cell_enrichment'][adata.obs['dataset']=='Reed' & adata.obs['FACS_status']=='live_sorted'] = 'DAPI+' #DAPI/live cells ###CHECK THIS
adata.obs['cell_enrichment'][adata.obs['dataset']=='Reed' & adata.obs['FACS_status']=='not_sorted'] = 'na' 
adata.obs['cell_enrichment'][adata.obs['dataset']=='Kumar'] = 'DAPI+' #DAPI/live cells
adata.obs['cell_enrichment'][adata.obs['dataset']=='Nee'] = adata.obs.cellID[adata.obs['dataset']=='Nee'].map(dict(zip(nee_extrametadata['cellID'], nee_extrametadata['celltype_enriched'])))
adata.obs['cell_enrichment'][adata.obs['dataset']=='Pal'] = adata.obs.sampleID[adata.obs['dataset']=='Pal'].str.split('_', n=1).str[1].map({'mix': 'DAPI+', 'epi': 'CL:0002327+'})
adata.obs['cell_enrichment'][adata.obs['dataset']=='Twigger'] = 'na'
adata.obs['cell_enrichment'][adata.obs['dataset']=='Murrow'] = adata.obs.cellID[adata.obs['dataset']=='Murrow'].map(dict(zip(murrow_extrametadata['cellID'], 
                                                                                                                             murrow_extrametadata['celltype_enriched'])))
adata.obs['cell_enrichment'][adata.obs['dataset']=='Gray'] = 'na'
adata.obs['assay_ontology_term_id'] = adata.obs.dataset.map({'Reed': 'EFO:0009922', 'Kumar': 'na', 'Nee': 'EFO:0009899', 'Pal': 'EFO:0009922', 
                                                             'Twigger': 'EFO:0009922', 'Murrow': 'na', 'Gray': 'EFO:0009922'})
adata.obs['assay_ontology_term_id'][adata.obs['dataset']=='Kumar'] = adata.obs.cellID[adata.obs['dataset']=='Kumar'].map(dict(zip(kumar_extrametadata['cellID'], 
                                                                                                                                  kumar_extrametadata['assay_ontology_term_id'])))
adata.obs['assay_ontology_term_id'][adata.obs['dataset']=='Murrow'] = adata.obs.cellID[adata.obs['dataset']=='Murrow'].map(dict(zip(murrow_extrametadata['cellID'], 
                                                                                                                                    murrow_extrametadata['assay_ontology_term_id'])))
adata.obs['sequencing_platform'] = adata.obs.dataset.map({'Reed': 'EFO:0008637', 'Kumar': 'EFO:0008637', 'Nee': 'EFO:0008637', # both HiSeq4000 and NovaSeq6000 quoted but not sure which for what cells
                                                          'Pal': 'EFO:0009173', 'Twigger': 'EFO:0008637', 'Murrow': 'na', 'Gray': 'EFO:0008567'})
adata.obs['sequencing_platform'][adata.obs['dataset']=='Murrow'] = adata.obs.cellID[adata.obs['dataset']=='Murrow'].map(dict(zip(murrow_extrametadata['cellID'], murrow_extrametadata['sequencing_platform'])))
adata.obs['reference_genome'] = adata.obs.dataset.map({'Reed': 'GRCh38', 'Kumar': 'GRCh38', 'Nee': 'GRCh38', 'Pal': 'GRCh38', 'Twigger': 'GRCh37', 'Murrow': 'GRCh37', 'Gray': 'GRCh38'})
adata.obs['alignment_software'] = adata.obs.dataset.map({'Reed': 'cell ranger 6.0.2', 'Kumar': 'cell ranger 3.1.0', 'Nee': 'cell ranger 3.1.0', 'Pal': 'cell ranger 3.02', 
                                                         'Twigger': 'cell ranger 3.0.2', 'Murrow': 'cell ranger 3.0.2', 'Gray': 'cell ranger 3.0.2'})

# adata.obs['age_exact'] = adata.obs['donor_age']
# adata.obs['age_binary'] = adata.obs['donor_age']
# adata.obs['age_binary'] = (adata.obs['donor_age'] >= 50).map({True: 'old', False: 'young'})
# adata.obs['age_binary'][adata.obs['dataset']=='Kumar'] = adata.obs.cellID[adata.obs['dataset']=='Kumar'].map(dict(zip(kumar_extrametadata['cellID'], kumar_extrametadata['donor_age'])))
adata.obs['parity_exact'] = adata.obs['parity'].astype(str).map({'0': '0', '0.0': '0', 
                                                                 '1': '1', '1.0': '1', 
                                                                 '2': '2', '2.0': '2', 
                                                                 '3': '3', '3.0': '3', 
                                                                 '4': '4', '4.0': '4', 
                                                                 'Nulliparous': '0', 
                                                                 'Nulliarous': '0', 
                                                                 'Parous': 'unknown', 
                                                                 'Unknown': 'unknown', 'unknown': 'unknown'}) #fixes some typos and mixed formatting
adata.obs['parity_binary'] = adata.obs['parity'].astype(str).map({'0': 'nulliparous', '0.0': 'nulliparous', 
                                                                  '1': 'parous', '1.0': 'parous',
                                                                  '2': 'parous', '2.0': 'parous',
                                                                  '3': 'parous', '3.0': 'parous',
                                                                  '4': 'parous', '4.0': 'parous',
                                                                  'Nulliparous': 'nulliparous', 
                                                                  'Nulliarous': 'nulliparous',
                                                                  'Parous': 'parous', 'Unknown': 'unknown',
                                                                  'unknown': 'unknown'}) #fixes some typos and mixed formatting
adata.obs['breast_feeding'] = 'unknown' # missing for most samples
adata.obs['breast_feeding'][adata.obs['parity_binary']=='nulliparous'] = 'na'
adata.obs['surgery_type'] = adata.obs['risk_status'].map({'AR': 'reduction mammoplasty', 
                                                          'HR-BR1': 'prophylactic mastectomy',
                                                          'HR-BR2': 'prophylactic mastectomy',
                                                          'HR-Unk': 'prophylactic mastectomy',
                                                          'HR-cBR1': 'contralateral mastectomy',
                                                          'HR-cUnk': 'contralateral mastectomy',
                                                          'HR-cBR2': 'contralateral mastectomy',
                                                          'HR-RAD': 'prophylactic mastectomy'})
adata.obs['BRCA1_status'] = adata.obs['risk_status'].map({'AR': 'WT',
                                                          'HR-BR1': 'BRCA1',
                                                          'HR-BR2': 'WT',
                                                          'HR-Unk': 'unknown',
                                                          'HR-cBR1': 'BRCA1',
                                                          'HR-cUnk': 'unknown',
                                                          'HR-cBR2': 'WT',
                                                          'HR-RAD': 'unknown'})
adata.obs['BRCA1_status'][(adata.obs['dataset'].astype(str)=='Reed') & (adata.obs['risk_status'].astype(str)=='HR-Unk')] = 'WT'
adata.obs['BRCA2_status'] = adata.obs['risk_status'].map({'AR': 'WT',
                                                          'HR-BR1': 'WT',
                                                          'HR-BR2': 'BRCA2',
                                                          'HR-Unk': 'unknown',
                                                          'HR-cBR1': 'WT',
                                                          'HR-cUnk': 'unknown',
                                                          'HR-cBR2': 'BRCA2',
                                                          'HR-RAD': 'unknown'})
adata.obs['BRCA2_status'][(adata.obs['dataset'].astype(str)=='Reed') & (adata.obs['risk_status'].astype(str)=='HR-Unk')] = 'WT'
adata.obs['other_germline_mutation'] = adata.obs['risk_status'].map({'AR': 'na',
                                                                     'HR-BR1': 'na',
                                                                     'HR-BR2': 'na',
                                                                     'HR-Unk': 'na',
                                                                     'HR-cBR1': 'na',
                                                                     'HR-cUnk': 'na',
                                                                     'HR-cBR2': 'na',
                                                                     'HR-RAD': 'RAD51C'})
# sample_type
# mixed          1193992
# Supernatant     406194
# Organoid        233709
# Organoid LP     163380
# Unknown         124790
adata.obs['sample_type_new'] = adata.obs.sample_type.map({'mixed': 'mixed', 
                                                          'Supernatant': 'stromal enriched',
                                                          'Organoid': 'epithelial enriched',
                                                          'Organoid LP': 'epithelial enriched',
                                                          'Unknown': 'mixed'})

# Metadata mapping dictionaries

protocol_url_dict = {'Reed': 'https://www.nature.com/articles/s41588-024-01688-9',
                     'Kumar': 'https://www.nature.com/articles/s41586-023-06252-9',
                     'Nee': 'https://www.nature.com/articles/s41588-023-01298-x',
                     'Pal': 'https://www.embopress.org/doi/full/10.15252/embj.2020107333',
                     'Twigger': 'https://www.nature.com/articles/s41467-021-27895-0',
                     'Murrow': 'https://www.sciencedirect.com/science/article/pii/S2405471222002757?via%3Dihub',
                     'Gray': 'https://www.sciencedirect.com/science/article/pii/S1534580722003318?via%3Dihub'}

#from
# ['Asian_or_Asian_British-Any_other_Asian_backgr..., 'Asian_or_Asian_British-Bangladeshi',
#                          'Black_or_Black_British-Any_other_Black_backgr..., 'Other_ethnic_group-Any_other_ethnic_group',
#                          'Other_ethnic_group-Not_stated', 'White-Any_other_background', 'White-British',
#                          'White-Irish']
#to Asian, British, Bangladeshi, Irish, Unknown
ethnicity_dict = {'Asian_or_Asian_British-Any_other_Asian_background': 'Asian',
                  'Asian_or_Asian_British-Bangladeshi': 'Bangladeshi',
                  'Black_or_Black_British-Any_other_Black_background': 'Unknown',
                  'Other_ethnic_group-Any_other_ethnic_group': 'Unknown',
                  'Other_ethnic_group-Not_stated': 'Unknown',
                  'White-Any_other_background': 'Unknown',
                  'White-British': 'British',
                  'White-Irish': 'Irish'}


# kumar_extrametadata['HBCA_donor_id'] = 'Kumar_' + kumar_extrametadata['donor_id'].astype(str)
kumar_institute_dict = dict(zip(kumar_extrametadata['donor_id'].astype(str), kumar_extrametadata['sample_source']))
other_institute_dict = {'Reed': 'University of Cambridge', 'Nee': 'University of California, Irvine', 'Pal': 'Walter and Eliza Hall Institute of Medical Research',
                        'Twigger': 'Helmholtz Zentrum München', 'Murrow': 'University of California, San Francisco', 'Gray': 'Broad Institute'}

# Create new metadata df
adata_obs = adata.obs[['donor_id', 'dataset']].copy()
adata_obs.columns = ['old_donor_id', 'old_dataset']

#Add HCA Tier 1 metadata
adata_obs['sample_id'] = adata.obs.donor_id.astype(str) + '_' + adata.obs.batch.astype(str) + '_' + adata.obs.sample_type.astype(str).str.replace(' ', '_')
adata_obs['donor_id'] = adata.dataset + '_' + adata.obs.donor_id
adata_obs['protocol_url'] = adata.obs.dataset.map(protocol_url_dict)
adata_obs['institute'] = adata.obs.old_donor_id.map(kumar_institute_dict)
adata_obs['institute'] = adata_obs['institute'].fillna(adata.obs['dataset'].map(other_institute_dict))
adata_obs['sample_collection_site'] = 'na' #not sure what to do on this one
adata_obs['sample_collection_relative_time_point'] = 'na' #not sure what to do on this one
adata_obs['library_id'] = adata.obs['sampleID']
adata_obs['library_id_repository'] = 'na' #not sure what to do on this one
adata_obs['author_batch_notes'] = 'na' #not sure what to do on this one 
adata_obs['organism_ontology_term_id'] = 'NCBITaxon:9606'
adata_obs['manner_of_death'] = 'not applicable' 
adata_obs['sample_source'] = 'surgical donor'
adata_obs['sex_ontology_term_id'] = 'PATO:0000383' #female
adata_obs['sample_collection_method'] = 'surgical resection'
adata_obs['tissue_type'] = 'tissue'
adata_obs['sampled_site_condition'] = 'healthy'
adata_obs['tissue_ontology_term_id'] = 'UBERON_0000310' #breast
adata_obs['tissue_free_text'] = 'breast'
adata_obs['sample_preservation_method'] = adata.obs['tissue_origin']
adata_obs['suspension_type'] = 'cell'
adata_obs['cell_enrichment'] = adata.obs['cell_enrichment']
adata_obs['cell_viability_percentage'] = 'na' # hard to collect
adata_obs['cell_number_loaded'] = 'na' # hard to collect
adata_obs['sample_collection_year'] = 'na' # hard to collect
adata_obs['assay_ontology_term_id'] = adata.obs['assay_ontology_term_id']
adata_obs['library_preparation_batch'] = adata.obs['dataset'].astype(str) + '_' + adata.obs['batch'].astype(str)
adata_obs['library_sequencing_run'] = 'na' # hard to collect
adata_obs['sequenced_fragment'] = '3 prime tag'
adata_obs['sequencing_platform'] = adata.obs['sequencing_platform']
adata_obs['is_primary_data'] = 'false'
adata_obs['reference_genome'] = adata.obs['reference_genome']
adata_obs['gene_annotation_version'] = 'na' # hard to collect
adata_obs['alignment_software'] = adata.obs['alignment_software']
adata_obs['intron_inclusion'] = 'no'
adata_obs['disease_ontology_term_id'] = 'PATO:0000461' #normal or healthy
ethnicity_development_stage_df = pd.read_csv(cellxgene_dir + 'donorID_to_ethnicity_plus_developmental_stage_mapping.csv', index_col=0)
ethnicity_development_stage_df['donor_id'] = ethnicity_development_stage_df.donor_id.str.split('_', n=1).str[1]
adata_obs['self_reported_ethnicity_ontology_term_id'] = adata.obs.donor_id.map(dict(zip(ethnicity_development_stage_df['donor_id'], ethnicity_development_stage_df['self_reported_ethnicity_ontology_term_id'])))
adata_obs['development_stage_ontology_term_id'] = adata.obs.donor_id.map(dict(zip(ethnicity_development_stage_df['donor_id'], ethnicity_development_stage_df['development_stage_ontology_term_id'])))

#Add HBCA Tier 2 metadata
tier2_cols = ['age_exact', 'age_binary', 'menarche_age', 'parity_exact', 'parity_binary', 'gravidity', 
              'breast_feeding', 'time_since_last_pregnancy', 'age_at_first_pregnancy', 'menopause_status',
              'hormone_replacement_therapy', 'contraceptive_history', 'body_mass_index', 'smoking_history',
              'surgery_type', 'BRCA1_status', 'BRCA2_status', 'other_germline_mutation', 'sample_type']
adata_obs['age_exact'] = pd.to_numeric(
    adata.obs["age"].astype(str),
    errors="coerce"
)
adata_obs["age_binary"] = np.where(
    adata_obs["age_exact"].notna(),
    np.where(adata_obs["age_exact"] >= 50, "O", "Y"),
    np.nan
)
adata_obs['menarche_age'] = 'unknown' # missing for most samples
adata_obs['parity_exact'] = adata.obs['parity_exact'].astype(str).replace({"nan": None})
adata_obs['parity_binary'] = adata.obs['parity_binary']
adata_obs['gravidity'] = 'unknown' # missing for most samples
adata_obs['breast_feeding'] = adata.obs['breast_feeding'] #missing for most samples
adata_obs['time_since_last_pregnancy'] = 'unknown' # missing for most samples
adata_obs['age_at_first_pregnancy'] = 'unknown' # missing for most samples
adata_obs['menopause_status'] = 'unknown' # missing for most samples
adata_obs['hormone_replacement_therapy'] = 'unknown' # missing for most samples
adata_obs['contraceptive_history'] = 'unknown' # missing for most samples
adata_obs['body_mass_index'] = 'unknown' # missing for most samples
adata_obs['smoking_history'] = 'unknown' # missing for most samples
adata_obs['surgery_type'] = adata.obs['surgery_type']
adata_obs['BRCA1_status'] = adata.obs['BRCA1_status']
adata_obs['BRCA2_status'] = adata.obs['BRCA2_status']
adata_obs['other_germline_mutation'] = adata.obs['other_germline_mutation']
adata_obs['sample_type'] = adata.obs['sample_type_new']

#Save updated metadata
adata_obs.to_csv(cellxgene_dir + 'integrated_HBCA_cellxgene_CAP_metadata_NEW.csv')
# adata_obs2 = pd.read_csv(cellxgene_dir + 'integrated_HBCA_cellxgene_CAP_metadata.csv', index_col=0)


# #Scraps
# df = adata_obs2[['donor_id', 'self_reported_ethnicity_ontology_term_id', 'development_stage_ontology_term_id']].drop_duplicates(subset=['donor_id'])
# df.to_csv(cellxgene_dir + 'donorID_to_ethnicity_plus_developmental_stage_mapping.csv')
