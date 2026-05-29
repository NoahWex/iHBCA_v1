#Validate CAP object for iHBCA

from cap_upload_validator import UploadValidator

h5ad_path = "/home/adr44/rds/rds-mammary-TaVjTYER47U/ADR44/data/hbca/DATA_FINAL/cellxgene/iHBCA_CAP_object.h5ad"

uv = UploadValidator(h5ad_path)
uv.validate()