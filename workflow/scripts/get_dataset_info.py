import os
import json
import encode_utils as eu
from encode_utils.connection import Connection

path_out, = snakemake.output
dcc_mode = snakemake.config["dcc_mode"]
sample = snakemake.wildcards["sample"]
series = snakemake.params["series"]
replicate_num = int(snakemake.params["replicate"])
log_dir, = snakemake.log

os.environ["DCC_API_KEY"] = snakemake.params["dcc_api_key"]
os.environ["DCC_SECRET_KEY"] = snakemake.params["dcc_secret_key"]

eu.connection.LOG_DIR = log_dir

conn = Connection(dcc_mode)
server = conn.dcc_url

series_data = conn.get(series)

for d in series_data["related_datasets"]:
    if "ATAC" in d["assay_term_name"]:
        experiment_atac = d["accession"]
        break

for d in series_data["related_datasets"]:
    if "RNA" in d["assay_term_name"]:
        experiment_rna = d["accession"]
        break

data = {sample: {
    "rna": experiment_rna,
    "atac": experiment_atac
}}

with open(path_out, "w") as f:
    json.dump(data, f, indent=4)