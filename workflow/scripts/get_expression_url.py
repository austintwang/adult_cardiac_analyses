import os
import json
import statistics
from urllib.parse import urljoin
import encode_utils as eu
from encode_utils.connection import Connection

path_out, = snakemake.output
dcc_mode = snakemake.config["dcc_mode"]
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
    if "RNA" in d["assay_term_name"]:
        experiment = d["accession"]
        break

data = conn.get(experiment)

path = None

files = data["files"]
for f in files:
    if f["biological_replicates"][0] != replicate_num:
        continue
    if (f["file_format"] == "tar") and (f["output_type"] == "sparse gene count matrix of unique reads"):
        path = f["cloud_metadata"]["url"]
        break

if path is None:
    raise ValueError("href not found for expression tarball")

# print(path) ####

with open(path_out, "w") as f:
    f.write(path)