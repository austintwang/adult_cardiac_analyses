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

if snakemake.params["dcc_api_key"] is not None:
    os.environ["DCC_API_KEY"] = snakemake.params["dcc_api_key"]
if snakemake.params["dcc_secret_key"] is not None:
    os.environ["DCC_SECRET_KEY"] = snakemake.params["dcc_secret_key"]

eu.connection.LOG_DIR = log_dir

conn = Connection(dcc_mode)
server = conn.dcc_url

series_data = conn.get(series)

for d in series_data["related_datasets"]:
    if "ATAC" in d["assay_term_name"]:
        experiment = d["accession"]
        break

data = conn.get(experiment)

path = None

files = data["files"]
for f in files:
    # print(f.keys()) ####
    # print(f["biological_replicates"]) ####
    if f["biological_replicates"][0] != replicate_num:
        continue
    # print(f["file_format"], f["output_type"]) ####
    if (f["file_format"] == "bam") and (f["analysis_step_version"]["aliases"][0] == "anshul-kundaje:scatac-seq-bam-filtering-step-1-0"):
        path = f["cloud_metadata"]["url"]
        break

if path is None:
    raise ValueError("href not found")

# print(path) ####

with open(path_out, "w") as f:
    f.write(path)