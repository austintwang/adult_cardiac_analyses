import os
import json
import statistics
from urllib.parse import urljoin
from urllib.request import urlopen
import gzip
import encode_utils as eu
from encode_utils.connection import Connection

def main(dataset_path):

    dcc_mode = "prod"

    conn = Connection(dcc_mode)
    server = conn.dcc_url

    with open("config/samples.tsv") as sample_file:
        h = sample_file.readline().rstrip('\n').split('\t')
        series_ind = h.index("Series")

        samples = []
        for line in sample_file:
            if line.startswith("#"):
                continue
            entries = line.rstrip('\n').split('\t')
            series = entries[series_ind]

            samples.append(series)

    for series in samples:
        replicate_num = 1

        series_data = conn.get(series)

        for d in series_data["related_datasets"]:
            if "ATAC" in d["assay_term_name"]:
                experiment = d["accession"]
                break

        data = conn.get(experiment)

        path = None

        files = data["files"]
        for f in files:

            id = f["@id"]
            if f["file_format"] != "fastq":
                continue
            if "replicate" not in f:
                continue
            if f["replicate"]["biological_replicate_number"] != replicate_num:
                continue
            if "derived_from" in f:
                continue
            if f["output_type"] == "index reads":
                continue

            if f["paired_end"] == "1":
                path = f["cloud_metadata"]["url"]
                break

        if path is None:
            raise ValueError("href not found for fastqs")

        with urlopen(path) as handle:
            with gzip.open(handle, "rt") as f:
                for line in f:
                    if line.startswith("@")
                        instrument, run, _ = line[1:].split(2)
                        batch = f"{instrument}:{run}"
                        break

        print(f"{dataset}\t{batch}")