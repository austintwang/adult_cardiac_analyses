import json
import gzip

import pysam
import polars as pl


def get_atac_dataset(dataset_info_json):
    with open(dataset_info_json, "r") as f:
        dataset_info = json.load(f)

    return list(dataset_info.values())[0]["atac"]


def build_bc_set(metadata_path, labels_path, target_cell_type, target_atac_sample):
    with gzip.open(metadata_path, "rb") as f, gzip.open(labels_path, "rb") as g:
        metadata_df = (
            pl.read_csv(f, separator="\t", comment_prefix="#", null_values="NA").lazy()
            .filter(
                (pl.col("passed_filtering") == 1)
                & (pl.col("atac_dataset") == target_atac_sample)
            )
            .select(["cell_id", "atac_barcode"])
        )
        labels_df = (
            pl.read_csv(g, separator="\t", comment_prefix="#", null_values="NA").lazy()
            .filter(pl.col("cell_type_id") == target_cell_type)
            .select(["cell_id"])
        )
        combined_df = metadata_df.join(labels_df, on="cell_id", how="inner").collect()
        bc_set = set(combined_df["atac_barcode"].to_list())

    return bc_set


def subset_bam(input_bam, output_bam, bc_set):
    with pysam.AlignmentFile(input_bam, "rb") as input_bamfile:
        with pysam.AlignmentFile(output_bam, "wb", header=input_bamfile.header) as output_bamfile:
            for read in input_bamfile:
                barcode = read.get_tag("CB")
                if barcode in bc_set:
                    output_bamfile.write(read)


def main(metadata_path, labels_path, dataset_info_json, input_bam, output_bam, cell_type):
    atac_dataset = get_atac_dataset(dataset_info_json)
    bc_set = build_bc_set(metadata_path, labels_path, cell_type, atac_dataset)
    subset_bam(input_bam, output_bam, bc_set)


metadata_path = snakemake.input.metadata
labels_path = snakemake.input.labels
dataset_info_json = snakemake.input.dataset_info
input_bam = snakemake.input.bam

output_bam = snakemake.output.bam

cell_type = snakemake.wildcards.label

main(metadata_path, labels_path, dataset_info_json, input_bam, output_bam, cell_type)