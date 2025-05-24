import polars as pl

def main(counts_path, out_path, samples):
    with open(counts_path, "r") as f:
        header = f.readline().strip().split("\t")[1:]
    missing_samples = set(samples) - set(header)
    
    df = (
        pl.scan_csv(counts_path, separator="\t")
        .rename({"": "Feature"})
        .with_columns([pl.lit(0).alias(s) for s in missing_samples])
        .select(["Feature", *samples])
        .sink_csv(out_path, include_header=False, separator="\t")
    )

counts_path = snakemake.input.counts
out_path = snakemake.output.counts
samples = snakemake.params.samples

main(counts_path, out_path, samples)