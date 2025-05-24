import subprocess
from subprocess import PIPE
import os
import io

import polars as pl

RASQUAL_SCHEMA = ["element", "variant", "chrom", "pos", "ref", "alt", "freq", "hwe_chisq", "imp_qual", "log10q", "chisq", 
                  "pi", "delta", "phi", "overdisp", "snp_idx", "n_fsnps", "n_snps", "iter_null", "iter_alt", "tie_loc", 
                  "ll_null", "convergence", "r2_fsnp", "r2_rsnp"]
RASQUAL_DTYPES = [pl.Utf8, pl.Utf8, pl.Utf8, pl.Int32, pl.Utf8, pl.Utf8, pl.Float32, pl.Float32, pl.Float32, pl.Float32, 
                  pl.Float32, pl.Float32, pl.Float32, pl.Float32, pl.Float32, pl.Float32, pl.Int32, pl.Int32, pl.Int32, 
                  pl.Int32, pl.Int32, pl.Float32, pl.Int32, pl.Float32, pl.Float32]

NARROWPEAK_SCHEMA = ["chr", "peak_start", "peak_end", "peak_name", "peak_score", 
                     "peak_strand", "peak_signal", "peak_pval", "peak_qval", "peak_summit"]
NARROWPEAK_DTYPES = [pl.Utf8, pl.UInt32, pl.UInt32, pl.Utf8, pl.UInt32, 
                     pl.Utf8, pl.Float32, pl.Float32, pl.Float32, pl.UInt32] 

def load_peaks(peaks_path, half_width):
    peaks = (
        pl.scan_csv(peaks_path, has_header=False, new_columns=NARROWPEAK_SCHEMA, 
                    separator='\t', quote_char=None, dtypes=NARROWPEAK_DTYPES)
        .select(
            chr=pl.col("chr"),
            peak_region_start=pl.col("peak_start") + pl.col("peak_summit") - half_width,
            peak_region_end=pl.col("peak_start") + pl.col("peak_summit") + half_width,
        )
        .collect()
    )
    return peaks

def run_rasqual_region(bin_path, vcf_path, total_counts_path, offsets_path, chrom, start, end, sample_size, ind):
    tabix_cmd = ["tabix", vcf_path, f"{chrom}:{start+1}-{end}"]
    tabix_cmd = [str(x) for x in tabix_cmd]
    tabix_process = subprocess.Popen(tabix_cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE, text=True)
    vcf_snippet, error = tabix_process.communicate()
    num_snps = vcf_snippet.count("\n")
    # print(num_snps) ####
    # print(" ".join(tabix_process.args)) ####

    rasqual_cmd = [bin_path, "-y", total_counts_path, "-k", offsets_path, "-n", sample_size, "-j", ind+1, 
                   "-l", num_snps, "-m", num_snps, "-s", start+1, "-e", end, "-f", f"{chrom}:{start}-{end}"]
    rasqual_cmd = [str(x) for x in rasqual_cmd]
    rasqual_env = {"LD_LIBRARY_PATH": os.path.join(os.environ["CONDA_PREFIX"], "lib")}
    rasqual_process = subprocess.Popen(rasqual_cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE, text=True, env=rasqual_env)
    output, error = rasqual_process.communicate(input=vcf_snippet)
    # print(output, error) ####
    if error:
        raise subprocess.CalledProcessError(error)
    
    rasqual_df = pl.read_csv(io.StringIO(output), has_header=False, separator='\t', new_columns=RASQUAL_SCHEMA, dtypes=RASQUAL_DTYPES)
    
    return rasqual_df

def main(bin_path, vcf_path, total_counts_path, offsets_path, peaks_path, out_path, samples):
    sample_size = len(samples)
    peaks = load_peaks(peaks_path, 250)
    dfs = []
    for ind, (chrom, start, end) in enumerate(peaks.iter_rows()):
        df = run_rasqual_region(bin_path, vcf_path, total_counts_path, offsets_path, chrom, start, end, sample_size, ind)
        dfs.append(df)

    rasqual_df = pl.concat(dfs)
    rasqual_df.write_csv(out_path, separator='\t')


bin_path = os.path.join(snakemake.input.rasqual_dir, "bin/rasqual")
vcf_path = snakemake.input.vcf
total_counts_path = snakemake.input.total_counts
offsets_path = snakemake.input.offsets
peaks_path = snakemake.input.peaks

out_path = snakemake.output[0]

samples = snakemake.params.samples

main(bin_path, vcf_path, total_counts_path, offsets_path, peaks_path, out_path, samples)