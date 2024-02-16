Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(Seurat)
library(Matrix)

params = snakemake@params 
input_paths = snakemake@input
output_paths = snakemake@output
log_paths = snakemake@log

set.seed(params[["seed"]])


proj <- readRDS(file = input_paths[["project_in"]])

counts <- GetAssayData(object = proj, assay.type = "RNA", slot = "counts")

writeMM(Matrix(counts, sparse = TRUE), output_paths[["counts_mat"]])
writeLines(rownames(counts), con = output_paths[["counts_rows"]])
writeLines(colnames(counts), con = output_paths[["counts_cols"]])
