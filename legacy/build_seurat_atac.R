Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

params = snakemake@params 
input_paths = snakemake@input
output_paths = snakemake@output
log_paths = snakemake@log

set.seed(params[["seed"]])

# Load the proj dataset
peak_matrix <- ReadMtx(
  mtx = input_paths[["mat"]], 
  features = input_paths[["features"]],
  cells = input_paths[["cells"]],
  feature.column = 1,
)
# Initialize the Seurat object with the raw (non-normalized data).

chrom_assay <- CreateChromatinAssay(
  counts = peak_matrix,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = input_paths[["fragments"]],
  min.cells = 10,
  min.features = 200
)

proj <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks"
)

saveRDS(proj, file = output_paths[["project_out"]])

sink(type = "message")
sink()