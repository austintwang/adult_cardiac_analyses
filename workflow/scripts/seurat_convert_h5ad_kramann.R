Sys.setenv(CONDA_BUILD_SYSROOT = "/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(reticulate)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(Matrix)

params <- snakemake@params
input_paths <- snakemake@input
output_paths <- snakemake@output
log_paths <- snakemake@log

set.seed(params[["seed"]])

sc <- import("scanpy")
adata <- sc$read_h5ad(input_paths[["h5ad"]])
adata ####

# counts <- t(adata$layers["counts"])
# typeof(adata$X) ####
counts <- t(adata$X)
colnames(counts) <- adata$obs_names$to_list()
rownames(counts) <- adata$var_names$to_list()
counts <- Matrix::Matrix(as.matrix(counts), sparse = TRUE)

proj <- CreateSeuratObject(counts)
proj <- AddMetaData(proj, adata$obs)

print(head(rownames(proj))) ####
head(proj@meta.data) ####

saveRDS(proj, file = output_paths[["project_out"]])

sink(type = "message")
sink()