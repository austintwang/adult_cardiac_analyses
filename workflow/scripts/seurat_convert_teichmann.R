Sys.setenv(CONDA_BUILD_SYSROOT = "/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

params <- snakemake@params
input_paths <- snakemake@input
output_paths <- snakemake@output
log_paths <- snakemake@log

set.seed(params[["seed"]])

Convert(input_paths[["h5ad"]], dest = output_paths[["h5seurat"]])

library(SeuratData, lib.loc=input_paths[["azimuth_library_dir"]])
library(SeuratDisk, lib.loc=input_paths[["seuratdisk_library_dir"]])

sink(type = "message")
sink()