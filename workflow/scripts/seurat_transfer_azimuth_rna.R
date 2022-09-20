Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

options(future.globals.maxSize = 8000 * 1024^2)

params = snakemake@params 
input_paths = snakemake@input
output_paths = snakemake@output
log_paths = snakemake@log

library(Azimuth, lib.loc=input_paths[["azimuth_library_dir"]])
library(SeuratData, lib.loc=input_paths[["azimuth_library_dir"]])

set.seed(params[["seed"]])

proj <- readRDS(file = input_paths[["project_rna"]])

proj_tmp <- RunAzimuth(query = input_paths[["project_rna"]], reference = "heartref")
proj$cell_type_l1 <- proj_tmp$predicted.celltype.l1
proj$cell_type_l2 <- proj_tmp$predicted.celltype.l2

plt <- DimPlot(proj, reduction = "umap", group.by = "cell_type_l1", label = TRUE)
ggsave(output_paths[["umap_azimuth_l1"]], plt, device = "pdf", width = 10, height = 7)

plt <- DimPlot(proj, reduction = "umap", group.by = "cell_type_l2", label = TRUE)
ggsave(output_paths[["umap_azimuth_l2"]], plt, device = "pdf", width = 10, height = 7)

label_data <- proj@meta.data[,c("cell_type_l1", "cell_type_l2"),drop=FALSE]
write.table(label_data, file= output_paths[["data_out"]], quote=FALSE, sep='\t', col.names = NA)
