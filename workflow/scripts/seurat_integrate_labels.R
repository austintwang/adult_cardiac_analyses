Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(harmony)

params = snakemake@params 
input_paths = snakemake@input
output_paths = snakemake@output
log_paths = snakemake@log

set.seed(params[["seed"]])

projs <- lapply(input_paths[["projects_labeled"]], readRDS)
proj_labeled <- merge(projs[[1]], projs[-1], project = "rna_labeled", add.cell.ids = unlist(params[["samples"]]))

proj <- readRDS(file = input_paths[["project_integrated"]])
label_cols <- c(
    "cell_type_kramann", 
    "cell_type_ellinor_coarse", 
    "cell_type_ellinor_fine"
)
label_data <- proj_labeled@meta.data[,label_cols]
proj_merged <- AddMetaData(proj, label_data)

plt <- DimPlot(proj_merged, reduction = "umap", group.by = "cell_type_kramann") + labs(title="Kramann labels")
ggsave(output_paths[["umap_kramann"]], plt, device = "pdf", width = 10, height = 7)

plt <- DimPlot(proj_merged, reduction = "umap", group.by = "cell_type_ellinor_coarse") + labs(title="Ellinor coarse labels")
ggsave(output_paths[["umap_ellinor_coarse"]], plt, device = "pdf", width = 10, height = 7)

plt <- DimPlot(proj_merged, reduction = "umap", group.by = "cell_type_ellinor_fine") + labs(title="Ellinor fine labels")
ggsave(output_paths[["umap_ellinor_fine"]], plt, device = "pdf", width = 10, height = 7)


saveRDS(proj_merged, file = output_paths[["project_out"]])

sink(type = "message")
sink()