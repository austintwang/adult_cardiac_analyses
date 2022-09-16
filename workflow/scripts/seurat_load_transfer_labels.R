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

read_fn <- function(path) {
   read.table(file = path, sep = '\t', header = TRUE, row.names = 1)
}

tables <- lapply(input_paths[["tables"]], read_fn)
head(tables[[1]])
label_data <- do.call(cbind, tables)
# head(label_data) ####

proj <- readRDS(file = input_paths[["project_integrated"]])
proj_merged <- AddMetaData(proj, label_data)

plt <- DimPlot(proj_merged, reduction = "umap", group.by = "cell_type_kramann", label = TRUE) + labs(title="Kramann labels")
ggsave(output_paths[["umap_kramann"]], plt, device = "pdf", width = 10, height = 7)

plt <- DimPlot(proj_merged, reduction = "umap", group.by = "cell_type_ellinor_coarse", label = TRUE) + labs(title="Ellinor coarse labels")
ggsave(output_paths[["umap_ellinor_coarse"]], plt, device = "pdf", width = 10, height = 7)

plt <- DimPlot(proj_merged, reduction = "umap", group.by = "cell_type_ellinor_fine", label = TRUE) + labs(title="Ellinor fine labels")
ggsave(output_paths[["umap_ellinor_fine"]], plt, device = "pdf", width = 10, height = 7)

plt <- DimPlot(proj_merged, reduction = "umap", group.by = "cell_type_teichmann_coarse", label = TRUE) + labs(title="Teichmann coarse labels")
ggsave(output_paths[["umap_teichmann_coarse"]], plt, device = "pdf", width = 10, height = 7)

plt <- DimPlot(proj_merged, reduction = "umap", group.by = "cell_type_teichmann_fine", label = TRUE) + labs(title="Teichmann fine labels")
ggsave(output_paths[["umap_teichmann_fine"]], plt, device = "pdf", width = 10, height = 7)


saveRDS(proj_merged, file = output_paths[["project_out"]])

sink(type = "message")
sink()