Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

params = snakemake@params 
wildcards = snakemake@wildcards 
input_paths = snakemake@input
output_paths = snakemake@output
log_paths = snakemake@log

set.seed(params[["seed"]])

stallion <- c("1"="#D51F26","2"="#272E6A","3"="#208A42","4"="#89288F","5"="#F47D2B", "6"="#FEE500","7"="#8A9FD1","8"="#C06CAB","19"="#E6C2DC",
               "10"="#90D5E4", "11"="#89C75F","12"="#F37B7D","13"="#9983BD","14"="#D24B27","15"="#3BBCA8", "16"="#6E4B9E","17"="#0C727C", "18"="#7E1416","9"="#D8A767","20"="#3D3D3D")
stallion <- unname(stallion)

plot_fn <- function(object, group, reduction, colors) {
    cats <- sort(unique(object@meta.data[[group]]))
    colors_out <- rep_len(colors, length(cats))
    names(colors_out) <- cats
    DimPlot(object, reduction = reduction, group.by = group, label = TRUE, cols = colors_out, pt.size=0.1)
}

proj <- readRDS(file = input_paths[["project_in"]])

proj_merged <- FindNeighbors(proj_merged, dims = 1:30, reduction = "scbasset_harmony", graph.name = "nn_scbasset")
proj_merged <- RunUMAP(proj_merged, dims = 1:30, reduction = "scbasset_harmony", nn.name = "nn_scbasset", reduction.name = "umap_scbasset")

proj_merged <- FindNeighbors(proj_merged, dims = 1:30, reduction = "cellspace_harmony", graph.name = "nn_cellspace")
proj_merged <- RunUMAP(proj_merged, dims = 1:30, reduction = "cellspace_harmony", nn.name = "nn_cellspace", reduction.name = "umap_cellspace")

plt <- plot_fn(proj, "cell_types_split", "umap_cellspace", stallion) + labs(title="Cell Types (CellSpace)")
ggsave(output_paths[["umap_cellspace"]], plt, device = "pdf", width = 10, height = 7)

plt <- plot_fn(proj, "cell_types_split", "umap_scbasset", stallion) + labs(title="Cell Types (scBasset)")
ggsave(output_paths[["umap_scbasset"]], plt, device = "pdf", width = 10, height = 7)

saveRDS(proj, file = output_paths[["project_out"]])