Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(harmony)

stallion <- c("1"="#D51F26","2"="#272E6A","3"="#208A42","4"="#89288F","5"="#F47D2B", "6"="#FEE500","7"="#8A9FD1","8"="#C06CAB","19"="#E6C2DC",
               "10"="#90D5E4", "11"="#89C75F","12"="#F37B7D","13"="#9983BD","14"="#D24B27","15"="#3BBCA8", "16"="#6E4B9E","17"="#0C727C", "18"="#7E1416","9"="#D8A767","20"="#3D3D3D")
stallion <- unname(stallion)

plot_fn <- function(object, group, reduction, colors) {
    cats <- sort(unique(object@meta.data[[group]]))
    colors_out <- rep_len(colors, length(cats))
    names(colors_out) <- cats
    DimPlot(object, reduction = reduction, group.by = group, label = TRUE, cols = colors_out, pt.size=0.1)
}

params = snakemake@params 
input_paths = snakemake@input
output_paths = snakemake@output
log_paths = snakemake@log

set.seed(params[["seed"]])

proj <- readRDS(file = input_paths[["project_in"]])

DefaultAssay(object = proj) <- "RNA"

proj <- NormalizeData(proj, normalization.method = "LogNormalize", scale.factor = 10000)
proj <- FindVariableFeatures(proj, selection.method = "vst", nfeatures = 2000)
proj <- ScaleData(proj)

proj <- RunPCA(proj, features = VariableFeatures(object = proj))

proj <- FindNeighbors(proj, dims = 1:30, reduction = "pca")
proj <- RunUMAP(proj, dims = 1:30, reduction = "pca")

proj <- RunHarmony(proj, "dataset", assay.use = "RNA")

proj <- FindNeighbors(proj, dims = 1:30, reduction = "harmony")
proj <- RunUMAP(proj, dims = 1:30, reduction = "harmony")


plt <- plot_fn(proj, "dataset", "umap", stallion) + labs(title="All counts embedding: datasets")
ggsave(output_paths[["umap_dataset"]], plt, device = "pdf", width = 10, height = 7)

plt <- plot_fn(proj, "cell_types_l1", "umap", stallion) + labs(title="All counts embedding: Level-1 cell types")
ggsave(output_paths[["umap_cell_types"]], plt, device = "pdf", width = 10, height = 7)


saveRDS(proj, file = output_paths[["project_out"]])

sink(type = "message")
sink()