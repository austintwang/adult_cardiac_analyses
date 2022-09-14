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


projs <- lapply(input_paths[["projects_in"]], readRDS)
lapply(projs, print) ####

proj_merged <- merge(projs[[1]], projs[-1], project = "merged_rna", add.cell.ids = unlist(params[["samples"]]))

proj_merged <- SCTransform(proj_merged, vars.to.regress = "percent.mt", verbose = FALSE)
proj_merged <- RunPCA(proj_merged)

proj_merged <- FindNeighbors(proj_merged, dims = 1:30, reduction = "pca")
proj_merged <- RunUMAP(proj_merged, dims = 1:30, reduction = "pca")

proj_merged$mixing_pca <- MixingMetric(
  proj_merged,
  "dataset",
  reduction = "pca",
  dims = 1:30
)

plt <- DimPlot(proj_merged, reduction = "umap", group.by = "cell_type_ref")
ggsave(output_paths[["umap_clusters_pre_harmony"]], plt, device = "pdf")

plt <- DimPlot(proj_merged, reduction = "umap", group.by = "dataset")
ggsave(output_paths[["umap_dataset_pre_harmony"]], plt, device = "pdf")

plt <- DimPlot(proj_merged, reduction = "umap", group.by = "mixing_pca")
ggsave(output_paths[["umap_mixing_pre_harmony"]], plt, device = "pdf")

proj_merged <- RunHarmony(proj_merged, "dataset")

proj_merged <- FindNeighbors(proj_merged, dims = 1:30, reduction = "harmony")
proj_merged <- RunUMAP(proj_merged, dims = 1:30, reduction = "harmony")

proj_merged$mixing_harmony <- MixingMetric(
  proj_merged,
  "dataset",
  reduction = "harmony",
  dims = 1:30
)

plt <- DimPlot(proj_merged, reduction = "umap", group.by = "cell_type_ref")
ggsave(output_paths[["umap_clusters_harmony"]], plt, device = "pdf")

plt <- DimPlot(proj_merged, reduction = "umap", group.by = "dataset")
ggsave(output_paths[["umap_dataset_harmony"]], plt, device = "pdf")

plt <- DimPlot(proj_merged, reduction = "umap", group.by = "mixing_harmony")
ggsave(output_paths[["umap_mixing_harmony"]], plt, device = "pdf")

saveRDS(proj_merged, file = output_paths[["project_out"]])

sink(type = "message")
sink()