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

plt <- DimPlot(proj_merged, reduction = "umap", group.by = "cell_type_kramann") + labs(title="Pre-Harmony, Kramann labels")
ggsave(output_paths[["umap_kramann_pre_harmony"]], plt, device = "pdf")

plt <- DimPlot(proj_merged, reduction = "umap", group.by = "cell_type_ellinor_coarse") + labs(title="Pre-Harmony, Ellinor coarse labels")
ggsave(output_paths[["umap_ellinor_coarse_pre_harmony"]], plt, device = "pdf")

plt <- DimPlot(proj_merged, reduction = "umap", group.by = "cell_type_ellinor_fine") + labs(title="Pre-Harmony, Ellinor fine labels")
ggsave(output_paths[["umap_ellinor_fine_pre_harmony"]], plt, device = "pdf")

plt <- DimPlot(proj_merged, reduction = "umap", group.by = "dataset") + labs(title="Pre-Harmony datasets")
ggsave(output_paths[["umap_dataset_pre_harmony"]], plt, device = "pdf")

plt <- DimPlot(proj_merged, reduction = "umap", group.by = "mixing_pca") + labs(title="Pre-Harmony batch mixing metric")
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

plt <- DimPlot(proj_merged, reduction = "umap", group.by = "cell_type_kramann") + labs(title="Post-Harmony, Kramann labels")
ggsave(output_paths[["umap_kramann_harmony"]], plt, device = "pdf")

plt <- DimPlot(proj_merged, reduction = "umap", group.by = "cell_type_ellinor_coarse") + labs(title="Post-Harmony, Ellinor coarse labels")
ggsave(output_paths[["umap_ellinor_coarse_harmony"]], plt, device = "pdf")

plt <- DimPlot(proj_merged, reduction = "umap", group.by = "cell_type_ellinor_fine") + labs(title="Post-Harmony, Ellinor fine labels")
ggsave(output_paths[["umap_ellinor_fine_harmony"]], plt, device = "pdf")

plt <- DimPlot(proj_merged, reduction = "umap", group.by = "dataset") + labs(title="Post-Harmony datasets")
ggsave(output_paths[["umap_dataset_harmony"]], plt, device = "pdf")

plt <- DimPlot(proj_merged, reduction = "umap", group.by = "mixing_harmony") + labs(title="Post-Harmony batch mixing metric")
ggsave(output_paths[["umap_mixing_harmony"]], plt, device = "pdf")

saveRDS(proj_merged, file = output_paths[["project_out"]])

sink(type = "message")
sink()