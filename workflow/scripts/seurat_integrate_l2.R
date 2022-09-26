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

print(params[["groups"]]) ####

projs <- lapply(input_paths[["projects_in"]], readRDS)
# lapply(projs, print) ####

group_metadata <- function(proj, group) {
    print(group) ####
    data <- strsplit(group, split = "-")
    print(data) ####
    # print(length(Cells(proj))) ####
    group <- rep(group, length(Cells(proj)))
    region <- rep(data[[1]], length(Cells(proj)))
    status <- rep(data[[2]], length(Cells(proj)))
    md <- data.frame(group, region, status, row.names = Cells(proj), stringsAsFactors = FALSE)
    print(head(md)) ####
    proj <- AddMetaData(
        object = proj,
        metadata = md,
    )
    proj
}

projs <- mapply(group_metadata, projs, params[["groups"]], SIMPLIFY = FALSE, USE.NAMES = FALSE)

proj_merged <- merge(projs[[1]], projs[-1], project = "merged_l2")

proj_merged <- NormalizeData(proj_merged, normalization.method = "LogNormalize", scale.factor = 10000)
proj_merged <- FindVariableFeatures(proj_merged, selection.method = "vst", nfeatures = 2000)
proj_merged <- ScaleData(proj_merged)

proj_merged <- RunPCA(proj_merged, features = VariableFeatures(object = proj_merged))

proj_merged <- FindNeighbors(proj_merged, dims = 1:30, reduction = "pca")
proj_merged <- RunUMAP(proj_merged, dims = 1:30, reduction = "pca")

proj_merged$mixing_pca <- MixingMetric(
  proj_merged,
  "dataset",
  reduction = "pca",
  dims = 1:30
)

plt <- DimPlot(proj_merged, reduction = "umap", group.by = "dataset") + labs(title="Pre-Harmony datasets")
ggsave(output_paths[["umap_dataset_pre_harmony"]], plt, device = "pdf", width = 10, height = 7)

plt <- FeaturePlot(proj_merged, reduction = "umap", features = "mixing_pca") + labs(title="Pre-Harmony batch mixing metric")
ggsave(output_paths[["umap_mixing_pre_harmony"]], plt, device = "pdf", width = 10, height = 7)

proj_merged <- RunHarmony(proj_merged, "dataset")

proj_merged <- FindNeighbors(proj_merged, dims = 1:30, reduction = "harmony")
proj_merged <- RunUMAP(proj_merged, dims = 1:30, reduction = "harmony")

proj_merged$mixing_harmony <- MixingMetric(
  proj_merged,
  "dataset",
  reduction = "harmony",
  dims = 1:30
)

plt <- DimPlot(proj_merged, reduction = "umap", group.by = "dataset") + labs(title="Post-Harmony datasets")
ggsave(output_paths[["umap_dataset_harmony"]], plt, device = "pdf", width = 10, height = 7)

plt <- FeaturePlot(proj_merged, reduction = "umap", features = "mixing_harmony") + labs(title="Post-Harmony batch mixing metric")
ggsave(output_paths[["umap_mixing_harmony"]], plt, device = "pdf", width = 10, height = 7)

plt <- DimPlot(proj_merged, reduction = "umap", group.by = "group") + labs(title="Post-Harmony groups")
ggsave(output_paths[["umap_group_harmony"]], plt, device = "pdf", width = 10, height = 7)

plt <- DimPlot(proj_merged, reduction = "umap", group.by = "region") + labs(title="Post-Harmony regions")
ggsave(output_paths[["umap_region_harmony"]], plt, device = "pdf", width = 10, height = 7)

plt <- DimPlot(proj_merged, reduction = "umap", group.by = "status") + labs(title="Post-Harmony health status")
ggsave(output_paths[["umap_status_harmony"]], plt, device = "pdf", width = 10, height = 7)

saveRDS(proj_merged, file = output_paths[["project_out"]])

sink(type = "message")
sink()