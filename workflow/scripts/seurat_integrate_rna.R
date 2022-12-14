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

if (length(projs) > 1) {
  proj_merged <- merge(projs[[1]], projs[-1], project = "merged_rna", add.cell.ids = unlist(params[["samples"]]))
} else {
  proj_merged <- projs[[1]]
  proj_merged <- RenameCells(proj_merged, add.cell.id = params[["samples"]][[1]], for.merge = TRUE)
}
DefaultAssay(object = proj_merged) <- "RNA_train"

proj_merged <- NormalizeData(proj_merged, normalization.method = "LogNormalize", scale.factor = 10000)
proj_merged <- FindVariableFeatures(proj_merged, selection.method = "vst", nfeatures = 2000)
proj_merged <- ScaleData(proj_merged)

proj_merged <- RunPCA(proj_merged, features = VariableFeatures(object = proj_merged))

proj_merged <- FindNeighbors(proj_merged, dims = 1:30, reduction = "pca")
proj_merged <- RunUMAP(proj_merged, dims = 1:30, reduction = "pca")

if (length(projs) > 1) {
  proj_merged$mixing_pca <- MixingMetric(
    proj_merged,
    "dataset",
    reduction = "pca",
    dims = 1:30
  )
} else {
  proj_merged <- AddMetaData(proj_merged, rep(0, length(Cells(proj_merged))), col.name = "mixing_pca")
}

plt <- DimPlot(proj_merged, reduction = "umap", group.by = "dataset") + labs(title="Pre-Harmony datasets")
ggsave(output_paths[["umap_dataset_pre_harmony"]], plt, device = "pdf", width = 10, height = 7)

plt <- FeaturePlot(proj_merged, reduction = "umap", features = "mixing_pca") + labs(title="Pre-Harmony batch mixing metric")
ggsave(output_paths[["umap_mixing_pre_harmony"]], plt, device = "pdf", width = 10, height = 7)

proj_merged <- RunHarmony(proj_merged, "dataset", assay.use = "RNA_train")

proj_merged <- FindNeighbors(proj_merged, dims = 1:30, reduction = "harmony")
proj_merged <- RunUMAP(proj_merged, dims = 1:30, reduction = "harmony")

if (length(projs) > 1) {
  proj_merged$mixing_harmony <- MixingMetric(
    proj_merged,
    "dataset",
    reduction = "harmony",
    dims = 1:30
  )
} else {
  proj_merged <- AddMetaData(proj_merged, rep(0, length(Cells(proj_merged))), col.name = "mixing_harmony")
}

plt <- DimPlot(proj_merged, reduction = "umap", group.by = "dataset") + labs(title="Post-Harmony datasets")
ggsave(output_paths[["umap_dataset_harmony"]], plt, device = "pdf", width = 10, height = 7)

plt <- FeaturePlot(proj_merged, reduction = "umap", features = "mixing_harmony") + labs(title="Post-Harmony batch mixing metric")
ggsave(output_paths[["umap_mixing_harmony"]], plt, device = "pdf", width = 10, height = 7)

saveRDS(proj_merged, file = output_paths[["project_out"]])

sink(type = "message")
sink()