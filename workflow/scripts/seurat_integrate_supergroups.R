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
wildcards = snakemake@wildcards
input_paths = snakemake@input
output_paths = snakemake@output
log_paths = snakemake@log

set.seed(params[["seed"]])

load_proj <- function(path) {
    proj <- readRDS(path)
    # proj <- subset(x = proj, subset = cell_types_1 == wildcards[["label"]])
    # if (length(WhichCells(proj, expression = cell_types_1 == wildcards[["label"]])) > 0) {
    #     proj <- subset(x = proj, subset = cell_types_1 == wildcards[["label"]])
    # } else {
    #     proj <- subset(x = proj, downsample = 1)
    # }
    # proj <- tryCatch(
    #     expr = subset(x = proj, subset = cell_types_1 == wildcards[["label"]]),
    #     error = function(e) {
    #         print(path)
    #         print(e)
    #         subset(x = proj, downsample = 1)
    #     }
    # )
    proj
}

projs <- lapply(input_paths[["projects_in"]], load_proj)
lapply(projs, print) ####

# group_metadata <- function(group, proj) {
#     group <- group[[1]]
#     print(group) ####
#     data <- unlist(strsplit(group, split = "-"))
#     print(data) ####
#     # print(length(Cells(proj))) ####
#     ncells <- length(Cells(proj))
#     group <- rep(group, ncells)
#     region <- rep(data[[1]], ncells)
#     status <- rep(data[[2]], ncells)
#     md <- data.frame(group, region, status, row.names = Cells(proj), stringsAsFactors = FALSE)
#     print(head(md)) ####
#     md
# }
# group_mds <- mapply(group_metadata, params[["groups"]], projs, SIMPLIFY = FALSE, USE.NAMES = FALSE)
# group_md <- do.call(rbind, group_mds)
# head(group_md)

if (length(projs) > 1) {
  proj <- merge(projs[[1]], projs[-1], project = "merged_rna_subcluster_unified")
} else {
  proj <- projs[[1]]
#   proj <- RenameCells(proj, add.cell.id = params[["samples"]][[1]], for.merge = TRUE)
}
# proj <- AddMetaData(proj, group_md)
# proj <- subset(x = proj, subset = cell_types_1 == wildcards[["label"]])
DefaultAssay(object = proj) <- "RNA_train"

proj <- NormalizeData(proj, normalization.method = "LogNormalize", scale.factor = 10000)
proj <- FindVariableFeatures(proj, selection.method = "vst", nfeatures = 2000)
proj <- ScaleData(proj)

proj <- RunPCA(proj, features = VariableFeatures(object = proj))

proj <- FindNeighbors(proj, dims = 1:30, reduction = "pca")
proj <- RunUMAP(proj, dims = 1:30, reduction = "pca")

if (length(projs) > 1) {
#   proj$mixing_pca <- MixingMetric(
#     proj,
#     "dataset",
#     reduction = "pca",
#     dims = 1:30
#   )
    proj$mixing_pca <- tryCatch(
        expr = MixingMetric(
            proj,
            "dataset",
            reduction = "pca",
            dims = 1:30
        ),
        error = function(e) {
            print(e)
            MixingMetric(
                proj,
                "dataset",
                reduction = "pca",
                dims = 1:30,
                max.k = 30
            )
        }
    )
} else {
  proj <- AddMetaData(proj, rep(0, length(Cells(proj))), col.name = "mixing_pca")
}

plt <- DimPlot(proj, reduction = "umap", group.by = "dataset") + labs(title="Pre-Harmony datasets")
ggsave(output_paths[["umap_dataset_pre_harmony"]], plt, device = "pdf", width = 10, height = 7)

plt <- FeaturePlot(proj, reduction = "umap", features = "mixing_pca") + labs(title="Pre-Harmony batch mixing metric")
ggsave(output_paths[["umap_mixing_pre_harmony"]], plt, device = "pdf", width = 10, height = 7)

proj <- RunHarmony(proj, "dataset", assay.use = "RNA_train")

proj <- FindNeighbors(proj, dims = 1:30, reduction = "harmony")
proj <- RunUMAP(proj, dims = 1:30, reduction = "harmony")

DefaultAssay(object = proj) <- "RNA_test"

proj <- NormalizeData(proj, normalization.method = "LogNormalize", scale.factor = 10000)
proj <- FindVariableFeatures(proj, selection.method = "vst", nfeatures = 2000)
proj <- ScaleData(proj)

proj <- RunPCA(proj, features = VariableFeatures(object = proj), reduction.name = "pca_test")

proj <- RunHarmony(proj, "dataset", reduction = "pca_test", assay.use = "RNA_test", reduction.save = "harmony_test")

proj <- FindNeighbors(proj, dims = 1:30, reduction = "harmony_test", graph.name = "nn_test")
proj <- RunUMAP(proj, dims = 1:30, reduction = "harmony_test", nn.name = "nn_test", reduction.name = "umap_test")

DefaultAssay(object = proj) <- "RNA_train"

if (length(projs) > 1) {
#   proj$mixing_harmony <- MixingMetric(
#     proj,
#     "dataset",
#     reduction = "harmony",
#     dims = 1:30
#   )
  proj$mixing_harmony <- tryCatch(
        expr = MixingMetric(
            proj,
            "dataset",
            reduction = "harmony",
            dims = 1:30
        ),
        error = function(e) {
            print(e)
            MixingMetric(
                proj,
                "dataset",
                reduction = "harmony",
                dims = 1:30,
                max.k = 30
            )
        }
    )
} else {
  proj <- AddMetaData(proj, rep(0, length(Cells(proj))), col.name = "mixing_harmony")
}

plt <- DimPlot(proj, reduction = "umap", group.by = "dataset") + labs(title="Post-Harmony datasets")
ggsave(output_paths[["umap_dataset_harmony"]], plt, device = "pdf", width = 10, height = 7)

plt <- FeaturePlot(proj, reduction = "umap", features = "mixing_harmony") + labs(title="Post-Harmony batch mixing metric")
ggsave(output_paths[["umap_mixing_harmony"]], plt, device = "pdf", width = 10, height = 7)

saveRDS(proj, file = output_paths[["project_out"]])

sink(type = "message")
sink()