Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(stats)

params = snakemake@params 
input_paths = snakemake@input
output_paths = snakemake@output
log_paths = snakemake@log

set.seed(params[["seed"]])

countsplit <- function(X, epsilon) {
    # from https://github.com/anna-neufeld/countsplit/blob/main/R/countsplit.R

    if (epsilon <= 0 | epsilon >= 1) {
    stop('Parameter epsilon must be in (0,1)')
    }

    Xtrain <-apply(X,2,function(u) rbinom(n=length(u), size=u, prob=epsilon))
    rownames(Xtrain) <- rownames(X)
    colnames(Xtrain) <- colnames(X)
    Xtest <- X-Xtrain
    rownames(Xtest) <- rownames(X)
    colnames(Xtest) <- colnames(X)

    return(list(train=Xtrain, test=Xtest))
}


proj <- readRDS(file = input_paths[["project_in"]])

counts <- proj@assays$RNA@counts
split <- countsplit(counts, epsilon=params[["split_frac"]])
train <- CreateAssayObject(counts = split$train)
test <- CreateAssayObject(counts = split$test)
proj[["RNA_train"]] <- train
proj[["RNA_test"]] <- test
DefaultAssay(object = proj) <- "RNA_train"

proj <- NormalizeData(proj, normalization.method = "LogNormalize", scale.factor = 10000)
proj <- FindVariableFeatures(proj, selection.method = "vst", nfeatures = 2000)
proj <- ScaleData(proj)

proj <- RunPCA(proj, features = VariableFeatures(object = proj))

proj <- FindNeighbors(proj, dims = 1:30)
proj <- FindClusters(object = proj) 

proj <- RunUMAP(proj, dims = 1:30, return.model = TRUE)

plt <- DimPlot(proj, reduction = "umap", group.by = "seurat_clusters")
ggsave(output_paths[["umap"]], plt, device = "pdf")

saveRDS(proj, file = output_paths[["project_out"]])

sink(type = "message")
sink()