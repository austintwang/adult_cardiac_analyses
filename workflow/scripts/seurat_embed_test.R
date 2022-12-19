Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cluster)
library(harmony)

params = snakemake@params 
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
print(proj) ####

DefaultAssay(object = proj_merged) <- "RNA_test"

proj_merged <- NormalizeData(proj_merged, normalization.method = "LogNormalize", scale.factor = 10000)
proj_merged <- FindVariableFeatures(proj_merged, selection.method = "vst", nfeatures = 2000)
proj_merged <- ScaleData(proj_merged)

proj_merged <- RunPCA(proj_merged, features = VariableFeatures(object = proj_merged), reduction.name = "pca_test")

proj_merged <- FindNeighbors(proj_merged, dims = 1:30, reduction = "pca_test")
proj_merged <- RunUMAP(proj_merged, dims = 1:30, reduction = "pca_test")

proj_merged <- RunHarmony(proj_merged, "dataset", assay.use = "RNA_test", reduction.save = "harmony_test")

proj_merged <- FindNeighbors(proj_merged, dims = 1:30, reduction = "harmony_test", graph.name = "nn_test")
proj_merged <- RunUMAP(proj_merged, dims = 1:30, reduction = "harmony_test", nn.name = "nn_test", reduction.name = "umap_test")

plt <- plot_fn(proj, "seurat_clusters", "umap_test", stallion) + labs(title="Seurat clustering (test split)")
ggsave(output_paths[["umap"]], plt, device = "pdf")

dims <- 1:30

dist.matrix <- dist(x = Embeddings(object = dataset[["harmony_test"]])[, dims])
clusters <- proj$seurat_clusters
sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
pdf(output_paths[["sil_test"]])
plot(sil)
dev.off()

dist.matrix <- dist(x = Embeddings(object = dataset[["harmony_train"]])[, dims])
clusters <- proj$seurat_clusters
sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
pdf(output_paths[["sil_train"]])
plot(sil)
dev.off()

DefaultAssay(object = proj_merged) <- "RNA_train"

saveRDS(proj, file = output_paths[["project_out"]])