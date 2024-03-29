Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(SoupX)
library(Matrix)
# library(MatrixExtra)

params = snakemake@params 
wildcards = snakemake@wildcards
input_paths = snakemake@input
output_paths = snakemake@output
log_paths = snakemake@log

set.seed(params[["seed"]])

# Load the proj dataset
expression_matrix <- ReadMtx(
  mtx = input_paths[["mat"]], 
  features = input_paths[["features"]],
  cells = input_paths[["cells"]]
)
# colnames(expression_matrix) <- noquote(colnames(expression_matrix)) ####
# rownames(expression_matrix) <- noquote(rownames(expression_matrix)) ####

expression_matrix_raw <- ReadMtx(
  mtx = input_paths[["mat_raw"]], 
  features = input_paths[["features_raw"]],
  cells = input_paths[["cells_raw"]]
)
# colnames(expression_matrix_raw) <- noquote(colnames(expression_matrix_raw)) ####
# rownames(expression_matrix_raw) <- noquote(rownames(expression_matrix_raw)) ####
# expression_matrix ####

metadata <- read.table(file = input_paths[["metadata"]], sep = '\t', header = TRUE)
print(head(metadata)) ####
metadata <- metadata[metadata$dataset == wildcards[["sample"]],]
rownames(metadata) <- metadata$barcode_rna
print(head(metadata)) ####
print(head(colnames(expression_matrix))) ####
print(head(rownames(expression_matrix))) ####

clusters <- metadata[colnames(expression_matrix),"cell_types_l1", drop=FALSE]
rownames(clusters) <- colnames(expression_matrix)
clusters$cell_types_l1[is.na(clusters$cell_types_l1)] <- -1
print(head(clusters)) ####
print(length(clusters)) ####
print(length(colnames(expression_matrix))) ####
clusters <- setNames(clusters$cell_types_l1, rownames(clusters))
# print(all(colnames(expression_matrix) %in% names(clusters))) ####
# clusters ####

sc <- SoupChannel(expression_matrix_raw, expression_matrix)
# sc ####
sc <- setClusters(sc, clusters)
# sc ####
pdf(output_paths[["rho"]])
sc <- autoEstCont(sc, forceAccept = TRUE)
title(wildcards[["sample"]])
dev.off()

sc ####

out <- adjustCounts(sc) 
# head(out) ####

# Initialize the Seurat object with the raw (non-normalized data).
proj <- CreateSeuratObject(counts = out, project = wildcards[["sample"]], min.cells = 3, min.features = 0, meta.data = metadata)
proj <- proj[,colnames(proj) %in% rownames(metadata)]
proj <- AddMetaData(object = proj, metadata = sc$metaData)
print(proj) ####

proj <- NormalizeData(proj, normalization.method = "LogNormalize", scale.factor = 10000)
proj <- FindVariableFeatures(proj, selection.method = "vst", nfeatures = 2000)
proj <- ScaleData(proj)

proj <- RunPCA(proj, features = VariableFeatures(object = proj))

proj <- FindNeighbors(proj, dims = 1:30)
proj <- FindClusters(object = proj) 

proj <- RunUMAP(proj, dims = 1:30, return.model = TRUE)

# plt <- FeaturePlot(proj, reduction = "umap", features = "rho")
# ggsave(output_paths[["umap"]], plt, device = "pdf")

saveRDS(proj, file = output_paths[["project_out"]])

sink(type = "message")
sink()