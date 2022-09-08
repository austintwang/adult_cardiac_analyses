Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(SoupX)

params = snakemake@params 
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
# expression_matrix ####
expression_matrix_raw <- ReadMtx(
  mtx = input_paths[["mat_raw"]], 
  features = input_paths[["features_raw"]],
  cells = input_paths[["cells_raw"]]
)
# expression_matrix ####

metadata <- read.table(file = input_paths[["metadata"]], sep = '\t', header = TRUE)
print(head(metadata)) ####
clusters <- metadata[colnames(expression_matrix),"seurat_clusters", drop=FALSE]
print(head(clusters)) ####
print(length(colnames(expression_matrix))) ####
print(all(colnames(expression_matrix) %in% names(clusters))) ####

sc <- SoupChannel(expression_matrix_raw, expression_matrix)
sc ####
sc <- setClusters(sc, clusters)
# sc ####
sc <- autoEstCont(sc)
# sc ####
out <- adjustCounts(sc)
# head(out) ####

# Initialize the Seurat object with the raw (non-normalized data).
proj <- CreateSeuratObject(counts = out, project = params[["sample_name"]], min.cells = 3, min.features = 0, meta.data = metadata)
print(proj) ####

proj <- NormalizeData(proj, normalization.method = "LogNormalize", scale.factor = 10000)
proj <- FindVariableFeatures(proj, selection.method = "vst", nfeatures = 2000)
proj <- ScaleData(proj)

proj <- RunPCA(proj, features = VariableFeatures(object = proj))

proj <- FindNeighbors(proj, dims = 1:30)
proj <- FindClusters(object = proj) 

proj <- RunUMAP(proj, dims = 1:30, return.model = TRUE)
ggsave(output_paths[["umap"]], plt, device = "pdf")

saveRDS(proj, file = output_paths[["project_out"]])

sink(type = "message")
sink()