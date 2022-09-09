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
library(MatrixExtra)

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
rownames(metadata) <- metadata$X
print(head(metadata)) ####
print(head(colnames(expression_matrix))) ####
print(head(rownames(expression_matrix))) ####

clusters <- metadata[colnames(expression_matrix),"seurat_clusters", drop=FALSE]
rownames(clusters) <- colnames(expression_matrix)
clusters$seurat_clusters[is.na(clusters$seurat_clusters)] <- -1
print(head(clusters)) ####
print(length(clusters)) ####
print(length(colnames(expression_matrix))) ####
clusters <- setNames(clusters$seurat_clusters, rownames(clusters))
print(all(colnames(expression_matrix) %in% names(clusters))) ####
# clusters ####

sc <- SoupChannel(expression_matrix_raw, expression_matrix)
sc ####
sc <- setClusters(sc, clusters)
sc ####
sc <- autoEstCont(sc)
sc ####
# s <- split(colnames(sc$toc),clusters[colnames(sc$toc)]) ####
# print(s) ####
# # print(sc$toc) ####
# print(s$`0`) ####
# b <- s$`0`[1:3]
# print(typeof(s$`0`))
# print(typeof(colnames(sc$toc))) ####
# print(all(s$`0` %in% colnames(sc$toc))) ####
# print(b) ####
# print(sc$toc[,b]) ####
# print(sc$toc[,s$`0`,drop=FALSE]) ####
# print(rowSums(sc$toc[,s[[2]],drop=FALSE])) ####
# a <- do.call(cbind,lapply(s,function(e) rowSums(sc$toc[,e,drop=FALSE]))) ####
# print(a) ####
print(head(sc$metaData)) ####
print(sc$soupProfile$est) ####
out <- as(sc$toc,'dgTMatrix') ####
expSoupCnts <- sc$metaData$nUMIs * sc$metaData$rho ####
soupFrac <- sc$soupProfile$est ####
out <- out - do.call(cbind,lapply(seq(ncol(out)),function(e) alloc(expSoupCnts[e],out[,e],soupFrac))) ####

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