Sys.setenv(CONDA_BUILD_SYSROOT = "/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(SoupX)

params <- snakemake@params
input_paths <- snakemake@input
output_paths <- snakemake@output
log_paths <- snakemake@log

set.seed(params[["seed"]])

expression_matrix <- ReadMtx(
  mtx = input_paths[["mat"]], 
  features = input_paths[["features"]],
  cells = input_paths[["cells"]]
)
# input_paths[is.na(input_paths)] <- 0

expression_matrix_raw <- ReadMtx(
  mtx = input_paths[["mat_raw"]], 
  features = input_paths[["features_raw"]],
  cells = input_paths[["cells_raw"]]
)
# expression_matrix_raw[is.na(expression_matrix_raw)] <- 0

metadata <- read.table(file = input_paths[["metadata"]], sep = '\t', header = TRUE, quote = "")
rownames(metadata) <- metadata$NAME
print(head(metadata)) ####

clusters <- metadata[colnames(expression_matrix),"cell_type_leiden06", drop=FALSE]
rownames(clusters) <- colnames(expression_matrix)
clusters$cell_type_leiden06[is.na(clusters$cell_type_leiden06)] <- "Unknown"
clusters <- setNames(clusters$cell_type_leiden06, rownames(clusters))

sc <- SoupChannel(expression_matrix_raw, expression_matrix)
sc <- setClusters(sc, clusters)
sc <- autoEstCont(sc)

out <- adjustCounts(sc) 

# Load the proj dataset
expression_matrix <- ReadMtx(
  mtx = input_paths[["mat"]], 
  features = input_paths[["features"]],
  cells = input_paths[["cells"]],
  feature.column = 1,
)

# Initialize the Seurat object with the raw (non-normalized data).
proj <- CreateSeuratObject(
    counts = expression_matrix, 
    project = "reference", 
    min.cells = 3, 
    min.features = 200, 
    meta.data = metadata
)
proj[["percent.mt"]] <- PercentageFeatureSet(proj, pattern = "^MT-")
proj$cell_type <- proj[["sub_cluster"]]

plt <- VlnPlot(proj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(output_paths[["qc_violin"]], plt, device = "pdf")

plot1 <- FeatureScatter(proj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(proj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plt <- plot1 + plot2
ggsave(output_paths[["qc_scatter"]], plt, device = "pdf")

proj <- subset(proj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

proj <- NormalizeData(proj, normalization.method = "LogNormalize", scale.factor = 10000)
proj <- FindVariableFeatures(proj, selection.method = "vst", nfeatures = 2000)
proj <- ScaleData(proj)

proj <- RunPCA(proj, features = VariableFeatures(object = proj))

proj <- FindNeighbors(proj, dims = 1:30)

proj <- RunUMAP(proj, dims = 1:30, return.model = TRUE)

plt <- DimPlot(proj, reduction = "umap", group.by = "cell_type")
ggsave(output_paths[["umap"]], plt, device = "pdf")

saveRDS(proj, file = output_paths[["project_out"]])

sink(type = "message")
sink()