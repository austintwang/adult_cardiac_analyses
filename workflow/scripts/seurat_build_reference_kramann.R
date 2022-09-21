Sys.setenv(CONDA_BUILD_SYSROOT = "/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
# library(biomaRt)
library(harmony)

params <- snakemake@params
input_paths <- snakemake@input
output_paths <- snakemake@output
log_paths <- snakemake@log

library(SeuratData, lib.loc=input_paths[["azimuth_library_dir"]])
library(SeuratDisk, lib.loc=input_paths[["seuratdisk_library_dir"]])

set.seed(params[["seed"]])

proj <- LoadH5Seurat(input_paths[["h5seurat"]])

print(head(rownames(proj))) ####
head(proj@meta.data) ####

# proj@assays$RNA@key <- "rna_"
# # print(Key(proj)) ####
# proj <- UpdateSeuratObject(proj)
Project(proj) <- "Kramann"

cell_type_fine <- proj@meta.data[, "cell_type_original", drop = FALSE]
names(cell_type_fine)["cell_type_original"] <- "cell_type_fine"
for (i in seq_along(params[["subtypes"]])) {
  sub_path <- input_paths[["subtypes"]][[i]]
  sub_name <- params[["subtypes"]][[i]]
  sub_proj <- readRDS(sub_path)
  
  shared <- intersect(Cells(sub_proj), Cells(proj))
  subtypes <- sub_proj@meta.data[shared, "annotation"]
  cell_type_fine[shared, "cell_type_fine"] <- subtypes
}




# ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
# res <- getBM(
#     attributes=c("ensembl_gene_id","hgnc_symbol"), 
#     filters = 'ensembl_gene_id', 
#     values = rownames(proj),
#     mart = ensembl
# )
# # head(res) ####
# res <- res[res[["hgnc_symbol"]] != "",]
# # head(res) ####
# # rownames(res) <- res[["ensembl_gene_id"]]
# # head(rownames(proj)) ####
# feats <- res[match(rownames(proj), res[["ensembl_gene_id"]]), "hgnc_symbol"]
# # head(feats) ####

# counts <- GetAssayData(proj, assay = "RNA", slot = "counts")
# rownames(counts) <- feats
# rownames(counts)[is.na(rownames(counts))] <- rownames(proj)[is.na(rownames(counts))]
# metadata <- proj@meta.data

# proj <- CreateSeuratObject(
#     counts = counts, 
#     project = "kramann", 
#     min.cells = 3, 
#     min.features = 200, 
#     meta.data = metadata
# )
proj <- AddMetaData(cell_type_fine)
proj$cell_type_coarse <- proj$cell_type_original

print(proj) ####
head(proj@meta.data) ####

# print(proj) ####
proj$nCount_RNA <- proj[["nCounts_RNA"]]
proj$nFeature_RNA <- proj[["nFeaturess_RNA"]]

proj[["percent.mt"]] <- PercentageFeatureSet(proj, pattern = "^MT-")

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

proj <- RunHarmony(proj, "sample")

proj[['harmony2']] <- CreateDimReducObject(
  embeddings = proj[['harmony']]@cell.embeddings,
  key = "harmony2_",
  loadings = proj[['pca']]@feature.loadings, 
  assay = "RNA"
)

proj <- FindNeighbors(proj, dims = 1:30, reduction = "harmony")

proj <- RunUMAP(proj, dims = 1:30, return.model = TRUE, reduction = "harmony")

plt <- DimPlot(proj, reduction = "umap", group.by = "cell_type")
ggsave(output_paths[["umap"]], plt, device = "pdf", width = 10, height = 7)

saveRDS(proj, file = output_paths[["project_out"]])

sink(type = "message")
sink()