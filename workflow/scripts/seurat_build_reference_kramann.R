Sys.setenv(CONDA_BUILD_SYSROOT = "/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(biomaRt)

params <- snakemake@params
input_paths <- snakemake@input
output_paths <- snakemake@output
log_paths <- snakemake@log

set.seed(params[["seed"]])

proj <- readRDS(input_paths[["seurat_object"]])
proj@assays$RNA@key <- "rna_"
# print(Key(proj)) ####
proj <- UpdateSeuratObject(proj)
Project(proj) <- "Kramann"

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
res <- getBM(
    attributes=c("ensembl_gene_id","hgnc_symbol"), 
    filters = 'ensembl_gene_id', 
    values = rownames(proj),
    mart = ensembl
)
res <- res[res[["hsapiens_gene_ensembl"]] != "",]
# rownames(res) <- res[["ensembl_gene_id"]]
feats <- res[match(rownames(proj), res[["ensembl_gene_id"]]), "hgnc_symbol"]
head(feats) ####

counts <- GetAssayData(proj, assay = "RNA", slot = "counts")
rownames(counts) <- feats

proj <- CreateSeuratObject(
    counts = counts, 
    project = "kramann", 
    min.cells = 3, 
    min.features = 200, 
    meta.data = metadata
)

print(proj) ####
head(proj@meta.data)
# proj <- DietSeurat(proj)

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

proj <- FindNeighbors(proj, dims = 1:30)

proj <- RunUMAP(proj, dims = 1:30, return.model = TRUE)

plt <- DimPlot(proj, reduction = "umap", group.by = "cell_type")
ggsave(output_paths[["umap"]], plt, device = "pdf")

saveRDS(proj, file = output_paths[["project_out"]])

sink(type = "message")
sink()