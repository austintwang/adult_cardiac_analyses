Sys.setenv(CONDA_BUILD_SYSROOT = "/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

params <- snakemake@params
input_paths <- snakemake@input
output_paths <- snakemake@output
log_paths <- snakemake@log

library(SeuratData, lib.loc=input_paths[["azimuth_library_dir"]])
library(SeuratDisk, lib.loc=input_paths[["seuratdisk_library_dir"]])

set.seed(params[["seed"]])

proj <- LoadH5Seurat(input_paths[["h5seurat"]])

print(proj) ####
head(proj@meta.data) ####
# proj <- DietSeurat(proj)

# print(proj) ####
proj$nCount_RNA <- proj[["n_counts"]]
proj$nFeature_RNA <- proj[["n_genes"]]

proj[["percent.mt"]] <- PercentageFeatureSet(proj, pattern = "^MT-")
proj$cell_type_coarse <- proj[["cell_type"]]
proj$cell_type_fine <- proj[["cell_states"]]
proj$cell_type_fine[is.na(proj$cell_type_fine)] <- proj$cell_type_coarse[is.na(proj$cell_type_fine)]
proj$cell_type_fine[proj$cell_type_fine == "nan"] <- proj$cell_type_coarse[proj$cell_type_fine == "nan"]
proj <- subset(proj, subset = cell_type_coarse != "NotAssigned")

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

plt <- DimPlot(proj, reduction = "umap", group.by = "cell_type_coarse")
ggsave(output_paths[["umap_coarse"]], plt, device = "pdf", width = 10, height = 7)
plt <- DimPlot(proj, reduction = "umap", group.by = "cell_type_fine")
ggsave(output_paths[["umap_fine"]], plt, device = "pdf", width = 10, height = 7)

saveRDS(proj, file = output_paths[["project_out"]])

sink(type = "message")
sink()