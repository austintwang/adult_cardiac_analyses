Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

options(future.globals.maxSize = 8000 * 1024^2)

params = snakemake@params 
input_paths = snakemake@input
output_paths = snakemake@output
log_paths = snakemake@log

set.seed(params[["seed"]])

proj <- readRDS(file = input_paths[["project_rna"]])

ellinor <- readRDS(file = input_paths[["project_ellinor"]])

# print(head(ref@meta.data)) ####
# print(ref) ####
# ref <- subset(x = ref, subset = `Factor.Value.inferred.cell.type...authors.labels.` != "")
# print(ref) ####
# print(proj) ####

anchors <- FindTransferAnchors(
  reference = ellinor,
  query = proj,
  dims = 1:30, 
  reference.reduction = "pca",
  reduction = "pcaproject"
)
# head(proj@meta.data) ####
# head(rownames(proj)) ####

proj_tmp <- MapQuery(
  anchorset = anchors,
  query = proj,
  reference = ellinor,
  refdata = "cell_type",
  reference.reduction = "pcaproject"
)
# head(proj_tmp@meta.data) ####
# head(rownames(proj_tmp)) ####
proj$cell_type_ellinor_fine <- proj_tmp$predicted.id
# head(proj@meta.data) ####
# head(rownames(proj)) ####

proj_tmp <- MapQuery(
  anchorset = anchors,
  query = proj,
  reference = ellinor,
  refdata = "cell_type_leiden06",
  reference.reduction = "pcaproject"
)
# head(proj_tmp@meta.data) ####
# head(rownames(proj_tmp)) ####
proj$cell_type_ellinor_coarse <- proj_tmp$predicted.id
# head(proj@meta.data) ####
# head(rownames(proj)) ####

proj <- FindNeighbors(proj, dims = 1:30, reduction = "pca")
proj <- RunUMAP(proj, dims = 1:30, reduction = "pca")

plt <- DimPlot(proj, reduction = "umap", group.by = "cell_type_ellinor_fine", label = TRUE, width = 10, height = 7)
ggsave(output_paths[["umap_ellinor_fine"]], plt, device = "pdf")

plt <- DimPlot(proj, reduction = "umap", group.by = "cell_type_ellinor_coarse", label = TRUE, width = 10, height = 7)
ggsave(output_paths[["umap_ellinor_coarse"]], plt, device = "pdf")

saveRDS(proj, file = output_paths[["project_out"]])