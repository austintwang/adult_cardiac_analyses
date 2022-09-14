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

kramann <- readRDS(file = input_paths[["project_kramann"]])
# head(rownames(kramann)) ####

anchors <- FindTransferAnchors(
  reference = kramann,
  query = proj,
  dims = 1:30, 
  reference.reduction  = "harmony2",
  reduction = "pcaproject"
)
proj_tmp <- MapQuery(
  anchorset = anchors,
  query = proj,
  reference = kramann,
  refdata = "cell_type",
  reference.reduction = "pcaproject"
)
proj$cell_type_kramann <- proj_tmp$predicted.id

proj <- FindNeighbors(proj, dims = 1:30, reduction = "pca")
proj <- RunUMAP(proj, dims = 1:30, reduction = "pca")

plt <- DimPlot(proj, reduction = "umap", group.by = "cell_type_kramann", label = TRUE)
ggsave(output_paths[["umap_kramann"]], plt, device = "pdf", width = 10, height = 7)

saveRDS(proj, file = output_paths[["project_out"]])