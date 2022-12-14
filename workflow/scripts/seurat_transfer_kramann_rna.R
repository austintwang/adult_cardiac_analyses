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
  reference.reduction  = "pca",
  reduction = "pcaproject"
)
# proj_tmp <- MapQuery(
#   anchorset = anchors,
#   query = proj,
#   reference = kramann,
#   refdata = "cell_type",
#   reference.reduction = "pcaproject"
# )
predictions <- TransferData(
  anchorset = anchors,
  refdata = kramann$cell_type_coarse
)
proj$cell_type_kramann_coarse <- predictions$predicted.id

predictions <- TransferData(
  anchorset = anchors,
  refdata = kramann$cell_type_fine
)
proj$cell_type_kramann_fine <- predictions$predicted.id

plt <- DimPlot(proj, reduction = "umap", group.by = "cell_type_kramann_fine", label = TRUE)
ggsave(output_paths[["umap_kramann_fine"]], plt, device = "pdf", width = 10, height = 7)

plt <- DimPlot(proj, reduction = "umap", group.by = "cell_type_kramann_coarse", label = TRUE)
ggsave(output_paths[["umap_kramann_coarse"]], plt, device = "pdf", width = 10, height = 7)

label_data <- proj@meta.data[,c("cell_type_kramann_coarse", "cell_type_kramann_fine"),drop=FALSE]
write.table(label_data, file= output_paths[["data_out"]], quote=FALSE, sep='\t', col.names = NA)
