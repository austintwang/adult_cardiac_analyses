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

teichmann <- readRDS(file = input_paths[["project_teichmann"]])

# print(head(ref@meta.data)) ####
# print(ref) ####
# ref <- subset(x = ref, subset = `Factor.Value.inferred.cell.type...authors.labels.` != "")
# print(ref) ####
# print(proj) ####

anchors <- FindTransferAnchors(
  reference = teichmann,
  query = proj,
  dims = 1:30, 
  reference.reduction = "pca",
  reduction = "pcaproject"
)
# head(proj@meta.data) ####
# head(rownames(proj)) ####

# proj_tmp <- MapQuery(
#   anchorset = anchors,
#   query = proj,
#   reference = ellinor,
#   refdata = "cell_type",
#   reference.reduction = "pcaproject"
# )
predictions <- TransferData(
  anchorset = teichmann,
  refdata = proj$cell_type_fine
)
# head(proj_tmp@meta.data) ####
# head(rownames(proj_tmp)) ####
proj$cell_type_teichmann_fine <- predictions$predicted.id
# head(proj@meta.data) ####
# head(rownames(proj)) ####

# proj_tmp <- MapQuery(
#   anchorset = anchors,
#   query = proj,
#   reference = ellinor,
#   refdata = "cell_type_leiden06",
#   reference.reduction = "pcaproject"
# )
predictions <- TransferData(
  anchorset = teichmann,
  refdata = proj$cell_type_coarse
)
# head(proj_tmp@meta.data) ####
# head(rownames(proj_tmp)) ####
proj$cell_type_teichmann_coarse <- predictions$predicted.id
# head(proj@meta.data) ####
# head(rownames(proj)) ####

plt <- DimPlot(proj, reduction = "umap", group.by = "cell_type_ellinor_fine", label = TRUE)
ggsave(output_paths[["umap_teichmann_fine"]], plt, device = "pdf", width = 10, height = 7)

plt <- DimPlot(proj, reduction = "umap", group.by = "cell_type_ellinor_coarse", label = TRUE)
ggsave(output_paths[["umap_teichmann_coarse"]], plt, device = "pdf", width = 10, height = 7)

label_data <- proj@meta.data[,c("cell_type_teichmann_coarse", "cell_type_teichmann_fine"),drop=FALSE]
write.table(label_data, file= output_paths[["data_out"]], quote=FALSE, sep='\t', col.names = NA)