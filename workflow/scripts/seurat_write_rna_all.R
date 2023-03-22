Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(Matrix)
library(harmony)

params = snakemake@params 
input_paths = snakemake@input
output_paths = snakemake@output
log_paths = snakemake@log

set.seed(params[["seed"]])

stallion <- c("1"="#D51F26","2"="#272E6A","3"="#208A42","4"="#89288F","5"="#F47D2B", "6"="#FEE500","7"="#8A9FD1","8"="#C06CAB","19"="#E6C2DC",
               "10"="#90D5E4", "11"="#89C75F","12"="#F37B7D","13"="#9983BD","14"="#D24B27","15"="#3BBCA8", "16"="#6E4B9E","17"="#0C727C", "18"="#7E1416","9"="#D8A767","20"="#3D3D3D")
stallion <- unname(stallion)

plot_fn <- function(object, group, reduction, colors) {
    cats <- sort(unique(object@meta.data[[group]]))
    colors_out <- rep_len(colors, length(cats))
    names(colors_out) <- cats
    DimPlot(object, reduction = reduction, group.by = group, label = TRUE, cols = colors_out, pt.size=0.1)
}

proj <- readRDS(file = input_paths[["project_in"]])

DefaultAssay(object = proj) <- "RNA"

rna_harmony <- Embeddings(proj, reduction = "harmony")
rna_pca <- Embeddings(proj, reduction = "pca")

proj@meta.data$seq_emb_bc <- paste0(proj@meta.data$dataset, "_", proj@meta.data$barcode_atac)
write.table(proj@meta.data, file=output_paths[["metadata"]], quote=FALSE, sep='\t', col.names = NA)

writeMM(Matrix(rna_harmony, sparse = TRUE), output_paths[["rna_harmony_mat"]])
writeLines(rownames(rna_harmony), con = output_paths[["rna_harmony_rows"]])
writeLines(colnames(rna_harmony), con = output_paths[["rna_harmony_cols"]])

writeMM(Matrix(rna_pca, sparse = TRUE), output_paths[["rna_pca_mat"]])
writeLines(rownames(rna_pca), con = output_paths[["rna_pca_rows"]])
writeLines(colnames(rna_pca), con = output_paths[["rna_pca_cols"]])