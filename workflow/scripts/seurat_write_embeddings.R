Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(Seurat)

params = snakemake@params 
input_paths = snakemake@input
output_paths = snakemake@output
log_paths = snakemake@log

set.seed(params[["seed"]])

proj <- readRDS(file = input_paths[["project_in"]])

DefaultAssay(object = proj) <- "RNA"

rna_pca <- Embeddings(proj, reduction = "pca")
write.table(rna_pca, file=output_paths[["rna_pca"]], quote=FALSE, sep='\t', col.names = NA)

rna_harmony <- Embeddings(proj, reduction = "harmony")
write.table(rna_harmony, file=output_paths[["rna_harmony"]], quote=FALSE, sep='\t', col.names = NA)

rna_umap <- Embeddings(proj, reduction = "umap")
write.table(rna_umap, file=output_paths[["rna_umap"]], quote=FALSE, sep='\t', col.names = NA)

sink(type = "message")
sink()