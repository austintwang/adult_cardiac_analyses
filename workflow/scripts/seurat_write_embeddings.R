Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(Seurat)

stallion <- c("1"="#D51F26","2"="#272E6A","3"="#208A42","4"="#89288F","5"="#F47D2B", "6"="#FEE500","7"="#8A9FD1","8"="#C06CAB","19"="#E6C2DC",
               "10"="#90D5E4", "11"="#89C75F","12"="#F37B7D","13"="#9983BD","14"="#D24B27","15"="#3BBCA8", "16"="#6E4B9E","17"="#0C727C", "18"="#7E1416","9"="#D8A767","20"="#3D3D3D")
stallion <- unname(stallion)

plot_fn <- function(object, group, reduction, colors) {
    cats <- sort(unique(object@meta.data[[group]]))
    colors_out <- rep_len(colors, length(cats))
    names(colors_out) <- cats
    DimPlot(object, reduction = reduction, group.by = group, label = TRUE, cols = colors_out, pt.size=0.1)
}


params = snakemake@params 
input_paths = snakemake@input
output_paths = snakemake@output
log_paths = snakemake@log

set.seed(params[["seed"]])

proj <- readRDS(file = input_paths[["project_in"]])

DefaultAssay(object = proj) <- "RNA"

write.table(proj@meta.data, file=output_paths[["metadata"]], quote=FALSE, sep='\t', col.names = NA)

rna_pca <- Embeddings(proj, reduction = "pca")
write.table(rna_pca, file=output_paths[["rna_pca"]], quote=FALSE, sep='\t', col.names = NA)

rna_harmony <- Embeddings(proj, reduction = "harmony")
write.table(rna_harmony, file=output_paths[["rna_harmony"]], quote=FALSE, sep='\t', col.names = NA)

rna_umap <- Embeddings(proj, reduction = "umap")
write.table(rna_umap, file=output_paths[["rna_umap"]], quote=FALSE, sep='\t', col.names = NA)

atac_lsi <- Embeddings(proj, reduction = "atac")
write.table(atac_lsi, file=output_paths[["atac_lsi"]], quote=FALSE, sep='\t', col.names = NA)

atac_harmony <- Embeddings(proj, reduction = "atac_harmony")
write.table(atac_harmony, file=output_paths[["atac_harmony"]], quote=FALSE, sep='\t', col.names = NA)

atac_umap <- Embeddings(proj, reduction = "umap_atac")
write.table(atac_umap, file=output_paths[["atac_umap"]], quote=FALSE, sep='\t', col.names = NA)

plt <- plot_fn(proj, "dataset", "umap", stallion) + labs(title="All counts embedding: datasets")
ggsave(output_paths[["umap_rna_dataset"]], plt, device = "pdf", width = 10, height = 7)

plt <- plot_fn(proj, "cell_types_l1", "umap", stallion) + labs(title="All counts embedding: Level-1 cell types")
ggsave(output_paths[["umap_rna_cell_types"]], plt, device = "pdf", width = 10, height = 7)

plt <- plot_fn(proj, "dataset", "umap_atac", stallion) + labs(title="All counts embedding: datasets")
ggsave(output_paths[["umap_atac_dataset"]], plt, device = "pdf", width = 10, height = 7)

plt <- plot_fn(proj, "cell_types_l1", "umap_atac", stallion) + labs(title="All counts embedding: Level-1 cell types")
ggsave(output_paths[["umap_atac_cell_types"]], plt, device = "pdf", width = 10, height = 7)

sink(type = "message")
sink()