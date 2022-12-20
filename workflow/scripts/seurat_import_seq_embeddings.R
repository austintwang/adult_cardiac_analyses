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

cellnames <- unlist(readLines(input_paths[["cells"]]))

mat_scbasset <- as.matrix(readMM(input_paths[["scbasset"]]))
rownames(mat_scbasset) <- cellnames

mat_cellspace <- as.matrix(readMM(input_paths[["cellspace"]]))
rownames(mat_cellspace) <- cellnames

head(mat_cellspace) ####

proj <- readRDS(file = input_paths[["project_in"]])

idx <- cellnames %in% Cells(proj)
head(idx) ####
mat_scbasset <- mat_scbasset[idx,,drop=FALSE]
mat_cellspace <- mat_cellspace[idx,,drop=FALSE]

head(mat_cellspace) ####

proj <- subset(proj, cells = cellnames)

a <- CreateDimReducObject(embeddings = mat_scbasset, key = "SCB_", assay = "RNA_train")
print(a) ####
proj[["scbasset"]] <- a
b <- CreateDimReducObject(embeddings = mat_cellspace, key = "CSP_", assay = "RNA_train")
proj[["cellspace"]] <- b

proj <- RunHarmony(proj, "dataset", assay.use = "RNA_train", reduction = "scbasset", reduction.save = "scbasset_harmony")
proj <- RunHarmony(proj, "dataset", assay.use = "RNA_train", reduction = "cellspace", reduction.save = "cellspace_harmony")

scbasset_int_data <- Embeddings(proj, reduction = "scbasset_harmony")
cellspace_int_data <- Embeddings(proj, reduction = "cellspace_harmony")

proj <- FindNeighbors(proj, dims = 1:30, reduction = "scbasset_harmony", graph.name = "nn_scbasset")
proj <- RunUMAP(proj, dims = 1:30, reduction = "scbasset_harmony", nn.name = "nn_scbasset", reduction.name = "umap_scbasset")

proj <- FindNeighbors(proj, dims = 1:30, reduction = "cellspace_harmony", graph.name = "nn_cellspace")
proj <- RunUMAP(proj, dims = 1:30, reduction = "cellspace_harmony", nn.name = "nn_cellspace", reduction.name = "umap_cellspace")

plt <- plot_fn(proj, "cell_types_l1", "umap_cellspace", stallion) + labs(title="Cell Types (CellSpace)")
ggsave(output_paths[["umap_cellspace"]], plt, device = "pdf", width = 10, height = 7)

plt <- plot_fn(proj, "cell_types_l1", "umap_scbasset", stallion) + labs(title="Cell Types (scBasset)")
ggsave(output_paths[["umap_scbasset"]], plt, device = "pdf", width = 10, height = 7)

proj@meta.data$seq_emb_bc <- paste0(proj@meta.data$dataset, "_", proj@meta.data$barcode_atac)
write.table(proj@meta.data, file=output_paths[["metadata"]], quote=FALSE, sep='\t', col.names = NA)

writeMM(Matrix(scbasset_int_data, sparse = TRUE), output_paths[["scbasset_harmony_mat"]])
writeLines(rownames(scbasset_int_data), con = output_paths[["scbasset_harmony_rows"]])
writeLines(colnames(scbasset_int_data), con = output_paths[["scbasset_harmony_cols"]])

writeMM(Matrix(cellspace_int_data, sparse = TRUE), output_paths[["cellspace_harmony_mat"]])
writeLines(rownames(cellspace_int_data), con = output_paths[["cellspace_harmony_rows"]])
writeLines(colnames(cellspace_int_data), con = output_paths[["cellspace_harmony_cols"]])

# plt <- plot_fn(proj, "bc_counts", "umap", stallion) + labs(title="Barcode counts")
# ggsave(output_paths[["umap"]], plt, device = "pdf")

saveRDS(proj, file = output_paths[["project_out"]])