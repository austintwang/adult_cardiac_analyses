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

cellnames_atac <- unlist(readLines(input_paths[["cells"]]))

mat_atac <- as.matrix(readMM(input_paths[["atac"]]))
# rownames(mat_atac) <- cellnames_atac

proj <- readRDS(file = input_paths[["project_in"]])

proj$cellnames_atac <- paste0(proj$dataset, "#", proj$barcode_atac)

idx_atac <- proj$cellnames_atac %in% cellnames_atac
cellnames_subset <- Cells(proj)[idx_atac]
proj <- subset(proj, cells = cellnames_subset)

idx <- na.omit(match(proj$cellnames_atac, cellnames_atac))
mat_atac <- mat_atac[idx,,drop=FALSE]
cellnames <- Cells(proj)[idx]
rownames(mat_atac) <- cellnames

a <- CreateDimReducObject(embeddings = mat_atac, key = "LSI_", assay = "RNA")
print(a) ####
proj[["atac"]] <- a


proj <- RunHarmony(proj, "dataset", assay.use = "RNA", reduction = "atac", reduction.save = "atac_harmony")

proj <- FindNeighbors(proj, dims = 1:30, reduction = "atac_harmony", graph.name = "nn_atac")
proj <- RunUMAP(proj, dims = 1:30, reduction = "atac_harmony", nn.name = "nn_atac", reduction.name = "umap_atac")


saveRDS(proj, file = output_paths[["project_out"]])