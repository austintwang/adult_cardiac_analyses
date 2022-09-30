Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(harmony)

params = snakemake@params 
input_paths = snakemake@input
output_paths = snakemake@output
log_paths = snakemake@log

set.seed(params[["seed"]])

read_fn <- function(path) {
    sub_proj <- readRDS(file = path)
    sub_proj@meta.data[, "cell_types_2", drop = FALSE]
}

plot_fn <- function(object, group, reduction, colors) {
    cats <- sort(unique(object@meta.data[[group]]))
    colors_out <- rep_len(colors, length(cats))
    names(colors_out) <- cats
    DimPlot(proj_merged, reduction = reduction, group.by = group, label = TRUE, cols = colors_out, pt.size=0.1, label.size=2)
}

tables <- lapply(input_paths[["projects_subcluster"]], read_fn)
head(tables[[1]])
sub_data <- do.call(rbind, tables)

proj <- readRDS(file = input_paths[["project_integrated"]])

label_data <- proj@meta.data[, "cell_types_1", drop = FALSE]
label_data[rownames(sub_data),] <- sub_data
head(label_data) ####

proj_merged <- AddMetaData(proj, label_data)

stallion <- c("1"="#D51F26","2"="#272E6A","3"="#208A42","4"="#89288F","5"="#F47D2B", "6"="#FEE500","7"="#8A9FD1","8"="#C06CAB","19"="#E6C2DC",
               "10"="#90D5E4", "11"="#89C75F","12"="#F37B7D","13"="#9983BD","14"="#D24B27","15"="#3BBCA8", "16"="#6E4B9E","17"="#0C727C", "18"="#7E1416","9"="#D8A767","20"="#3D3D3D")
stallion <- unname(stallion)

plt <- plot_fn(proj_merged, "cell_types_2", "umap", stallion) + labs(title="Cell Types (Round 2)")
ggsave(output_paths[["umap"]], plt, device = "pdf", width = 13, height = 7)

saveRDS(proj_merged, file = output_paths[["project_out"]])

sink(type = "message")
sink()