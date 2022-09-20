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
   read.table(file = path, sep = '\t', header = TRUE, row.names = 1)
}

plot_fn <- function(object, group, reduction, colors) {
    cats <- sort(unique(object@meta.data[[group]]))
    colors_out <- rep_len(colors, length(cats))
    names(colors_out) <- cats
    print(colors_out) ####
    DimPlot(proj_merged, reduction = reduction, group.by = group, label = TRUE, cols = colors_out, pt.size=0.1)
}

tables <- lapply(input_paths[["tables"]], read_fn)
head(tables[[1]])
label_data <- do.call(cbind, tables)
# head(label_data) ####

proj <- readRDS(file = input_paths[["project_integrated"]])
proj_merged <- AddMetaData(proj, label_data)


stallion <- c("1"="#D51F26","2"="#272E6A","3"="#208A42","4"="#89288F","5"="#F47D2B", "6"="#FEE500","7"="#8A9FD1","8"="#C06CAB","19"="#E6C2DC",
               "10"="#90D5E4", "11"="#89C75F","12"="#F37B7D","13"="#9983BD","14"="#D24B27","15"="#3BBCA8", "16"="#6E4B9E","17"="#0C727C", "18"="#7E1416","9"="#D8A767","20"="#3D3D3D")
stallion <- unname(stallion)

# plt <- DimPlot(proj_merged, reduction = "umap", group.by = "cell_type_kramann", label = TRUE) + labs(title="Kramann labels")
plt <- plot_fn(proj_merged, "cell_type_kramann", "umap", stallion) + labs(title="Kramann labels")
ggsave(output_paths[["umap_kramann"]], plt, device = "pdf", width = 10, height = 7)

plt <-  plot_fn(proj_merged, "cell_type_ellinor_coarse", "umap", stallion) + labs(title="Ellinor coarse labels")
# plt <- DimPlot(proj_merged, reduction = "umap", group.by = "cell_type_ellinor_coarse", label = TRUE) + labs(title="Ellinor coarse labels")
ggsave(output_paths[["umap_ellinor_coarse"]], plt, device = "pdf", width = 10, height = 7)

plt <- plot_fn(proj_merged, "cell_type_ellinor_fine", "umap", stallion) + labs(title="Ellinor fine labels")
# plt <- DimPlot(proj_merged, reduction = "umap", group.by = "cell_type_ellinor_fine", label = TRUE) + labs(title="Ellinor fine labels")
ggsave(output_paths[["umap_ellinor_fine"]], plt, device = "pdf", width = 10, height = 7)

plt <- plot_fn(proj_merged, "cell_type_teichmann_coarse", "umap", stallion) + labs(title="Teichmann coarse labels")
# plt <- DimPlot(proj_merged, reduction = "umap", group.by = "cell_type_teichmann_coarse", label = TRUE) + labs(title="Teichmann coarse labels")
ggsave(output_paths[["umap_teichmann_coarse"]], plt, device = "pdf", width = 10, height = 7)

plt <- plot_fn(proj_merged, "cell_type_teichmann_fine", "umap", stallion) + labs(title="Teichmann fine labels")
# plt <- DimPlot(proj_merged, reduction = "umap", group.by = "cell_type_teichmann_fine", label = TRUE) + labs(title="Teichmann fine labels")
ggsave(output_paths[["umap_teichmann_fine"]], plt, device = "pdf", width = 10, height = 7)

plt <- plot_fn(proj_merged, "cell_type_azimuth_l1", "umap", stallion) + labs(title="Azimuth coarse labels")
ggsave(output_paths[["umap_azimuth_coarse"]], plt, device = "pdf", width = 10, height = 7)

plt <- plot_fn(proj_merged, "cell_type_azimuth_l2", "umap", stallion) + labs(title="Azimuth fine labels")
ggsave(output_paths[["umap_azimuth_fine"]], plt, device = "pdf", width = 10, height = 7)


saveRDS(proj_merged, file = output_paths[["project_out"]])

sink(type = "message")
sink()