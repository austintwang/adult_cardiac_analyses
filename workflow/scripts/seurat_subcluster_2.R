Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

params = snakemake@params 
wildcards = snakemake@wildcards
input_paths = snakemake@input
output_paths = snakemake@output
log_paths = snakemake@log

set.seed(params[["seed"]])

stallion <- c("1"="#D51F26","2"="#272E6A","3"="#208A42","4"="#89288F","5"="#F47D2B", "6"="#FEE500","7"="#8A9FD1","8"="#C06CAB","19"="#E6C2DC",
               "10"="#90D5E4", "11"="#89C75F","12"="#F37B7D","13"="#9983BD","14"="#D24B27","15"="#3BBCA8", "16"="#6E4B9E","17"="#0C727C", "18"="#7E1416","9"="#D8A767","20"="#3D3D3D")
stallion <- unname(stallion)

plot_fn <- function(object, group, reduction, colors) {
    group_data <- as.data.frame(table(object@meta.data[[group]]))
    cats <- group_data[order(group_data$Freq, decreasing = TRUE), 1]

    # cats <- sort(unique(object@meta.data[[group]]))
    colors_out <- rep_len(colors, length(cats))
    names(colors_out) <- cats
    DimPlot(object, reduction = reduction, group.by = group, label = TRUE, cols = colors_out, pt.size=0.2)
}

proj <- readRDS(file = input_paths[["project_in"]])

print(proj) ####

proj <- FindNeighbors(proj, dims = 1:30, reduction = "harmony")
proj <- RunUMAP(proj, dims = 1:30, reduction = "harmony")
proj <- FindClusters(object = proj, resolution = 0.5) 

plt <- plot_fn(proj, "seurat_clusters", "umap", stallion) + labs(title="Seurat subclustering, Train set")
ggsave(output_paths[["umap"]], plt, device = "pdf")

plt <- plot_fn(proj, "seurat_clusters", "umap_test", stallion) + labs(title="Seurat subclustering, Test Set")
ggsave(output_paths[["umap_test"]], plt, device = "pdf", width = 10, height = 7)

plt <- DimPlot(proj, reduction = "umap", group.by = "group") + labs(title="Post-Harmony covariate groups")
ggsave(output_paths[["umap_group_harmony"]], plt, device = "pdf", width = 10, height = 7)

plt <- DimPlot(proj, reduction = "umap", group.by = "region") + labs(title="Post-Harmony heart regions")
ggsave(output_paths[["umap_region_harmony"]], plt, device = "pdf", width = 10, height = 7)

plt <- DimPlot(proj, reduction = "umap", group.by = "status") + labs(title="Post-Harmony health status")
ggsave(output_paths[["umap_status_harmony"]], plt, device = "pdf", width = 10, height = 7)

plt <- plot_fn(proj, "cell_type_kramann_coarse", "umap", stallion) + labs(title="Kramann coarse labels")
ggsave(output_paths[["umap_kramann_coarse"]], plt, device = "pdf", width = 10, height = 7)

plt <- plot_fn(proj, "cell_type_kramann_fine", "umap", stallion) + labs(title="Kramann fine labels")
ggsave(output_paths[["umap_kramann_fine"]], plt, device = "pdf", width = 10, height = 7)

plt <-  plot_fn(proj, "cell_type_ellinor_coarse", "umap", stallion) + labs(title="Ellinor coarse labels")
# plt <- DimPlot(proj, reduction = "umap", group.by = "cell_type_ellinor_coarse", label = TRUE) + labs(title="Ellinor coarse labels")
ggsave(output_paths[["umap_ellinor_coarse"]], plt, device = "pdf", width = 10, height = 7)

plt <- plot_fn(proj, "cell_type_ellinor_fine", "umap", stallion) + labs(title="Ellinor fine labels")
# plt <- DimPlot(proj, reduction = "umap", group.by = "cell_type_ellinor_fine", label = TRUE) + labs(title="Ellinor fine labels")
ggsave(output_paths[["umap_ellinor_fine"]], plt, device = "pdf", width = 10, height = 7)

plt <- plot_fn(proj, "cell_type_teichmann_coarse", "umap", stallion) + labs(title="Teichmann coarse labels")
# plt <- DimPlot(proj, reduction = "umap", group.by = "cell_type_teichmann_coarse", label = TRUE) + labs(title="Teichmann coarse labels")
ggsave(output_paths[["umap_teichmann_coarse"]], plt, device = "pdf", width = 10, height = 7)

plt <- plot_fn(proj, "cell_type_teichmann_fine", "umap", stallion) + labs(title="Teichmann fine labels")
# plt <- DimPlot(proj, reduction = "umap", group.by = "cell_type_teichmann_fine", label = TRUE) + labs(title="Teichmann fine labels")
ggsave(output_paths[["umap_teichmann_fine"]], plt, device = "pdf", width = 10, height = 7)

plt <- plot_fn(proj, "cell_type_azimuth_coarse", "umap", stallion) + labs(title="Azimuth coarse labels")
ggsave(output_paths[["umap_azimuth_coarse"]], plt, device = "pdf", width = 10, height = 7)

plt <- plot_fn(proj, "cell_type_azimuth_fine", "umap", stallion) + labs(title="Azimuth fine labels")
ggsave(output_paths[["umap_azimuth_fine"]], plt, device = "pdf", width = 10, height = 7)

saveRDS(proj, file = output_paths[["project_out"]])