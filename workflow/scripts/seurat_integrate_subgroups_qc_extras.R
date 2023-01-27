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

plt <- DimPlot(proj, reduction = "umap", group.by = "group") + labs(title="Post-Harmony covariate groups")
ggsave(output_paths[["umap_group_harmony"]], plt, device = "pdf", width = 10, height = 7)

plt <- DimPlot(proj, reduction = "umap", group.by = "region") + labs(title="Post-Harmony heart regions")
ggsave(output_paths[["umap_region_harmony"]], plt, device = "pdf", width = 10, height = 7)

plt <- DimPlot(proj, reduction = "umap", group.by = "status") + labs(title="Post-Harmony health status")
ggsave(output_paths[["umap_status_harmony"]], plt, device = "pdf", width = 10, height = 7)