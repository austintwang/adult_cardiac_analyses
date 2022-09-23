Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(plyr)
library(pheatmap)

confusionMatrix <- function(
  i = NULL, 
  j = NULL
  ){
  ui <- unique(i)
  uj <- unique(j)
  m <- Matrix::sparseMatrix(
    i = match(i, ui),
    j = match(j, uj),
    x = rep(1, length(i)),
    dims = c(length(ui), length(uj))
  )
  rownames(m) <- ui
  colnames(m) <- uj
  m
}

plot_cm <- function(clusters, ref_labels) {
  proj <- readRDS(file = input_paths[["project_in"]])

  cM <- confusionMatrix(proj$seurat_clusters, proj$cell_type_ref)
  cM <- as.matrix(cM)
  cM <- cbind(cM, Unknown = 0.01)
  labelOld <- rownames(cM)
  labelNew <- colnames(cM)[apply(cM, 1, which.max)]
  proj$labels_named <- plyr::mapvalues(x = proj$seurat_clusters, from = labelOld, to = labelNew)

  plt <- pheatmap::pheatmap(
      mat = cM / rowSums(cM), 
      border_color = "black"
  )
  plt
}

params = snakemake@params 
input_paths = snakemake@input
output_paths = snakemake@output
log_paths = snakemake@log

set.seed(params[["seed"]])

proj <- readRDS(file = input_paths[["project_in"]])

plt <- plot_cm(proj$seurat_clusters, proj$cell_type_kramann_coarse) + labs(title="Kramann coarse labels")
ggsave(output_paths[["mat_kramann_coarse"]], plt, device = "pdf", width = 10, height = 7)

plt <- plot_cm(proj$seurat_clusters, proj$cell_type_kramann_fine) + labs(title="Kramann fine labels")
ggsave(output_paths[["mat_kramann_fine"]], plt, device = "pdf", width = 10, height = 7)

plt <- plot_cm(proj$seurat_clusters, proj$cell_type_ellinor_coarse) + labs(title="Ellinor coarse labels")
ggsave(output_paths[["mat_ellinor_coarse"]], plt, device = "pdf", width = 10, height = 7)

plt <- plot_cm(proj$seurat_clusters, proj$cell_type_ellinor_fine) + labs(title="Ellinor fine labels")
ggsave(output_paths[["mat_ellinor_fine"]], plt, device = "pdf", width = 10, height = 7)

plt <- plot_cm(proj$seurat_clusters, proj$cell_type_teichmann_coarse) + labs(title="Teichmann coarse labels")
ggsave(output_paths[["mat_teichmann_coarse"]], plt, device = "pdf", width = 10, height = 7)

plt <- plot_cm(proj$seurat_clusters, proj$cell_type_teichmann_fine) + labs(title="Teichmann fine labels")
ggsave(output_paths[["mat_teichmann_fine"]], plt, device = "pdf", width = 10, height = 7)

plt <- plot_cm(proj$seurat_clusters, proj$cell_type_azimuth_coarse) + labs(title="Azimuth coarse labels")
ggsave(output_paths[["mat_azimuth_coarse"]], plt, device = "pdf", width = 10, height = 7)

plt <- plot_cm(proj$seurat_clusters, proj$cell_type_azimuth_fine) + labs(title="Azimuth fine labels")
ggsave(output_paths[["mat_azimuth_fine"]], plt, device = "pdf", width = 10, height = 7)

