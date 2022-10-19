Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(Seurat)
library(ggplot2)
library(patchwork)

params = snakemake@params 
input_paths = snakemake@input
output_paths = snakemake@output
log_paths = snakemake@log

set.seed(params[["seed"]])

proj <- readRDS(file = input_paths[["project_in"]])
head(proj@meta.data) ####

proj$log_counts <- log10(proj$nCount_RNA)
plt <- FeaturePlot(proj, reduction = "umap", features = "log_counts") + labs(title="Post-Harmony log10 UMIs")
ggsave(output_paths[["umap_counts"]], plt, device = "pdf", width = 10, height = 7)

doubletfinder_cols <- grep("pANN", names(proj@meta.data), value = TRUE)
d <- proj@meta.data[, doubletfinder_cols, drop = FALSE ]
d[is.na(d)] <- 0
d$sum <- rowSums(d)
proj$pANN <- d$sum

plt <- FeaturePlot(proj, reduction = "umap", features = "pANN") +  labs(title = "Post-Harmony Doubletfinder pANN")
ggsave(output_paths[["umap_doubletfinder"]], plt, device = "pdf", width = 10, height = 7)

proj$amulet_nlp <- -log10(proj$amulet_pval)
plt <- FeaturePlot(proj, reduction = "umap", features = "amulet_nlp") + labs(title="Post-Harmony Amulet -log10 p-Value")
ggsave(output_paths[["umap_amulet"]], plt, device = "pdf", width = 10, height = 7)

sink(type = "message")
sink()