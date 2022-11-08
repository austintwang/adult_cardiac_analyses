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

proj$log_frags <- log10(proj$frag_count)
plt <- FeaturePlot(proj, reduction = "umap", features = "log_frags") + labs(title="Post-Harmony log10 Fragment Count")
ggsave(output_paths[["umap_frag"]], plt, device = "pdf", width = 10, height = 7)

proj$log_ratio <- proj$log_counts - proj$log_frags
plt <- FeaturePlot(proj, reduction = "umap", features = "log_ratio") + labs(title="Post-Harmony log10 RNA-ATAC Signal Ratio")
ggsave(output_paths[["umap_ratio"]], plt, device = "pdf", width = 10, height = 7)

plt <- FeaturePlot(proj, reduction = "umap", features = "percent.mt") + labs(title="Post-Harmony RNA Mito Frac")
ggsave(output_paths[["umap_mito"]], plt, device = "pdf", width = 10, height = 7)

plt <- FeaturePlot(proj, reduction = "umap", features = "tss_enr") + labs(title="Post-Harmony ATAC TSS Enrichment")
ggsave(output_paths[["umap_tss"]], plt, device = "pdf", width = 10, height = 7)

proj$log_frag_tss <- proj$log_frags + log10(proj$tss_enr)
plt <- FeaturePlot(proj, reduction = "umap", features = "log_frag_tss") + labs(title="Post-Harmony Total TSS Fragment Score")
ggsave(output_paths[["umap_frag_tss"]], plt, device = "pdf", width = 10, height = 7)

proj$tss_score <- proj$log_counts - proj$log_frag_tss
plt <- FeaturePlot(proj, reduction = "umap", features = "tss_score") + labs(title="Post-Harmony RNA-ATAC TSS Score")
ggsave(output_paths[["umap_score"]], plt, device = "pdf", width = 10, height = 7)

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