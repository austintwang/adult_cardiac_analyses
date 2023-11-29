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

metadata_read <- function(table_path) {
    metadata <- read.table(file = table_path, sep = '\t', header = TRUE)
    rownames(metadata) <- metadata$cell_id_seurat
    metadata[, "num_overlaps", drop = FALSE]
}

tables <- lapply(input_paths[["atac_qc_extended"]], metadata_read)
overlap_data <- do.call("rbind", tables) 

proj <- readRDS(file = input_paths[["project_in"]])
proj <- AddMetaData(proj, overlap_data)
head(proj@meta.data) ####

proj$log_counts <- log10(proj$nCount_RNA)
plt <- FeaturePlot(proj, reduction = "umap", features = "log_counts") + labs(title="Post-Harmony log10 UMIs")
ggsave(output_paths[["umap_counts"]], plt, device = "pdf", width = 10, height = 7)

proj$log_frags <- log10(proj$frag_count)
plt <- FeaturePlot(proj, reduction = "umap", features = "log_frags") + labs(title="Post-Harmony log10 Fragment Count")
ggsave(output_paths[["umap_frag"]], plt, device = "pdf", width = 10, height = 7)

# proj$log_ratio <- proj$log_counts - proj$log_frags
# plt <- FeaturePlot(proj, reduction = "umap", features = "log_ratio") + labs(title="Post-Harmony log10 RNA-ATAC Signal Ratio")
# ggsave(output_paths[["umap_ratio"]], plt, device = "pdf", width = 10, height = 7)

plt <- FeaturePlot(proj, reduction = "umap", features = "percent.mt") + labs(title="Post-Harmony RNA Mito Frac")
ggsave(output_paths[["umap_mito"]], plt, device = "pdf", width = 10, height = 7)

plt <- FeaturePlot(proj, reduction = "umap", features = "tss_enr") + labs(title="Post-Harmony ATAC TSS Enrichment")
ggsave(output_paths[["umap_tss"]], plt, device = "pdf", width = 10, height = 7)

# proj$log_frag_tss <- proj$log_frags + log10(proj$tss_enr)
# plt <- FeaturePlot(proj, reduction = "umap", features = "log_frag_tss") + labs(title="Post-Harmony Total TSS Fragment Score")
# ggsave(output_paths[["umap_frag_tss"]], plt, device = "pdf", width = 10, height = 7)

# proj$tss_score <- proj$log_counts - proj$log_frag_tss
# plt <- FeaturePlot(proj, reduction = "umap", features = "tss_score") + labs(title="Post-Harmony RNA-ATAC TSS Score")
# ggsave(output_paths[["umap_score"]], plt, device = "pdf", width = 10, height = 7)

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

plt <- FeaturePlot(proj, reduction = "umap", features = "num_overlaps") + labs(title="Post-Harmony AMULET overlap count")
ggsave(output_paths[["umap_overlap_count"]], plt, device = "pdf", width = 10, height = 7)

proj$frac_overlap <- proj$num_overlaps / proj$frag_count
plt <- FeaturePlot(proj, reduction = "umap", features = "frac_overlap") + labs(title="Post-Harmony AMULET overlap fraction")
ggsave(output_paths[["umap_overlap_frac"]], plt, device = "pdf", width = 10, height = 7)

cm_nuc_genes <- c("RBM20", "TECRL", "MLIP", "CHRM2", "TRDN", "PALLD", "SGCD", "CMYA5", "MYOM2", "TBX5", "ESRRG",
"LINC02248", "KCNJ3", "TACC2", "CORIN", "DPY19L2", "WNK2", "MITF", "OBSCN", "FHOD3", "MYLK3",
"DAPK2", "NEXN")
# proj$cm_nuc_logcounts <- log10(colSums(GetAssayData(object = proj, slot = 'counts')[cm_nuc_genes, , drop = FALSE]) + 1)
proj$cm_nuc_pct <- PercentageFeatureSet(proj, features = cm_nuc_genes)
plt <- FeaturePlot(proj, reduction = "umap", features = "cm_nuc_pct") + labs(title="Post-Harmony Nuclear CM Marker read percent")
ggsave(output_paths[["umap_cm_nuc"]], plt, device = "pdf", width = 10, height = 7)

cm_cyto_genes <- c("TTN", "RYR2", "PAM", "TNNT2",
"RABGAP1L", "PDLIM5", "MYL7", "MYH6")
# proj$cm_cyto_logcounts <- log10(colSums(GetAssayData(object = proj, slot = 'counts')[cm_cyto_genes, , drop = FALSE]) + 1)
proj$cm_cyto_pct <- PercentageFeatureSet(proj, features = cm_cyto_genes)
plt <- FeaturePlot(proj, reduction = "umap", features = "cm_cyto_pct") + labs(title="Post-Harmony Cytoplasmic CM Marker read percent")
ggsave(output_paths[["umap_cm_cyto"]], plt, device = "pdf", width = 10, height = 7)

proj$cm_ratio <- log10(proj$cm_cyto_pct + 0.01) - log10(proj$cm_nuc_pct + 0.01)
plt <- FeaturePlot(proj, reduction = "umap", features = "cm_ratio") + labs(title="Post-Harmony Cytoplasmic-Nuclear CM Marker log ratio")
ggsave(output_paths[["umap_cm_ratio"]], plt, device = "pdf", width = 10, height = 7)

plt <- FeatureScatter(object = proj, feature1 = 'tss_enr', feature2 = 'cm_cyto_pct')
ggsave(output_paths[["scatter_tss_cm_cyto"]], plt, device = "pdf")

plt <- FeatureScatter(object = proj, feature1 = 'tss_enr', feature2 = 'cm_nuc_pct')
ggsave(output_paths[["scatter_tss_cm_nuc"]], plt, device = "pdf")


sink(type = "message")
sink()