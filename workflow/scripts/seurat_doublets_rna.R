Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(remotes)

remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(DoubletFinder)

params = snakemake@params 
input_paths = snakemake@input
output_paths = snakemake@output
log_paths = snakemake@log

set.seed(params[["seed"]])

proj <- readRDS(file = input_paths[["project_in"]])

print(proj) ####

# pK identification (no ground-truth)
sweep.list <- paramSweep_v3(proj, PCs = 1:30, sct = TRUE, num.cores = snakemake@threads)
sweep.stats <- summarizeSweep(sweep.list)
bcmvn <- find.pK(sweep.stats)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic doublet proportion estimate
annotations <- proj@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 
nExp.poi <- round(params[["doublet_rate"]] * nrow(proj@meta.data)) 
nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))

# run DoubletFinder
proj <- doubletFinder_v3(
    seu = proj, 
    PCs = 1:30, 
    pK = optimal.pk,
    nExp = nExp.poi.adj
)
print(head(proj@meta.data)) ####
doublets_col <- grep("DF.classifications", names(proj@meta.data), value = TRUE)
proj$doublet_doubletfinder <- proj[[doublets_col]] != "Singlet"
proj$doublet_amulet <- proj$amulet_qval > params[["amulet_fdr"]]
proj$doublet_union <- proj$doublet_doubletfinder | proj$doublet_amulet
proj$doublet_intersect <- proj$doublet_doubletfinder & proj$doublet_amulet

print(head(proj@meta.data)) ####
proj$doublet_status <- rep(NA, length(Cells(proj)))
proj$doublet_status[proj$doublet_amulet] <- "Amulet only"
proj$doublet_status[proj$doublet_doubletfinder] <- "DoubletFinder only"
proj$doublet_status[proj$doublet_intersect] <- "Amulet and DoubletFinder"
print(head(proj@meta.data)) ####

plt <- DimPlot(proj, reduction = "umap", group.by = "doublet_status")
ggsave(output_paths[["umap"]], plt, device = "pdf")

proj_filtered <- subset(proj, subset = doublet_union == FALSE)

proj_filtered <- NormalizeData(proj_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
proj_filtered <- FindVariableFeatures(proj_filtered, selection.method = "vst", nfeatures = 2000)
proj_filtered <- ScaleData(proj_filtered)

proj_filtered <- RunPCA(proj_filtered, features = VariableFeatures(object = proj))

proj_filtered <- FindNeighbors(proj_filtered, dims = 1:30)
proj_filtered <- FindClusters(object = proj_filtered) 

proj_filtered <- RunUMAP(proj_filtered, dims = 1:30, return.model = TRUE)

plt <- DimPlot(proj_filtered, reduction = "umap", group.by = "seurat_clusters")
ggsave(output_paths[["umap_filtered"]], plt, device = "pdf")

print(proj)
print(proj_filtered)
saveRDS(proj, file = output_paths[["project_out_all"]])
saveRDS(proj_filtered, file = output_paths[["project_out_filtered"]])

write.table(proj_filtered@meta.data, file=output_paths[["metadata"]], quote=FALSE, sep='\t', col.names = NA)

sink(type = "message")
sink()