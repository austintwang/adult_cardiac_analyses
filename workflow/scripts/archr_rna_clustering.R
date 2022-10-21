Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(parallel)

library(ArchR)
library(pheatmap)

addArchRVerbose(verbose = FALSE)
addArchRChrPrefix(chrPrefix = FALSE)

# Disable HDF5 file locking
# Workaround for HDF5 I/O issues on NFS
# https://github.com/GreenleafLab/ArchR/issues/248#issuecomment-789453997
Sys.setenv("HDF5_USE_FILE_LOCKING" = "FALSE")
Sys.setenv("RHDF5_USE_FILE_LOCKING" = "FALSE")

##########

params = snakemake@params 
input_paths = snakemake@input
output_paths = snakemake@output
threads = snakemake@threads
log_paths = snakemake@log

seed <- params[["seed"]]
set.seed(seed)

addArchRThreads(threads = threads)

##########

# Load and move project
proj <- loadArchRProject(path = input_paths[["project_in"]], force = FALSE, showLogo = FALSE)
proj <- saveArchRProject(
    ArchRProj = proj,
    outputDirectory = output_paths[["project_out"]],
    overwrite = TRUE,
    load = TRUE,
    logFile = log_paths[["move"]],
)

##########
seurat_data <- read.table(file = input_paths[["seurat_data"]], sep = '\t', header = TRUE)
bc_names <- paste0(seurat_data[["dataset"]], "#", seurat_data[["barcode_atac"]])
clusters <- as.character(seurat_data[["seurat_clusters"]])


proj <- addCellColData(
    ArchRProj = proj,
    data = clusters,
    name = "seurat_clusters",
    cells = bc_names,
    force = FALSE
)

# Plot ATAC clusters
p1 <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "cellColData", 
    name = "seurat_clusters", 
    embedding = "UMAP_Harmony",
    logFile = log_paths[["umap_plot"]]
)

plotPDF(p1, name = "umap_seurat_clusters.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

# Cluster cells by ATAC LSI values
proj <- addClusters(
    input = proj,
    reducedDims = "LSI_ATAC",
    method = "Seurat",
    name = "Clusters_ATAC",
    resolution = 0.8,
    logFile = log_paths[["cluster_atac"]]
)


cM <- confusionMatrix(proj$Clusters_ATAC, proj$seurat_clusters)
cM <- as.matrix(cM)
print(cM) ####
print(colnames(cM)) ####
cM <- cM[, !is.na(colnames(cM))]
cM <- cbind(cM, Unknown = 0.01)
print(cM) ####

p <- pheatmap::pheatmap(
    mat = cM / rowSums(cM), 
    color = paletteContinuous("whiteBlue"), 
    border_color = "black"
)
plotPDF(p, name = "seurat_label_cm.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

##########

saveArchRProject(
    ArchRProj = proj,
    outputDirectory = output_paths[["project_out"]],
    overwrite = TRUE,
    load = FALSE,
    logFile = log_paths[["save"]],
)

sink(type = "message")
sink()