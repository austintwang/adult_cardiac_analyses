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

genes <- read.table(file = params[["genes"]], sep = '\t', header = TRUE, comment.char = "")
head(genes) ####

proj <- loadArchRProject(path = input_paths[["project_in"]], force = FALSE, showLogo = FALSE)

dir.create(output_paths[["umaps"]])
for (index in seq_len(nrow(genes))) { 
    gene <- genes[index, "Name"]
    desc <- genes[index, "Description"]
    tryCatch(
        expr = {
            p <- plotEmbedding(
                ArchRProj = proj, 
                colorBy = "GeneScoreMatrix", 
                name = gene, 
                embedding = "UMAP_Harmony",
                quantCut = c(0.01, 0.95),
                imputeWeights = NULL
            )
            plt <- p + labs(title=paste0(gene, ": ", desc))
            ggsave(file.path(output_paths[["umaps"]], paste0(gene, ".pdf")), plt, device = "pdf", width = 5, height = 5)
        },
        error = function(e){
            print(gene)
            print(e)
        }
    )
} 

sink(type = "message")
sink()