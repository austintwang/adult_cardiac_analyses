Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

params = snakemake@params 
input_paths = snakemake@input
output_paths = snakemake@output
log_paths = snakemake@log

set.seed(params[["seed"]])

genes <- read.table(file = params[["genes"]], sep = '\t', header = TRUE, comment.char = "")
head(genes) ####
proj <- readRDS(file = input_paths[["project_in"]])

dir.create(output_paths[["umaps"]])
for (index in seq_len(nrow(genes))) { 
    gene <- genes[index, "Name"]
    desc <- genes[index, "Description"]
    tryCatch(
        expr = {
            plt <- FeaturePlot(proj, features = gene, reduction = "umap") + labs(title=paste0(gene, ": ", desc))
            ggsave(file.path(output_paths[["umaps"]], paste0(gene, ".pdf")), plt, device = "pdf", width = 8, height = 7)
        },
        error = function(e){
            print(gene)
            print(e)
        }
    )
} 

plt <- DotPlot(proj, features = genes[["Name"]]) + RotatedAxis()
ggsave(output_paths[["dotplot"]], plt, device = "pdf", width = 11, height = 7)