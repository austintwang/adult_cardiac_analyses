Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(harmony)

params = snakemake@params 
input_paths = snakemake@input
output_paths = snakemake@output
log_paths = snakemake@log

set.seed(params[["seed"]])

read_fn <- function(path) {
    t <- read.table(file = path, sep = '\t', header = TRUE)
    t[, "barcode_rna", drop = FALSE]
}

tables <- lapply(input_paths[["data"]], read_fn)
head(tables[[1]])
data_all <- do.call(rbind, tables)

a <- data_all[["barcode_rna"]]
bc_counts <- as.data.frame((table(a)))
rownames(bc_counts) <- bc_counts[["a"]]

inds <- match(data_all[["barcode_rna"]], bc_counts[["a"]])
counts <- bc_counts[inds, "Freq", drop = FALSE]
rownames(counts) <- rownames(data_all)

write.table(proj@meta.data, file=output_paths[["data"]], quote=FALSE, sep='\t', col.names = NA)

sink(type = "message")
sink()