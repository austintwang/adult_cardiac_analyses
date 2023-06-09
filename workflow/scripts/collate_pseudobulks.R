Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library("plyr")

params = snakemake@params 
wildcards = snakemake@wildcards
input_paths = snakemake@input
output_paths = snakemake@output
log_paths = snakemake@log

set.seed(params[["seed"]])


mat_paths <- input_paths[["mats"]]
metadata_paths <- input_paths[["metadata"]]

mat_lst <- vector("list", length(mat_paths))
# print(length(mat_paths)) ####
# print(length(mat_paths)[[1]]) ####
# print(seq_along(length(mat_paths))) ####
# print(seq_along(length(mat_paths)[[1]])) ####
gene_names <- vector()
for (i in seq_along(mat_paths)) {
    print(i) ####
    mat <- read.table(mat_paths[[i]], header = TRUE, sep = "\t", row.names = 1)
    mat_lst[[i]] <- mat
    print(head(mat)) ####
    gene_names <- union(gene_names, rownames(mat))
}

mat_lst_aligned <- vector("list", length(mat_paths))
for (i in seq_along(mat_paths)) {
    mat_aln <- matrix(0, nrow = length(gene_names), ncol = ncol(mat_lst[[i]]))
    rownames(mat_aln) <- gene_names
    colnames(mat_aln) <- colnames(mat_lst[[i]])
    mat_aln[rownames(mat_lst[[i]]),] <- mat_lst[[i]]
    mat_lst_aligned[[i]] <- mat_aln
}

# print(mat_lst) ####

# cbind_op <- function(x, y) {cbind.fill(x, y, fill = 0)}
mat_merged <- do.call(cbind, mat_lst_aligned)
write.table(mat_merged, file=output_paths[["mat"]], quote=FALSE, sep='\t', col.names = NA)

meta_lst <- vector("list", length(metadata_paths))
for (i in length(metadata_paths)) {
    metadata <- read.table(metadata_paths[[i]], header = TRUE, sep = "\t", row.names = 1)
    meta_lst[[i]] <- metadata
}
meta_merged <- do.call(rbind, meta_lst)
write.table(meta_merged, file=output_paths[["metadata"]], quote=FALSE, sep='\t', col.names = NA)
