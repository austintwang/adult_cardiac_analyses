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
for (i in seq_along(length(mat_paths))) {
    mat <- read.table(i, header = TRUE, sep = "\t", row.names = 1)
    mat_lst[[i]] <- mat
}

cbind_op <- function(x, y) {cbind.fill(x, y, fill = 0)}
mat_merged <- do.call(cbind_op, mat_lst)
write.table(mat_merged, file=output_paths[["mat"]], quote=FALSE, sep='\t', col.names = NA)

meta_lst <- vector("list", length(metadata_paths))
for (i in length(metadata_paths)) {
    metadata <- read.table(i, header = TRUE, sep = "\t", row.names = 1)
    meta_lst[[i]] <- metadata
}
meta_merged <- do.call(rbind, met_lst)
write.table(meta_merged, file=output_paths[["metadata"]], quote=FALSE, sep='\t', col.names = NA)
