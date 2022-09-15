Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

params = snakemake@params 
input_paths = snakemake@input
output_paths = snakemake@output
log_paths = snakemake@log

set.seed(params[["seed"]])

read_fn <- function(path, sample) {
    data <- read.table(file = path, sep = '\t', quote = "", header = TRUE, row.names=1)
    rownames(data) <- paste0(sample, "_", rownames(data))
    # print(head(data)) ####
    data
}

tables <- mapply(read_fn, input_paths[["tables"]], params[["samples"]], SIMPLIFY = FALSE)
print(tables) ####
# head(tables[[1]]) ####
data_stacked <- do.call(rbind, tables)
head(data_stacked) ####

write.table(data_stacked, file = output_paths[["data_out"]], quote=FALSE, sep='\t', col.names = NA)

sink(type = "message")
sink()