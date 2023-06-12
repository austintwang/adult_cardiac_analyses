Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library("DESeq2")

params = snakemake@params 
wildcards = snakemake@wildcards
input_paths = snakemake@input
output_paths = snakemake@output
log_paths = snakemake@log

set.seed(params[["seed"]])

mat_path <- input_paths[["mat"]]
metadata_path <- input_paths[["metadata"]]

mat <- as.matrix(read.table(mat_path, header = TRUE, sep = "\t", row.names = 1))
metadata_chr <- read.table(metadata_path, header = TRUE, sep = "\t", row.names = 1)
metadata <- as.data.frame(unclass(metadata_chr), stringsAsFactors = TRUE)

dds <- DESeqDataSetFromMatrix(
    countData = mat,
    colData = metadata,
    design = ~ region + status
    # design = ~ region + status + sex
)

cts <- counts(dds)
geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
print(geoMeans) ####
dds <- estimateSizeFactors(dds, geoMeans=geoMeans)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)

print(resultsNames(dds) )# lists the coefficients





# res <- results(dds, name="condition_trt_vs_untrt")
# # or to shrink log fold changes association with condition:
# res <- lfcShrink(dds, coef="condition_trt_vs_untrt", type="apeglm")