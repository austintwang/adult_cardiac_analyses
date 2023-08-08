Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library("DESeq2")
library(sva)
library(ggplot2)

params = snakemake@params 
wildcards = snakemake@wildcards
input_paths = snakemake@input
output_paths = snakemake@output
log_paths = snakemake@log

set.seed(params[["seed"]])

mat_path <- input_paths[["mat"]]
metadata_path <- input_paths[["metadata"]]
coefficients_dir <- output_paths[["coefficients"]]
dispersion_plot_path <- output_paths[["dispersion_plot"]]
sv_plot_path <- output_paths[["sv_plot"]]


mat <- as.matrix(read.table(mat_path, header = TRUE, sep = "\t", row.names = 1))
metadata_chr <- read.table(metadata_path, header = TRUE, sep = "\t", row.names = 1)
metadata <- as.data.frame(unclass(metadata_chr), stringsAsFactors = TRUE)

print(colSums(mat)) ####

data <- tryCatch({
    dds <- DESeqDataSetFromMatrix(
        countData = mat,
        colData = metadata,
        design = ~ status + sex
    )
    list(dds=dds, incl_sex=TRUE)
}, error = function(e) {
    dds <- DESeqDataSetFromMatrix(
        countData = mat,
        colData = metadata,
        design = ~ status
    )
    list(dds=dds, incl_sex=FALSE)
})
dds <- data[["dds"]]
incl_sex <- data[["incl_sex"]]

dds <- dds[, dds$region %in% c("lv")]
dds$region <- droplevels(dds$region)

dds$status <- relevel(dds$status, ref = "healthy")

dds <- estimateSizeFactors(dds, type = "poscounts")

norm.cts <- counts(dds, normalized=TRUE)
if (incl_sex) {
    mm <- model.matrix(~ status + sex, colData(dds))
} else {
    mm <- model.matrix(~ status, colData(dds))
}
mm0 <- model.matrix(~ 1, colData(dds))
norm.cts <- norm.cts[rowSums(norm.cts) > 0,]
fit <- svaseq(norm.cts, mod=mm, mod0=mm0, n.sv=2)

dds$SV1 <- fit$sv[,1]
dds$SV2 <- fit$sv[,2]

dds_df <- as.data.frame(colData(dds))
if (incl_sex) {
    plt <- ggplot(dds_df, aes(x=SV1, y=SV2, color=status, shape=sex)) + geom_point()
} else {
    plt <- ggplot(dds_df, aes(x=SV1, y=SV2, color=status)) + geom_point()
}
ggsave(sv_plot_path, plt, device = "pdf")

if (incl_region) {
    design(dds) <- ~ status + sex + SV1 + SV2
} else {
    design(dds) <- ~ status + SV1 + SV2
}
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)

pdf(dispersion_plot_path)
plotDispEsts(dds)
dev.off()

coef <- resultsNames(dds) # lists the coefficients

for (c in coef) {
    if (c == "Intercept") {
        next
    }
    res_noshrink <- results(dds, name=c)
    res_shrink <- lfcShrink(dds, coef=c, type="apeglm")
    # res_shrink <- lfcShrink(dds, coef=c, type="normal")

    out_dir <- file.path(coefficients_dir, c)
    dir.create(out_dir, recursive = TRUE)

    res_noshrink_ordered <- as.data.frame(res_noshrink[order(res_noshrink$pvalue),])
    write.table(res_noshrink_ordered, file=file.path(out_dir, "noshrink.tsv"), quote=FALSE, sep='\t', col.names = NA)

    res_shrink_ordered <- as.data.frame(res_shrink[order(res_shrink$pvalue),])
    write.table(res_shrink_ordered, file=file.path(out_dir, "shrink.tsv"), quote=FALSE, sep='\t', col.names = NA)

    pdf(file.path(out_dir, "lfc_noshrink.pdf"))
    plotMA(res_noshrink, ylim=c(-2,2))
    dev.off() 

    pdf(file.path(out_dir, "lfc_shrink.pdf"))
    plotMA(res_shrink, ylim=c(-2,2))
    dev.off() 

    pdf(file.path(out_dir, "lfc_comparison.pdf"))
    plot(res_noshrink$log2FoldChange, res_shrink$log2FoldChange); abline(0,1)
    dev.off() 

}

