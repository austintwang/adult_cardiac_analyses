Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(Seurat)

params = snakemake@params 
wildcards = snakemake@wildcards
input_paths = snakemake@input
output_paths = snakemake@output
log_paths = snakemake@log

set.seed(params[["seed"]])

cell_types <- params[["cell_types"]]

out_mat_paths <- output_paths[["out_mats"]]
names(out_mat_paths) <- cell_types

out_meta_path <- output_paths[["out_metadata"]]

group <- strsplit(wildcards[["group"]], "-")[[1]]
print(group) ####
region <- group[[1]]
status <- group[[2]]
sex <- group[[3]]

proj <- readRDS(file = input_paths[["project_in"]])
DefaultAssay(object = proj) <- "RNA_test"

datasets <- unique(proj@meta.data[["dataset"]])
n_datasets <- length(datasets)

pseudobulk_metadata <- data.frame(
    dataset = datasets,
    region = rep(region, n_datasets),
    status = rep(status, n_datasets),
    sex = rep(sex, n_datasets),
    stringsAsFactors = FALSE
)
write.table(pseudobulk_metadata, file=out_meta_path, quote=FALSE, sep='\t', col.names = NA)

num_genes <- nrow(proj)
out_mats <- list()
for (c in cell_types) {
    out_mat <- matrix(, nrow = num_genes, ncol = n_datasets)
    rownames(out_mat) <- rownames(proj)
    colnames(out_mat) <- datasets
    for (d in datasets) {
        sel <- (proj@meta.data[["dataset"]] == d) & (proj@meta.data[["cell_types_l1"]] == c)
        out_mat[,d] <- rowSums(GetAssayData(object = proj, slot = "counts")[,sel,drop=FALSE])
    }
    out_mats[[c]] <- out_mat
}

for (c in cell_types) {
    write.table(out_mats[[c]], file=out_mat_paths[[c]], quote=FALSE, sep='\t', col.names = NA)
}