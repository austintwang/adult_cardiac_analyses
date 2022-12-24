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

stallion <- c("1"="#D51F26","2"="#272E6A","3"="#208A42","4"="#89288F","5"="#F47D2B", "6"="#FEE500","7"="#8A9FD1","8"="#C06CAB","19"="#E6C2DC",
               "10"="#90D5E4", "11"="#89C75F","12"="#F37B7D","13"="#9983BD","14"="#D24B27","15"="#3BBCA8", "16"="#6E4B9E","17"="#0C727C", "18"="#7E1416","9"="#D8A767","20"="#3D3D3D")
stallion <- unname(stallion)

plot_fn <- function(object, group, reduction, colors) {
    cats <- sort(unique(object@meta.data[[group]]))
    colors_out <- rep_len(colors, length(cats))
    names(colors_out) <- cats
    DimPlot(object, reduction = reduction, group.by = group, label = TRUE, cols = colors_out, pt.size=0.1)
}

genes <- read.table(file = params[["genes"]], sep = '\t', header = TRUE, comment.char = "")
head(genes) ####
proj <- readRDS(file = input_paths[["project_in"]])
DefaultAssay(object = proj) <- "RNA_test"

plt <- plot_fn(proj, "cell_types_1", "umap_test", stallion) + labs(title="Cell Types (Level-1), Test Set")
ggsave(output_paths[["umap_clusters"]], plt, device = "pdf", width = 10, height = 7)

dir.create(output_paths[["umaps"]])
for (index in seq_len(nrow(genes))) { 
    gene <- genes[index, "Name"]
    desc <- genes[index, "Description"]
    tryCatch(
        expr = {
            plt <- FeaturePlot(proj, features = gene, reduction = "umap_test") + labs(title=paste0(gene, ": ", desc))
            ggsave(file.path(output_paths[["umaps"]], paste0(gene, ".pdf")), plt, device = "pdf", width = 8, height = 7)
        },
        error = function(e){
            print(gene)
            print(e)
        }
    )
} 

plt <- DotPlot(proj, features = genes[["Name"]], group.by = "cell_types_1") + RotatedAxis()
ggsave(output_paths[["dotplot"]], plt, device = "pdf", width = 11, height = 7)