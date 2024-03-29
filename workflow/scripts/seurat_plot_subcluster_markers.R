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

# genes <- read.table(file = params[["genes"]], sep = '\t', header = TRUE, comment.char = "")
# head(genes) ####
proj <- readRDS(file = input_paths[["project_in"]])
DefaultAssay(object = proj) <- "RNA_test"

# markers <- FindAllMarkers(proj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, return.thresh = 0.05)
markers <- FindAllMarkers(proj, min.pct = 0.25, logfc.threshold = 0.25, return.thresh = 0.05)

write.table(markers, file=output_paths[["markers"]], quote=FALSE, sep='\t', col.names = NA)
# markers %>%
#     group_by(cluster) %>%
#     slice_max(n = 2, order_by = avg_log2FC)

if (length(markers) == 0) {
    dir.create(output_paths[["umaps"]])
    pdf(output_paths[["dotplot"]])
    # plot(x=c(1,2,4,2,5))
    dev.off()
    pdf(output_paths[["heatmap"]])
    # plot(x=c(1,2,4,2,5))
    dev.off()
} else {
    markers %>%
        group_by(cluster) %>%
        top_n(n = 5, wt = avg_log2FC) -> top10

    genes <- unique(top10$gene)
    head(genes) ####

    dir.create(output_paths[["umaps"]])
    for (index in seq_len(length(genes))) { 
        gene <- genes[[index]]
        tryCatch(
            expr = {
                plt <- FeaturePlot(proj, features = gene, reduction = "umap_test") + labs(title=gene)
                ggsave(file.path(output_paths[["umaps"]], paste0(gene, ".pdf")), plt, device = "pdf", width = 8, height = 7)
            },
            error = function(e){
                print(gene)
                print(e)
            }
        )
    } 

    plt <- DotPlot(proj, features = genes) + RotatedAxis()
    ggsave(output_paths[["dotplot"]], plt, device = "pdf", width = 30, height = 7)

    proj <- ScaleData(object = proj, features = rownames(proj))
    plt <- DoHeatmap(subset(proj, downsample = 100), features = genes, size = 3)
    ggsave(output_paths[["heatmap"]], plt, device = "pdf", width = 10, height = 10, limitsize = FALSE)
}
