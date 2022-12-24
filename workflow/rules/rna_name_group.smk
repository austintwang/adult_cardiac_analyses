rule seurat_name_group_1:
    """
    Seurat first-round RNA cluster naming
    """
    input:
        project_in =  "results_groups/{group}/rna/seurat_add_batch_data/proj.rds"
    output:
        project_out = "results_groups/{group}/rna/seurat_name_group_1/proj.rds",
        umap = "results_groups/{group}/rna/seurat_name_group_1/umap.pdf",
    params:
        seed = config["seurat_seed"],
        label_data = workflow.source_path("../files/cluster_names_group_1.txt")
    log:
        console = "logs/merged/{group}/rna/seurat_name_group_1/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_name_group_1.R"

# rule seurat_plot_genes:
#     """
#     Plot external markers
#     """
#     input:
#         project_in = "results_groups/{group}/rna/seurat_name_group_1/proj.rds"
#     output:
#         umaps = directory("results_groups/{group}/rna/seurat_plot_genes/umaps"),
#         dotplot = "results_groups/{group}/rna/seurat_plot_genes/dotplot.pdf",
#         umap_clusters = "results_groups/{group}/rna/seurat_plot_genes/umap_clusters.pdf"
#     params:
#         seed = config["seurat_seed"],
#         genes = workflow.source_path("../files/plot_genes.tsv")
#     log:
#         console = "logs/merged/{group}/rna/seurat_plot_genes/console.log"
#     conda:
#         "../envs/seurat.yaml"
#     script:
#         "../scripts/seurat_plot_genes.R"

# rule seurat_visualize_seq_emb:
#     """
#     Visualize labels in ATAC sequence space
#     """
#     input:
#         project_in =  "results_groups/{group}/rna/seurat_add_batch_data/proj.rds"
#     output:
#         project_out = "results_groups/{group}/rna/seurat_name_group/proj.rds",
#         umap_full = "results_groups/{group}/rna/seurat_name_group/umap_full.pdf",
#         umap_merged = "results_groups/{group}/rna/seurat_name_group/umap_merged.pdf"
#     params:
#         seed = config["seurat_seed"],
#         label_data = workflow.source_path("../files/cluster_names_group.txt")
#     log:
#         console = "logs/merged/{group}/rna/seurat_name_group/console.log"
#     conda:
#         "../envs/seurat.yaml"
#     script:
#         "../scripts/seurat_name_group.R"