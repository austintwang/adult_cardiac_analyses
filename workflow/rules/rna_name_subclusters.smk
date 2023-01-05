rule seurat_name_subclusters:
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