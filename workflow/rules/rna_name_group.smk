rule seurat_name_1:
    """
    Seurat first-round RNA cluster naming
    """
    input:
        project_in =  "results_groups/{group}/rna/seurat_add_batch_data/proj.rds"
    output:
        project_out = "results_merged/{group}/rna/seurat_name_group/proj.rds",
        umap_full = "results_merged/all/{group}/seurat_name_group/umap_full.pdf",
        umap_merged = "results_merged/{group}/rna/seurat_name_group/umap_merged.pdf"
    params:
        seed = config["seurat_seed"],
        label_data = workflow.source_path("../files/cluster_names_group.txt")
    log:
        console = "logs/merged/{group}/rna/seurat_name_group/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_name_group.R"