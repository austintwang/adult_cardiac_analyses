
rule seurat_name_1:
    """
    Seurat first-round RNA cluster naming
    """
    input:
        project_in = "results_merged/all/rna/seurat_cluster_rna/proj.rds"
    output:
        project_out = "results_merged/all/rna/seurat_name_1/proj.rds",
        umap = "results_merged/all/rna/seurat_name_1/umap.pdf"
    params:
        seed = config["seurat_seed"],
        label_data = workflow.source_path("../files/cluster_names_1.tsv")
    log:
        console = "logs/merged/all/rna/seurat_name_1/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_name_1.R"