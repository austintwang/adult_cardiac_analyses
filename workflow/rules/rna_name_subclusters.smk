rule seurat_name_subclusters:
    """
    Seurat first-round RNA cluster naming
    """
    input:
        project_in = "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_2/proj.rds"
    output:
        project_out = "results_subcluster/{supergroup}/{label}/rna/seurat_name_subclusters/proj.rds",
        umap = "results_subcluster/{supergroup}/{label}/rna/seurat_name_subclusters/umap.pdf"
    params:
        seed = config["seurat_seed"],
        label_data = workflow.source_path("../files/subcluster_names.tsv")
    log:
        console = "logs/subcluster/{supergroup}/{label}/rna/seurat_name_subclusters/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_name_subclusters.R"