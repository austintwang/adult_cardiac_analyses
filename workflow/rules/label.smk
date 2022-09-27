
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

rule seurat_subcluster:
    """
    Seurat RNA clustering
    """
    input:
        project_in = "results_merged/all/rna/seurat_name_1/proj.rds"
    output:
        project_out = "results_merged/all/rna_subcluster/{cluster}/seurat_subcluster/proj.rds",
        umap = "results_merged/all/rna_subcluster/{cluster}/seurat_subcluster/umap_clusters.pdf",
        umap_kramann_coarse = "results_merged/all/rna_subcluster/{cluster}/seurat_subcluster/umap_kramann_fine.pdf",
        umap_ellinor_coarse = "results_merged/all/rna_subcluster/{cluster}/seurat_subcluster/umap_ellinor_fine.pdf",
        umap_teichmann_coarse = "results_merged/all/rna_subcluster/{cluster}/seurat_subcluster/umap_teichmann_fine.pdf",
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/merged/all/rna_subcluster/{cluster}/seurat_subcluster/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_subcluster.R"

rule seurat_merge_label_plots_subcluster:
    """
    Merge reference projection plot pdf's
    """
    input:
        umap = "results_merged/all/rna_subcluster/{cluster}/seurat_subcluster/umap_clusters.pdf",
        umap_kramann_coarse = "results_merged/all/rna_subcluster/{cluster}/seurat_subcluster/umap_kramann_fine.pdf",
        umap_ellinor_coarse = "results_merged/all/rna_subcluster/{cluster}/seurat_subcluster/umap_ellinor_fine.pdf",
        umap_teichmann_coarse = "results_merged/all/rna_subcluster/{cluster}/seurat_subcluster/umap_teichmann_fine.pdf",
    output:
        "results_merged/all/rna_subcluster/{cluster}/seurat_subcluster/label_plots.pdf"
    conda:
        "../envs/fetch.yaml"
    shell:
        "pdfunite {input} {output}; "
