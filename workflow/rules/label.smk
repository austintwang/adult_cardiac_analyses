
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
        umap_kramann_fine = "results_merged/all/rna_subcluster/{cluster}/seurat_subcluster/umap_kramann_fine.pdf",
        umap_ellinor_fine = "results_merged/all/rna_subcluster/{cluster}/seurat_subcluster/umap_ellinor_fine.pdf",
        umap_teichmann_fine = "results_merged/all/rna_subcluster/{cluster}/seurat_subcluster/umap_teichmann_fine.pdf",
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
        "results_merged/all/rna_subcluster/{cluster}/seurat_subcluster/umap_clusters.pdf",
        "results_merged/all/rna_subcluster/{cluster}/seurat_subcluster/umap_kramann_fine.pdf",
        "results_merged/all/rna_subcluster/{cluster}/seurat_subcluster/umap_ellinor_fine.pdf",
        "results_merged/all/rna_subcluster/{cluster}/seurat_subcluster/umap_teichmann_fine.pdf",
    output:
        "results_merged/all/rna_subcluster/{cluster}/seurat_subcluster/label_plots.pdf"
    conda:
        "../envs/fetch.yaml"
    shell:
        "pdfunite {input} {output}; "

rule seurat_name_2:
    """
    Seurat first-round RNA cluster naming
    """
    input:
        project_in = "results_merged/all/rna_subcluster/{cluster}/seurat_subcluster/proj.rds"
    output:
        project_out = "results_merged/all/rna_subcluster/{cluster}/seurat_name_2/proj.rds",
        umap = "results_merged/all/rna_subcluster/{cluster}/seurat_name_2/umap.pdf"
    params:
        seed = config["seurat_seed"],
        label_data = workflow.source_path("../files/cluster_names_2.tsv")
    log:
        console = "logs/merged/all/rna_subcluster/{cluster}/seurat_name_2/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_name_2.R"

rule seurat_load_subclusters:
    """
    Add subclusters to integrated analyses
    """
    input:
        project_integrated = "results_merged/all/rna/seurat_integrate_l2/proj.rds",
        projects_subcluster = expand("results_merged/all/rna_subcluster/{cluster}/seurat_name_2/proj.rds", cluster=subcluster_jobs),
    output:
        project_out = "results_merged/all/rna/seurat_load_subclusters/proj.rds",
        umap = "results_merged/all/rna/seurat_load_subclusters/umap.pdf",
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/merged/all/rna/seurat_load_subclusters/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_load_subclusters.R"

rule seurat_l2_labels_to_subgroups:
    """
    Visualize integrated clustering in subgroups
    """
    input:
        project_integrated = "results_merged/all/rna/seurat_load_subclusters/proj.rds",
        project_subgroup = "results_merged/{group}/rna/seurat_cluster_rna/proj.rds"
    output:
        project_out = "results_merged/{group}/rna/seurat_l2_labels_to_subgroups/proj.rds",
        umap = "results_merged/{group}/rna/seurat_l2_labels_to_subgroups/umap.pdf",
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/merged/{group}/rna/seurat_l2_labels_to_subgroups/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_l2_labels_to_subgroups.R"