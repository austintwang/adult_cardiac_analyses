rule seurat_subcluster_cleanup:
    """
    Seurat clean up subclusters
    """
    input:
        project_in = "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_supergroup/proj.rds"
    output:
        project_out = "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_cleanup/proj.rds"
    params:
        seed = config["seurat_seed"],
        label_data = workflow.source_path("../files/subcluster_cleanup.tsv")
    log:
        console = "logs/subcluster/{supergroup}/{label}/rna/seurat_subcluster_cleanup/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_subcluster_cleanup.R"

rule seurat_subcluster_2:
    """
    Seurat RNA subclustering
    """
    input:
        project_in = "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_cleanup/proj.rds"
    output:
        project_out = "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_2/proj.rds",
        umap = "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_2/umap_clusters.pdf",
        umap_test = "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_2/umap_clusters_test.pdf",
        umap_kramann_coarse = "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_2/umap_kramann_coarse.pdf",
        umap_kramann_fine = "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_2/umap_kramann_fine.pdf",
        umap_ellinor_coarse = "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_2/umap_ellinor_coarse.pdf",
        umap_ellinor_fine = "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_2/umap_ellinor_fine.pdf",
        umap_teichmann_coarse = "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_2/umap_teichmann_coarse.pdf",
        umap_teichmann_fine = "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_2/umap_teichmann_fine.pdf",
        umap_azimuth_coarse = "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_2/umap_azimuth_coarse.pdf",
        umap_azimuth_fine = "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_2/umap_azimuth_fine.pdf",
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/subcluster/{supergroup}/{label}/rna/seurat_subcluster_2/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_subcluster_supergroup.R"

rule seurat_merge_label_plots_subcluster_2:
    """
    Merge reference projection plot pdf's
    """
    input:
        "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_2/umap_kramann_coarse.pdf",
        "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_2/umap_kramann_fine.pdf",
        "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_2/umap_ellinor_coarse.pdf",
        "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_2/umap_ellinor_fine.pdf",
        "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_2/umap_teichmann_coarse.pdf",
        "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_2/umap_teichmann_fine.pdf",
        "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_2/umap_azimuth_coarse.pdf",
        "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_2/umap_azimuth_fine.pdf",
    output:
        "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_2/label_plots.pdf"
    conda:
        "../envs/fetch.yaml"
    shell:
        "pdfunite {input} {output}; "

rule seurat_plot_subcluster_markers:
    """
    Plot subclustering markers
    """
    input:
        project_in = "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_2/proj.rds"
    output:
        umaps = directory("results_subcluster/{supergroup}/{label}/rna/seurat_plot_subcluster_markers/umaps"),
        dotplot = "results_subcluster/{supergroup}/{label}/rna/seurat_plot_subcluster_markers/dotplot.pdf",
        markers = "results_subcluster/{supergroup}/{label}/rna/seurat_plot_subcluster_markers/markers.tsv"
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/subcluster/{supergroup}/{label}/rna/seurat_plot_subcluster_markers/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_plot_subcluster_markers.R"