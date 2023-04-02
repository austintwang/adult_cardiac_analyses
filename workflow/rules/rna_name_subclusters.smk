rule seurat_name_subclusters:
    """
    Seurat first-round RNA cluster naming
    """
    input:
        project_in = "results_subcluster_unified/{label}/rna/seurat_subcluster_2/proj.rds"
    output:
        project_out = "results_subcluster_unified/{label}/rna/seurat_name_subclusters/proj.rds",
        umap = "results_subcluster_unified/{label}/rna/seurat_name_subclusters/umap.pdf"
    params:
        seed = config["seurat_seed"],
        label_data = workflow.source_path("../files/subcluster_names.tsv")
    log:
        console = "logs/subcluster_unified/{label}/rna/seurat_name_subclusters/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_name_subclusters.R"

rule seurat_named_subcluster_markers:
    """
    Call subcluster markers
    """
    input:
        project_in = "results_subcluster_unified/{label}/rna/seurat_name_subclusters/proj.rds"
    output:
        umaps = directory("results_subcluster_unified/{label}/rna/seurat_named_subcluster_markers/umaps"),
        dotplot = "results_subcluster_unified/{label}/rna/seurat_named_subcluster_markers/dotplot.pdf",
        heatmap = "results_subcluster_unified/{label}/rna/seurat_named_subcluster_markers/heatmap.pdf",
        markers = "results_subcluster_unified/{label}/rna/seurat_named_subcluster_markers/markers.tsv"
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/subcluster_unified/{label}/rna/seurat_named_subcluster_markers/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_named_subcluster_markers.R"

rule seurat_subclusters_to_groups:
    """
    Add subclusters to integrated analyses
    """
    input:
        project_integrated = "results_groups/{group}/rna/seurat_name_group_1/proj.rds",
        projects_subcluster = expand("results_subcluster_unified/{label}/rna/seurat_name_subclusters/proj.rds", label=config["l1_labels"])
    output:
        project_out = "results_groups/{group}/rna/seurat_subclusters_to_groups/proj.rds",
        umap_l1 = "results_groups/{group}/rna/seurat_subclusters_to_groups/umap_l1.pdf",
        umap_l2 = "results_groups/{group}/rna/seurat_subclusters_to_groups/umap_l2.pdf"
    params:
        seed = config["seurat_seed"],
    log:
        console =  "logs/merged/{group}/rna/seurat_subclusters_to_groups/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_subclusters_to_groups.R"

rule seurat_named_cluster_markers:
    """
    Call subcluster markers
    """
    input:
        project_in = "results_groups/{group}/rna/seurat_subclusters_to_groups/proj.rds"
    output:
        umaps = directory("results_groups/{group}/rna/seurat_named_cluster_markers/umaps"),
        dotplot = "results_groups/{group}/rna/seurat_named_cluster_markers/dotplot.pdf",
        heatmap = "results_groups/{group}/rna/seurat_named_cluster_markers/heatmap.pdf",
        markers = "results_groups/{group}/rna/seurat_named_cluster_markers/markers.tsv"
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/merged/{group}/rna/seurat_named_cluster_markers/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_named_cluster_markers.R"

rule seurat_integrate_supergroup:
    """
    Integrate cardiomyocytes in region-specific way
    """
    input:
        projects_in = lambda w: [f"results_groups/{group}/rna/seurat_subclusters_to_groups/proj.rds" for group in supergroups[w.supergroup]]
    output:
        project_out = "results_supergroups/{supergroup}/{label}/rna/seurat_integrate_supergroup/proj.rds",
        umap_dataset_pre_harmony = "results_supergroups/{supergroup}/{label}/rna/seurat_integrate_supergroup/umap_dataset_pre_harmony.pdf",
        umap_mixing_pre_harmony = "results_supergroups/{supergroup}/{label}/rna/seurat_integrate_supergroup/umap_mixing_pre_harmony.pdf",
        umap_dataset_harmony = "results_supergroups/{supergroup}/{label}/rna/seurat_integrate_supergroup/umap_dataset_harmony.pdf",
        umap_mixing_harmony = "results_supergroups/{supergroup}/{label}/rna/seurat_integrate_supergroup/umap_mixing_harmony.pdf",
    params:
        seed = config["seurat_seed"],
        groups = lambda w: supergroups[w.supergroup]
    log:
        console = "logs/subcluster_supergroups/{supergroup}/{label}/rna/seurat_integrate_supergroup/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_integrate_subgroups.R"

rule seurat_subcluster_supergroup:
    """
    Seurat RNA subclustering
    """
    input:
        project_in = "results_supergroups/{supergroup}/{label}/rna/seurat_integrate_supergroup/proj.rds"
    output:
        project_out = "results_supergroups/{supergroup}/{label}/rna/seurat_subcluster_supergroup/proj.rds",
        umap = "results_supergroups/{supergroup}/{label}/rna/seurat_subcluster_supergroup/umap_clusters.pdf",
        umap_test = "results_supergroups/{supergroup}/{label}/rna/seurat_subcluster_supergroup/umap_clusters_test.pdf",
        umap_kramann_coarse = "results_supergroups/{supergroup}/{label}/rna/seurat_subcluster_supergroup/umap_kramann_coarse.pdf",
        umap_kramann_fine = "results_supergroups/{supergroup}/{label}/rna/seurat_subcluster_supergroup/umap_kramann_fine.pdf",
        umap_ellinor_coarse = "results_supergroups/{supergroup}/{label}/rna/seurat_subcluster_supergroup/umap_ellinor_coarse.pdf",
        umap_ellinor_fine = "results_supergroups/{supergroup}/{label}/rna/seurat_subcluster_supergroup/umap_ellinor_fine.pdf",
        umap_teichmann_coarse = "results_supergroups/{supergroup}/{label}/rna/seurat_subcluster_supergroup/umap_teichmann_coarse.pdf",
        umap_teichmann_fine = "results_supergroups/{supergroup}/{label}/rna/seurat_subcluster_supergroup/umap_teichmann_fine.pdf",
        umap_azimuth_coarse = "results_supergroups/{supergroup}/{label}/rna/seurat_subcluster_supergroup/umap_azimuth_coarse.pdf",
        umap_azimuth_fine = "results_supergroups/{supergroup}/{label}/rna/seurat_subcluster_supergroup/umap_azimuth_fine.pdf",
    params:
        seed = config["seurat_seed"],
        resolution = config["rna_cluster_resolution_cleanup"]
    log:
        console = "logs/subcluster_supergroups/{supergroup}/{label}/rna/seurat_subcluster_supergroup/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_subcluster_supergroup.R"

rule seurat_merge_label_plots_subcluster_supergroup:
    """
    Merge reference projection plot pdf's
    """
    input:
        "results_supergroups/{supergroup}/{label}/rna/seurat_subcluster_supergroup/umap_kramann_coarse.pdf",
        "results_supergroups/{supergroup}/{label}/rna/seurat_subcluster_supergroup/umap_kramann_fine.pdf",
        "results_supergroups/{supergroup}/{label}/rna/seurat_subcluster_supergroup/umap_ellinor_coarse.pdf",
        "results_supergroups/{supergroup}/{label}/rna/seurat_subcluster_supergroup/umap_ellinor_fine.pdf",
        "results_supergroups/{supergroup}/{label}/rna/seurat_subcluster_supergroup/umap_teichmann_coarse.pdf",
        "results_supergroups/{supergroup}/{label}/rna/seurat_subcluster_supergroup/umap_teichmann_fine.pdf",
        "results_supergroups/{supergroup}/{label}/rna/seurat_subcluster_supergroup/umap_azimuth_coarse.pdf",
        "results_supergroups/{supergroup}/{label}/rna/seurat_subcluster_supergroup/umap_azimuth_fine.pdf",
    output:
        "results_supergroups/{supergroup}/{label}/rna/seurat_subcluster_supergroup/label_plots.pdf"
    conda:
        "../envs/fetch.yaml"
    shell:
        "pdfunite {input} {output}; "

rule seurat_plot_subcluster_markers_supergroup:
    """
    Plot subclustering markers
    """
    input:
        project_in = "results_supergroups/{supergroup}/{label}/rna/seurat_subcluster_supergroup/proj.rds"
    output:
        umaps = directory("results_supergroups/{supergroup}/{label}/rna/seurat_plot_subcluster_markers_supergroup/umaps"),
        dotplot = "results_supergroups/{supergroup}/{label}/rna/seurat_plot_subcluster_markers_supergroup/dotplot.pdf",
        heatmap = "results_supergroups/{supergroup}/{supergroup}_unified/{label}/rna/seurat_plot_subcluster_markers_supergroup/heatmap.pdf",
        markers = "results_supergroups/{supergroup}/{label}/rna/seurat_plot_subcluster_markers_supergroup/markers.tsv"
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/subcluster_supergroups/{supergroup}/{label}/rna/seurat_plot_subcluster_markers_supergroup/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_plot_subcluster_markers.R"

rule seurat_plot_subcluster_genes_supergroup:
    """
    Plot external markers
    """
    input:
        project_in ="results_supergroups/{supergroup}/{label}/rna/seurat_subcluster_supergroup/proj.rds"
    output:
        umaps = directory("results_supergroups/{supergroup}/{label}/rna/seurat_plot_subcluster_genes_supergroup/umaps"),
        dotplot = "results_supergroups/{supergroup}/{label}/rna/seurat_plot_subcluster_genes_supergroup/dotplot.pdf",
        heatmap = "results_supergroups/{supergroup}/{label}/rna/seurat_plot_subcluster_genes_supergroup/heatmap.pdf",
    params:
        seed = config["seurat_seed"],
        genes = workflow.source_path("../files/plot_genes_subcluster.tsv")
    log:
        console = "logs/subcluster_supergroups/{supergroup}/{label}/rna/seurat_plot_genes/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_plot_subcluster_genes.R"

# rule seurat_integrate_supergroups:
#     """
#     Integrate RNA samples using Harmony
#     """
#     input:
#         projects_in = expand("results_subcluster/{supergroup}/{label}/rna/seurat_name_subclusters/proj.rds", supergroup=supergroup_names, allow_missing=True)
#     output:
#         project_out = "results_subcluster_unified/{label}/rna/seurat_integrate_supergroups/proj.rds",
#         umap_dataset_pre_harmony = "results_subcluster_unified/{label}/rna/seurat_integrate_supergroups/umap_dataset_pre_harmony.pdf",
#         umap_mixing_pre_harmony = "results_subcluster_unified/{label}/rna/seurat_integrate_supergroups/umap_mixing_pre_harmony.pdf",
#         umap_dataset_harmony = "results_subcluster_unified/{label}/rna/seurat_integrate_supergroups/umap_dataset_harmony.pdf",
#         umap_mixing_harmony = "results_subcluster_unified/{label}/rna/seurat_integrate_supergroups/umap_mixing_harmony.pdf",
#     params:
#         seed = config["seurat_seed"],
#     log:
#         console = "logs/subcluster_unified/{label}/rna/seurat_integrate_supergroups/console.log"
#     conda:
#         "../envs/seurat.yaml"
#     script:
#         "../scripts/seurat_integrate_supergroups.R"

# rule seurat_subcluster_unified:
#     """
#     Seurat RNA subclustering
#     """
#     input:
#         project_in = "results_subcluster_unified/{label}/rna/seurat_integrate_supergroups/proj.rds"
#     output:
#         project_out = "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/proj.rds",
#         umap = "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/umap_clusters.pdf",
#         umap_test = "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/umap_clusters_test.pdf",
#         umap_supergroup = "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/umap_labels_supergroup.pdf",
#         umap_kramann_coarse = "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/umap_kramann_coarse.pdf",
#         umap_kramann_fine = "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/umap_kramann_fine.pdf",
#         umap_ellinor_coarse = "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/umap_ellinor_coarse.pdf",
#         umap_ellinor_fine = "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/umap_ellinor_fine.pdf",
#         umap_teichmann_coarse = "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/umap_teichmann_coarse.pdf",
#         umap_teichmann_fine = "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/umap_teichmann_fine.pdf",
#         umap_azimuth_coarse = "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/umap_azimuth_coarse.pdf",
#         umap_azimuth_fine = "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/umap_azimuth_fine.pdf",
#     params:
#         seed = config["seurat_seed"],
#     log:
#         console = "logs/subcluster_unified/{label}/rna/seurat_subcluster_unified/console.log"
#     conda:
#         "../envs/seurat.yaml"
#     script:
#         "../scripts/seurat_subcluster_unified.R"

# rule seurat_merge_label_plots_subcluster_unified:
#     """
#     Merge reference projection plot pdf's
#     """
#     input:
#         "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/umap_kramann_coarse.pdf",
#         "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/umap_kramann_fine.pdf",
#         "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/umap_ellinor_coarse.pdf",
#         "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/umap_ellinor_fine.pdf",
#         "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/umap_teichmann_coarse.pdf",
#         "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/umap_teichmann_fine.pdf",
#         "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/umap_azimuth_coarse.pdf",
#         "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/umap_azimuth_fine.pdf",
#     output:
#         "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/label_plots.pdf"
#     conda:
#         "../envs/fetch.yaml"
#     shell:
#         "pdfunite {input} {output}; "
