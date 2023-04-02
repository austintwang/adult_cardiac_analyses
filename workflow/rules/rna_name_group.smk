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

rule seurat_plot_genes:
    """
    Plot external markers
    """
    input:
        project_in = "results_groups/{group}/rna/seurat_name_group_1/proj.rds"
    output:
        umaps = directory("results_groups/{group}/rna/seurat_plot_genes/umaps"),
        dotplot = "results_groups/{group}/rna/seurat_plot_genes/dotplot.pdf",
        umap_clusters = "results_groups/{group}/rna/seurat_plot_genes/umap_clusters.pdf"
    params:
        seed = config["seurat_seed"],
        genes = workflow.source_path("../files/plot_genes.tsv")
    log:
        console = "logs/merged/{group}/rna/seurat_plot_genes/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_plot_genes.R"

# rule seurat_integrate_subgroups:
#     """
#     Integrate RNA samples using Harmony
#     """
#     input:
#         projects_in = lambda w: [f"results_groups/{group}/rna/seurat_name_group_1/proj.rds" for group in supergroups[w.supergroup]]
#     output:
#         project_out = "results_subcluster/{supergroup}/{label}/rna/seurat_integrate_subgroups/proj.rds",
#         umap_dataset_pre_harmony = "results_subcluster/{supergroup}/{label}/rna/seurat_integrate_subgroups/umap_dataset_pre_harmony.pdf",
#         umap_mixing_pre_harmony = "results_subcluster/{supergroup}/{label}/rna/seurat_integrate_subgroups/umap_mixing_pre_harmony.pdf",
#         umap_dataset_harmony = "results_subcluster/{supergroup}/{label}/rna/seurat_integrate_subgroups/umap_dataset_harmony.pdf",
#         umap_mixing_harmony = "results_subcluster/{supergroup}/{label}/rna/seurat_integrate_subgroups/umap_mixing_harmony.pdf",
#     params:
#         seed = config["seurat_seed"],
#         groups = lambda w: supergroups[w.supergroup]
#     log:
#         console = "logs/subcluster/{supergroup}/{label}/rna/seurat_integrate_subgroups/console.log"
#     conda:
#         "../envs/seurat.yaml"
#     script:
#         "../scripts/seurat_integrate_subgroups.R"

# rule seurat_subcluster_supergroup:
#     """
#     Seurat RNA subclustering
#     """
#     input:
#         project_in = "results_subcluster/{supergroup}/{label}/rna/seurat_integrate_subgroups/proj.rds"
#     output:
#         project_out = "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_supergroup/proj.rds",
#         umap = "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_supergroup/umap_clusters.pdf",
#         umap_test = "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_supergroup/umap_clusters_test.pdf",
#         umap_kramann_coarse = "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_supergroup/umap_kramann_coarse.pdf",
#         umap_kramann_fine = "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_supergroup/umap_kramann_fine.pdf",
#         umap_ellinor_coarse = "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_supergroup/umap_ellinor_coarse.pdf",
#         umap_ellinor_fine = "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_supergroup/umap_ellinor_fine.pdf",
#         umap_teichmann_coarse = "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_supergroup/umap_teichmann_coarse.pdf",
#         umap_teichmann_fine = "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_supergroup/umap_teichmann_fine.pdf",
#         umap_azimuth_coarse = "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_supergroup/umap_azimuth_coarse.pdf",
#         umap_azimuth_fine = "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_supergroup/umap_azimuth_fine.pdf",
#     params:
#         seed = config["seurat_seed"],
#     log:
#         console = "logs/subcluster/{supergroup}/{label}/rna/seurat_subcluster_supergroup/console.log"
#     conda:
#         "../envs/seurat.yaml"
#     script:
#         "../scripts/seurat_subcluster_supergroup.R"

# rule seurat_merge_label_plots_subcluster_supergroup:
#     """
#     Merge reference projection plot pdf's
#     """
#     input:
#         "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_supergroup/umap_kramann_coarse.pdf",
#         "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_supergroup/umap_kramann_fine.pdf",
#         "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_supergroup/umap_ellinor_coarse.pdf",
#         "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_supergroup/umap_ellinor_fine.pdf",
#         "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_supergroup/umap_teichmann_coarse.pdf",
#         "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_supergroup/umap_teichmann_fine.pdf",
#         "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_supergroup/umap_azimuth_coarse.pdf",
#         "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_supergroup/umap_azimuth_fine.pdf",
#     output:
#         "results_subcluster/{supergroup}/{label}/rna/seurat_subcluster_supergroup/label_plots.pdf"
#     conda:
#         "../envs/fetch.yaml"
#     shell:
#         "pdfunite {input} {output}; "

rule seurat_integrate_subgroups_unified:
    """
    Integrate RNA samples using Harmony
    """
    input:
        projects_in = expand("results_groups/{group}/rna/seurat_name_group_1/proj.rds", group=group_names)
    output:
        project_out = "results_subcluster_unified/{label}/rna/seurat_integrate_subgroups_unified/proj.rds",
        umap_dataset_pre_harmony = "results_subcluster_unified/{label}/rna/seurat_integrate_subgroups_unified/umap_dataset_pre_harmony.pdf",
        umap_mixing_pre_harmony = "results_subcluster_unified/{label}/rna/seurat_integrate_subgroups_unified/umap_mixing_pre_harmony.pdf",
        umap_dataset_harmony = "results_subcluster_unified/{label}/rna/seurat_integrate_subgroups_unified/umap_dataset_harmony.pdf",
        umap_mixing_harmony = "results_subcluster_unified/{label}/rna/seurat_integrate_subgroups_unified/umap_mixing_harmony.pdf",
    params:
        seed = config["seurat_seed"],
        groups = lambda w: group_names
    log:
        console = "logs/subcluster_unified/{label}/rna/seurat_integrate_subgroups_unified/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_integrate_subgroups_unified.R"

rule seurat_integrate_subgroups_qc_extras:
    """
    Integrate RNA samples using Harmony
    """
    input:
        project_in = "results_subcluster_unified/{label}/rna/seurat_integrate_subgroups_unified/proj.rds"
    output:
        umap_group_harmony = "results_subcluster_unified/{label}/rna/seurat_integrate_subgroups_unified/umap_group_harmony.pdf",
        umap_region_harmony = "results_subcluster_unified/{label}/rna/seurat_integrate_subgroups_unified/umap_region_harmony.pdf",
        umap_status_harmony = "results_subcluster_unified/{label}/rna/seurat_integrate_subgroups_unified/umap_status_harmony.pdf",
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/subcluster_unified/{label}/rna/seurat_integrate_subgroups_qc_extras/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_integrate_subgroups_qc_extras.R"

rule seurat_subcluster_unified:
    """
    Seurat RNA subclustering
    """
    input:
        project_in = "results_subcluster_unified/{label}/rna/seurat_integrate_subgroups_unified/proj.rds"
    output:
        project_out = "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/proj.rds",
        umap = "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/umap_clusters.pdf",
        umap_test = "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/umap_clusters_test.pdf",
        umap_kramann_coarse = "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/umap_kramann_coarse.pdf",
        umap_kramann_fine = "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/umap_kramann_fine.pdf",
        umap_ellinor_coarse = "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/umap_ellinor_coarse.pdf",
        umap_ellinor_fine = "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/umap_ellinor_fine.pdf",
        umap_teichmann_coarse = "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/umap_teichmann_coarse.pdf",
        umap_teichmann_fine = "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/umap_teichmann_fine.pdf",
        umap_azimuth_coarse = "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/umap_azimuth_coarse.pdf",
        umap_azimuth_fine = "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/umap_azimuth_fine.pdf",
    params:
        seed = config["seurat_seed"]
    log:
        console = "logs/subcluster_unified/{label}/rna/seurat_subcluster_unified/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_subcluster_unified.R"

rule seurat_merge_label_plots_subcluster_unified:
    """
    Merge reference projection plot pdf's
    """
    input:
        "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/umap_kramann_coarse.pdf",
        "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/umap_kramann_fine.pdf",
        "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/umap_ellinor_coarse.pdf",
        "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/umap_ellinor_fine.pdf",
        "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/umap_teichmann_coarse.pdf",
        "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/umap_teichmann_fine.pdf",
        "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/umap_azimuth_coarse.pdf",
        "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/umap_azimuth_fine.pdf",
    output:
        "results_subcluster_unified/{label}/rna/seurat_subcluster_unified/label_plots.pdf"
    conda:
        "../envs/fetch.yaml"
    shell:
        "pdfunite {input} {output}; "

