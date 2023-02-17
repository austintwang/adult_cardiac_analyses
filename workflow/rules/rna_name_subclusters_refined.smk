rule seurat_name_subclusters_supergroup:
    """
    Seurat RNA cluster naming
    """
    input:
        project_in = "results_supergroups/{supergroup}/{label}/rna/seurat_subcluster_supergroup/proj.rds"
    output:
        project_out = "results_supergroups/{supergroup}/{label}/rna/seurat_name_subclusters_supergroup/proj.rds",
        umap = "results_supergroups/{supergroup}/{label}/rna/seurat_name_subclusters_supergroup/umap.pdf"
    params:
        seed = config["seurat_seed"],
        label_data = workflow.source_path("../files/supergroup_subcluster_names.tsv")
    log:
        console = "logs/subcluster_supergroups/{supergroup}/rna/seurat_name_subclusters_supergroup/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_name_subclusters_supergroup.R"

rule seurat_named_subcluster_markers_supergroup:
    """
    Call subcluster markers
    """
    input:
        project_in = "results_supergroups/{supergroup}/{label}/rna/seurat_name_subclusters_supergroup/proj.rds"
    output:
        umaps = directory("results_supergroups/{supergroup}/{label}/rna/seurat_named_subcluster_markers_supergroup/umaps"),
        dotplot = "results_supergroups/{supergroup}/{label}/rna/seurat_named_subcluster_markers_supergroup/dotplot.pdf",
        heatmap = "results_supergroups/{supergroup}/{label}/rna/seurat_named_subcluster_markers_supergroup/heatmap.pdf",
        markers = "results_supergroups/{supergroup}/{label}/rna/seurat_named_subcluster_markers_supergroup/markers.tsv"
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/subcluster_supergroups/{supergroup}/{label}/rna/seurat_named_subcluster_markers_supergroup/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_named_subcluster_markers.R"

rule seurat_subclusters_supergroups_to_groups:
    """
    Add subclusters to integrated analyses
    """
    input:
        project_integrated = "results_groups/{group}/rna/seurat_subclusters_to_groups/proj.rds",
        projects_subcluster = expand("results_subcluster_unified/{label}/rna/seurat_name_subclusters/proj.rds", label=config["l1_labels"])
    output:
        project_out = "results_groups/{group}/rna/seurat_subclusters_supergroups_to_groups/proj.rds",
        umap_l1 = "results_groups/{group}/rna/seurat_subclusters_supergroups_to_groups/umap_l1.pdf",
        umap_l2 = "results_groups/{group}/rna/seurat_subclusters_supergroups_to_groups/umap_l2.pdf"
    params:
        seed = config["seurat_seed"],
    log:
        console =  "logs/merged/{group}/rna/seurat_subclusters_supergroups_to_groups/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_subclusters_supregroup_to_groups.R"

rule seurat_named_cluster_markers_refined:
    """
    Call subcluster markers
    """
    input:
        project_in = "results_groups/{group}/rna/seurat_subclusters_supergroups_to_groups/proj.rds"
    output:
        umaps = directory("results_groups/{group}/rna/seurat_named_cluster_markers_refined/umaps"),
        dotplot = "results_groups/{group}/rna/seurat_named_cluster_markers_refined/dotplot.pdf",
        heatmap = "results_groups/{group}/rna/seurat_named_cluster_markers_refined/heatmap.pdf",
        markers = "results_groups/{group}/rna/seurat_named_cluster_markers_refined/markers.tsv"
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/merged/{group}/rna/seurat_named_cluster_markers_refined/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_named_cluster_markers.R"