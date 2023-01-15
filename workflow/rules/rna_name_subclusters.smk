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

rule seurat_named_subcluster_markers:
    """
    Call subcluster markers
    """
    input:
        project_in = "results_subcluster/{supergroup}/{label}/rna/seurat_name_subclusters/proj.rds"
    output:
        umaps = directory("results_subcluster/{supergroup}/{label}/rna/seurat_named_subcluster_markers/umaps"),
        dotplot = "results_subcluster/{supergroup}/{label}/rna/seurat_named_subcluster_markers/dotplot.pdf",
        heatmap = "results_subcluster/{supergroup}/{label}/rna/seurat_named_subcluster_markers/heatmap.pdf",
        markers = "results_subcluster/{supergroup}/{label}/rna/seurat_named_subcluster_markers/markers.tsv"
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/subcluster/{supergroup}/{label}/rna/seurat_named_subcluster_markers/console.log"
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
        projects_subcluster = lambda w: [f"results_subcluster/{group_to_supergroup[w.group]}/{label}/rna/seurat_name_subclusters/proj.rds" for label in config["l1_labels"]]
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
        markers = "rresults_groups/{group}/rna/seurat_named_cluster_markers/markers.tsv"
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/merged/{group}/rna/seurat_named_cluster_markers/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_named_cluster_markers.R"