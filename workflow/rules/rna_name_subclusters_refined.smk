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
        console = "logs/subcluster_supergroups/{supergroup}/{label}/rna/seurat_name_subclusters_supergroup/console.log"
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
        projects_subcluster = expand("results_supergroups/{supergroup}/{label}/rna/seurat_name_subclusters_supergroup/proj.rds", supergroup=supergroup_names, label=config["supergroup_labels"])
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
        "../scripts/seurat_subclusters_supergroup_to_groups.R"

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

rule seurat_write_l1_markers:
    """
    Write level-1 cell type marker genes
    """
    input:
        project_in = "results_groups/{group}/rna/seurat_write_l1_markers/proj.rds"
    output:
        markers = directory("results_groups/{group}/rna/seurat_write_l1_markers/markers")
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/merged/{group}/rna/seurat_write_l1_markers/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_write_l1_markers.R"

rule seurat_embed_all:
    """
    Build RNA embeddings with all data 
    """
    input:
        project_in = "results_groups/{group}/rna/seurat_subclusters_supergroups_to_groups/proj.rds"
    output:
        project_out = "results_groups/{group}/rna/seurat_embed_all/proj.rds",
        umap_dataset = "results_groups/{group}/rna/seurat_embed_all/umap_dataset.pdf",
        umap_cell_types = "results_groups/{group}/rna/seurat_embed_all/umap_cell_types.pdf",
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/merged/{group}/rna/umap_cell_types/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_embed_all.R"

rule seurat_write_embeddings:
    """
    Seurat save RNA embeddings
    """
    input:
        project_in = "results_groups/{group}/rna/seurat_embed_all/proj.rds"
    output:
        rna_pca = "results_merged/{group}/rna/seurat_write_embeddings/rna_pca.tsv",
        rna_harmony = "results_merged/{group}/rna/seurat_write_embeddings/rna_harmony.tsv",
        rna_umap = "results_merged/{group}/rna/seurat_write_embeddings/rna_umap.tsv",
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/merged/{group}/rna/seurat_write_embeddings/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_write_embeddings.R"