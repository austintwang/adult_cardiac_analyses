rule archr_build_group:
    """
    Build subgroup-level ArchR project
    """
    input:
        frags = lambda w: [f"results/{sample}/fetch/fragments.tsv.gz" for sample in groups[w.group]],
        frag_inds = lambda w: [f"results/{sample}/fetch/fragments.tsv.gz.tbi" for sample in groups[w.group]],
        seurat_data = "results_groups/{group}/rna/seurat_integrate_rna/metadata.tsv",
        blacklist = "resources/blacklist.bed",
        bsgenome_library_dir = "resources/bsgenome_lib",
        gene_anno = "resources/gene_annotation.rda"
    output:
        qc_dir = directory("results_groups/{group}/atac/archr_qc"),
        project_dir = directory("results_groups/{group}/atac/archr_init"),
        arrow_dir = directory("results_groups/{group}/atac/archr_arrows"),
    params:
        sample_names = lambda w: groups[w.group],
        seed = config["archr_seed"],
        min_frags = config["archr_min_frags"],
        min_tss_enr = config["archr_min_tss_enr"],
        bsgenome = config["bsgenome_name"],
        gene_anno = config["gene_anno_name"]
    log:
        console = "logs/merged/{group}/atac/archr_build/console.log",
        arrow_create = "logs/merged/{group}/atac/archr_build/arrow_create.log",
        save = "logs/merged/{group}/atac/archr_build/save.log"
    threads:
        max_threads
    conda:
        "../envs/archr.yaml"
    script:
        "../scripts/archr_build_group.R"

rule archr_embed:
    """
    ArchR dimensionality reduction and integration
    """
    input:
        project_in = "results_groups/{group}/atac/archr_init"
    output:
        project_out = directory("results_groups/{group}/atac/archr_embed")
    params:
        seed = config["archr_seed"]
    log:
        console = "logs/merged/{group}/atac/archr_embed/console.log",
        move = "logs/merged/{group}/atac/archr_embed/move.log",
        lsi_atac = "logs/merged/{group}/atac/archr_embed/lsi_atac.log",
        umap_plot = "logs/merged/{group}/atac/archr_embed/umap_plot.log",
        save = "logs/merged/{group}/atac/archr_embed/save.log"
    threads:
        max_threads
    conda:
        "../envs/archr.yaml"
    script:
        "../scripts/archr_embed.R"

rule archr_rna_clustering:
    """
    ArchR visualize RNA clustering in ATAC
    """
    input:
        project_in = "results_groups/{group}/atac/archr_embed",
        seurat_data = "results_groups/{group}/rna/seurat_cluster_rna/metadata.tsv"
    output:
        project_out = directory("results_groups/{group}/atac/archr_rna_clustering")
    params:
        seed = config["archr_seed"]
    log:
        console = "logs/merged/{group}/atac/archr_rna_clustering/console.log",
        move = "logs/merged/{group}/atac/archr_rna_clustering/move.log",
        umap_plot = "logs/merged/{group}/atac/archr_rna_clustering/umap_plot.log",
        cluster_atac = "logs/merged/{group}/atac/archr_rna_clustering/cluster_atac.log",
        doublets = "logs/merged/{group}/atac/archr_rna_clustering/doublets.log",
        save = "logs/merged/{group}/atac/archr_rna_clustering/save.log"
    threads:
        max_threads
    conda:
        "../envs/archr.yaml"
    script:
        "../scripts/archr_rna_clustering.R"

rule archr_plot_genes:
    """
    Plot external markers
    """
    input:
        project_in = "results_groups/{group}/atac/archr_rna_clustering"
    output:
        umaps = directory("results_groups/{group}/atac/archr_plot_genes/umaps"),
    params:
        seed = config["seurat_seed"],
        genes = workflow.source_path("../files/plot_genes_atac.tsv")
    log:
        console = "logs/merged/{group}/atac/archr_plot_genes/console.log"
    conda:
        "../envs/archr.yaml"
    script:
        "../scripts/archr_plot_genes.R"

rule archr_rna_labeling:
    """
    ArchR visualize RNA labeling in ATAC
    """
    input:
        project_in = "results_groups/{group}/atac/archr_rna_clustering",
        seurat_data = "results_groups/{group}/rna/seurat_label_rna/metadata.tsv"
    output:
        project_out = directory("results_groups/{group}/atac/archr_rna_labeling")
    params:
        seed = config["archr_seed"]
    log:
        console = "logs/merged/{group}/atac/archr_rna_labeling/console.log",
        move = "logs/merged/{group}/atac/archr_rna_labeling/move.log",
        umap_plot = "logs/merged/{group}/atac/archr_rna_labeling/umap_plot.log",
        cluster_atac = "logs/merged/{group}/atac/archr_rna_labeling/cluster_atac.log",
        doublets = "logs/merged/{group}/atac/archr_rna_labeling/doublets.log",
        save = "logs/merged/{group}/atac/archr_rna_labeling/save.log"
    threads:
        max_threads
    conda:
        "../envs/archr.yaml"
    script:
        "../scripts/archr_rna_labeling.R"