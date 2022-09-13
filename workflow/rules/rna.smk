rule seurat_install_doubletfinder:
    """
    Install doubletfinder R package
    """
    output:
        doubletfinder_library_dir = directory("resources/doubletfinder_lib")
    log:
        console = "logs/resources/seurat_install_doubletfinder/console.log"
    conda:
        "../envs/archr.yaml"
    script:
        "../scripts/seurat_install_doubletfinder.R"

rule seurat_build_reference_kramann:
    """
    Build Seurat reference dataset
    """
    input:
        seurat_object = "reference/kramann/fetch/seurat.rds",
    output:
        project_out = "reference/kramann/seurat_build_reference/proj.rds",
        qc_violin = "reference/kramann/seurat_build_reference/qc_violin.pdf",
        qc_scatter = "reference/kramann/seurat_build_reference/qc_scatter.pdf",
        umap = "reference/kramann/seurat_build_reference/umap.pdf"
    params:
        seed = config["seurat_seed"]
    log:
        console = "logs/reference/kramann/seurat_build_reference/seurat_build_rna/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_build_reference_kramann.R"

rule seurat_build_reference_ellinor:
    """
    Build Seurat reference dataset
    """
    input:
        mat = "reference/ellinor/fetch/matrix.mtx.gz",
        features = "reference/ellinor/fetch/features.tsv.gz",
        cells = "reference/ellinor/fetch/barcodes.tsv.gz",
        # mat_raw = "reference/ellinor/fetch/matrix_raw.mtx.gz",
        # features_raw = "reference/ellinor/fetch/features_raw.tsv.gz",
        # cells_raw = "reference/ellinor/fetch/barcodes_raw.tsv.gz",
        metadata = "reference/ellinor/fetch/metadata.tsv.gz"
    output:
        project_out = "reference/ellinor/seurat_build_reference/proj.rds",
        qc_violin = "reference/ellinor/seurat_build_reference/qc_violin.pdf",
        qc_scatter = "reference/ellinor/seurat_build_reference/qc_scatter.pdf",
        umap = "reference/ellinor/seurat_build_reference/umap.pdf"
    params:
        seed = config["seurat_seed"]
    log:
        console = "logs/reference/ellinor/seurat_build_reference/seurat_build_rna/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_build_reference_ellinor.R"

rule seurat_build_rna:
    """
    Build Seurat project
    """
    input:
        mat = "results/{sample}/fetch/matrix.mtx",
        features = "results/{sample}/fetch/features.tsv",
        cells = "results/{sample}/fetch/barcodes.tsv",
        metadata = "results/{sample}/atac/atac_qc.tsv"
    output:
        project_out = "results/{sample}/rna/seurat_build_rna/proj.rds",
        metadata = "results/{sample}/rna/seurat_build_rna/metadata.tsv",
        metadata_filtered = "results/{sample}/rna/seurat_build_rna/metadata_filtered.tsv",
        qc_violin = "results/{sample}/rna/seurat_build_rna/qc_violin.pdf",
        qc_scatter = "results/{sample}/rna/seurat_build_rna/qc_scatter.pdf",
    params:
        sample_name = lambda w: w.sample,
        seed = config["seurat_seed"],
        min_frags = config["archr_min_frags"],
        min_tss_enr = config["archr_min_tss_enr"],
        min_count_rna = config["seurat_min_count"],
        max_pct_mito_rna = config["seurat_max_pct_mito"]
    log:
        console = "logs/{sample}/rna/seurat_build_rna/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_build.R"

rule seurat_soupx_rna:
    """
    Filter ambient RNA with SoupX
    """
    input:
        mat = "results/{sample}/fetch/matrix.mtx",
        features = "results/{sample}/fetch/features.tsv",
        cells = "results/{sample}/fetch/barcodes.tsv",
        mat_raw = "results/{sample}/fetch/matrix_raw.mtx",
        features_raw = "results/{sample}/fetch/features_raw.tsv",
        cells_raw = "results/{sample}/fetch/barcodes_raw.tsv",
        metadata = "results/{sample}/rna/seurat_build_rna/metadata_filtered.tsv"
    output:
        project_out = "results/{sample}/rna/seurat_soupx_rna/proj.rds",
        umap = "results/{sample}/rna/seurat_soupx_rna/umap.pdf"
    params:
        sample_name = lambda w: w.sample,
        seed = config["seurat_seed"],
        min_count_rna = config["seurat_min_count"],
        max_pct_mito_rna = config["seurat_max_pct_mito"]
    log:
        console = "logs/{sample}/rna/seurat_soupx_rna/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_soupx_rna.R"

rule seurat_merge_soupx_plots:
    """
    Merge soupx plot pdf's
    """
    input:
        expand("results/{sample}/rna/seurat_soupx_rna/umap.pdf", sample=samples)
    output:
        "qc_all/seurat_soupx_umaps.pdf"
    conda:
        "../envs/fetch.yaml"
    shell:
        "pdfunite {input} {output}"

rule seurat_doublets_rna:
    """
    Filter doublets
    """
    input:
        project_in = "results/{sample}/rna/seurat_soupx_rna/proj.rds",
        doubletfinder_library_dir = "resources/doubletfinder_lib"
    output:
        project_out_all = "results/{sample}/rna/seurat_doublets_rna/proj_all.rds",
        project_out_filtered = "results/{sample}/rna/seurat_doublets_rna/proj_filtered.rds",
        metadata = "results/{sample}/rna/seurat_doublets_rna/metadata.tsv",
        umap = "results/{sample}/rna/seurat_doublets_rna/umap.pdf",
        umap_filtered = "results/{sample}/rna/seurat_doublets_rna/umap_filtered.pdf"
    params:
        seed = config["seurat_seed"],
        doublet_rate = config["doublet_formation_rate"],
        amulet_fdr = config["amulet_fdr"]
    log:
        console = "logs/{sample}/rna/seurat_doublets_rna/console.log"
    threads:
        config["max_threads_per_rule"]
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_doublets_rna.R"

rule seurat_merge_doublet_plots:
    """
    Merge doublet plot pdf's
    """
    input:
        expand("results/{sample}/rna/seurat_doublets_rna/umap.pdf", sample=samples)
    output:
        "qc_all/seurat_doublet_umaps.pdf"
    conda:
        "../envs/fetch.yaml"
    shell:
        "pdfunite {input} {output}"

rule seurat_transfer_rna:
    """
    Seurat RNA reference labeling
    """
    input:
        project_rna = "results/{sample}/rna/seurat_doublets_rna/proj_filtered.rds",
        project_ellinor = "reference/ellinor/seurat_build_reference/proj.rds",
        project_kramann = "reference/kramann/seurat_build_reference/proj.rds"
    output:
        project_out = "results/{sample}/rna/seurat_transfer_rna/proj.rds",
        umap_ellinor = "results/{sample}/rna/seurat_transfer_rna/umap_ellinor.pdf",
        umap_kramann = "results/{sample}/rna/seurat_transfer_rna/umap_kramann.pdf"
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/{sample}/rna/seurat_transfer_rna/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_transfer_rna.R"

rule seurat_merge_rna:
    """
    Merge RNA samples
    """
    input:
        projects_in = lambda w: f"results/{i}/rna/seurat_transfer_rna/proj.rds" for i in groups[w.group]
    output:
        project_out = "results_merged/{group}/rna/seurat_merge_rna/proj.rds",
        umap_dataset_pre_harmony = "results_merged/{group}/rna/seurat_merge_rna/umap_dataset_pre_harmony.pdf",
        umap_pre_harmony = "results_merged/{group}/rna/seurat_merge_rna/umap_pre_harmony.pdf",
        umap_dataset = "results_merged/{group}/rna/seurat_merge_rna/umap_dataset.pdf",
        umap = "results_merged/{group}/rna/seurat_merge_rna/umap.pdf",
    params:
        seed = config["seurat_seed"],
        samples = samples,
    log:
        console = "logs/merged/{group}/rna/seurat_merge_rna/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_merge_rna.R"

rule seurat_cluster_rna:
    """
    Seurat RNA clustering
    """
    input:
        project_in = "results_merged/rna/seurat_merge_rna/proj.rds",
        # project_ref = "reference/seurat_build_reference/proj.rds"
    output:
        project_out = "results_merged/rna/seurat_cluster_rna/proj.rds",
        umap = "results_merged/rna/seurat_cluster_rna/umap_clusters.pdf",
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/merged/rna/seurat_cluster_rna/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_cluster_rna.R"

rule seurat_name_rna:
    """
    Seurat RNA cluster naming
    """
    input:
        project_in = "results_merged/rna/seurat_cluster_rna/proj.rds",
    output:
        project_out = "results_merged/rna/seurat_name_rna/proj.rds",
        umap = "results_merged/rna/seurat_name_rna/umap_clusters.pdf",
        umap_qc = "results_merged/rna/seurat_name_rna/umap_qc.pdf",
        mat = "results_merged/rna/seurat_name_rna/confusion_mat.pdf",
        metadata = "results_merged/rna/seurat_name_rna/metadata.tsv",
    params:
        seed = config["seurat_seed"],
        # cluster_names = rna_cluster_names
    log:
        console = "logs/merged/rna/seurat_name_rna/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_name_rna.R"

rule seurat_write_embeddings:
    """
    Seurat save RNA embeddings
    """
    input:
        project_in = "results_merged/rna/seurat_name_rna/proj.rds",
    output:
        emb_coords = "results_merged/rna/seurat_write_embeddings/emb_coords.tsv",
        umap_coords = "results_merged/rna/seurat_write_embeddings/umap_coords.tsv",
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/merged/rna/seurat_write_embeddings/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_write_embeddings.R"

rule seurat_write_markers:
    """
    Seurat save RNA marker genes
    """
    input:
        project_in = "results_merged/rna/seurat_name_rna/proj.rds",
    output:
        markers = directory("results_merged/rna/seurat_write_markers/markers"),
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/merged/rna/seurat_write_markers/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_write_markers.R"

rule export_rna_embeddings:
    """
    Export RNA embeddings
    """
    input:
        emb = "results_merged/rna/seurat_write_embeddings/emb_coords.tsv",
        umap = "results_merged/rna/seurat_write_embeddings/umap_coords.tsv",
    output:
        emb = "export/rna/embeddings/harmony.tsv.gz",
        umap = "export/rna/embeddings/umap.tsv.gz",
    conda:
        "../envs/fetch.yaml"
    script:
        "../scripts/export_rna_embeddings.py"

rule export_rna_labels:
    """
    Export RNA cell types
    """
    input:
        "results_merged/rna/seurat_name_rna/metadata.tsv"
    output:
        "export/rna/labels/cell_types.tsv.gz",
    conda:
        "../envs/fetch.yaml"
    script:
        "../scripts/export_rna_labels.py"

rule export_rna_markers:
    """
    Export RNA markers
    """
    input:
        markers = "results_merged/rna/seurat_write_markers/markers",
        genes = expand("results/{sample}/fetch/features.tsv", sample=samples)
    output:
        directory("export/rna/markers")
    conda:
        "../envs/fetch.yaml"
    script:
        "../scripts/export_rna_markers.py"

rule export_rna_metadata:
    """
    Export RNA metadata
    """
    input:
        metadata = expand("results/{sample}/rna/seurat_build_rna/metadata.tsv", sample=samples),
        final_data = "results_merged/rna/seurat_name_rna/metadata.tsv"
    output:
        "export/rna/metadata.tsv.gz",
    conda:
        "../envs/fetch.yaml"
    script:
        "../scripts/export_rna_metadata.py"

rule export_rna_figures:
    """
    Export RNA figures
    """
    input:
        umap_labels = "results_merged/rna/seurat_name_rna/umap_clusters.pdf",
        umap_samples = "results_merged/rna/seurat_merge_rna/umap_dataset.pdf"
    output:
        scratch = directory("results_merged/rna/export_figures"),
        tarball = "export/rna/figures.tar.gz"
    params:
        readme = workflow.source_path("../resources/rna_figures_readme.txt")
    conda:
        "../envs/fetch.yaml"
    shell:
        "mkdir -p {output.scratch}; "
        "cp {input.umap_labels} {output.scratch}/umap_labels.pdf; "
        "cp {input.umap_samples} {output.scratch}/umap_samples.pdf; "
        "cp {params.readme} {output.scratch}/README.txt; "
        "tar -zcvf {output.tarball} {output.scratch}"

rule export_rna_dataset_names:
    """
    Export RNA dataset names used in analysis
    """
    output:
        "export/rna/datasets.txt"
    params:
        datasets = samples
    conda:
        "../envs/fetch.yaml"
    shell:
        "echo {params.datasets} | tr ' ' '\\n' > {output}"