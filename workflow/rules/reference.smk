rule seurat_build_reference_kramann:
    """
    Build Kramann reference dataset
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
    Build Ellinor reference dataset
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

rule seurat_transfer_kramann_rna:
    """
    Seurat RNA reference labeling
    """
    input:
        project_rna = "results/{sample}/rna/seurat_doublets_rna/proj_filtered.rds",
        project_kramann = "reference/kramann/seurat_build_reference/proj.rds"
    output:
        data_out = "results/{sample}/rna/seurat_transfer_rna/labels_kramann.tsv",
        umap_kramann = "results/{sample}/rna/seurat_transfer_rna/umap_kramann.pdf"
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/{sample}/rna/seurat_transfer_rna/kramann.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_transfer_kramann_rna.R"

rule seurat_transfer_ellinor_rna:
    """
    Seurat RNA reference labeling
    """
    input:
        project_rna = "results/{sample}/rna/seurat_doublets_rna/proj_filtered.rds",
        project_ellinor = "reference/ellinor/seurat_build_reference/proj.rds",
    output:
        data_out = "results/{sample}/rna/seurat_transfer_rna/labels_ellinor.tsv",
        umap_ellinor_coarse = "results/{sample}/rna/seurat_transfer_rna/umap_ellinor_coarse.pdf",
        umap_ellinor_fine = "results/{sample}/rna/seurat_transfer_rna/umap_ellinor_fine.pdf",
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/{sample}/rna/seurat_transfer_rna/ellinor.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_transfer_ellinor_rna.R"

rule seurat_transfer_azimuth_rna:
    """
    Seurat RNA reference labeling
    """
    input:
        project_rna = "results/{sample}/rna/seurat_doublets_rna/proj_filtered.rds"
    output:
        data_out = "results/{sample}/rna/seurat_transfer_rna/labels_azimuth.tsv",
        umap_azimuth = "results/{sample}/rna/seurat_transfer_rna/umap_azimuth.pdf",
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/{sample}/rna/seurat_transfer_rna/ellinor.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_transfer_azimuth_rna.R"

rule seurat_merge_kramann_plots:
    """
    Merge kramann projection plot pdf's
    """
    input:
        expand("results/{sample}/rna/seurat_transfer_rna/umap_kramann.pdf", sample=samples)
    output:
        "qc_all/seurat_kramann_umaps.pdf"
    conda:
        "../envs/fetch.yaml"
    shell:
        "pdfunite {input} {output}"

rule seurat_merge_ellinor_plots:
    """
    Merge ellinor projection plot pdf's
    """
    input:
        coarse = expand("results/{sample}/rna/seurat_transfer_rna/umap_ellinor_coarse.pdf", sample=samples),
        fine = expand("results/{sample}/rna/seurat_transfer_rna/umap_ellinor_fine.pdf", sample=samples)
    output:
        coarse = "qc_all/seurat_ellinor_coarse_umaps.pdf",
        fine = "qc_all/seurat_ellinor_fine_umaps.pdf"
    conda:
        "../envs/fetch.yaml"
    shell:
        "pdfunite {input.coarse} {output.coarse}; "
        "pdfunite {input.fine} {output.fine}"

