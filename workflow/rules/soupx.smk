rule seurat_soupx_filtered:
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
        metadata = lambda w: f"results_groups/{sample_to_group[w.sample]}/rna/seurat_write_embeddings/metadata.tsv"
    output:
        project_out = "results/{sample}/rna/seurat_soupx_filtered/proj.rds",
        rho = "results/{sample}/rna/seurat_soupx_filtered/rho.pdf"
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/{sample}/rna/seurat_soupx_filtered/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_soupx_filtered.R"

rule seurat_merge_soupx_plots_filtered:
    """
    Merge soupx plot pdf's
    """
    input:
        expand("results/{sample}/rna/seurat_soupx_filtered/rho.pdf", sample=samples)
    output:
        "qc_all/seurat_soupx_filtered.pdf"
    conda:
        "../envs/fetch.yaml"
    shell:
        "pdfunite {input} {output}"

rule seurat_write_rna_soupx:
    """
    Write SoupX-corrected count matrices
    """
    input:
        project_in = "results/{sample}/rna/seurat_soupx_filtered/proj.rds"
    output:
        counts_mat = "results/{sample}/rna/seurat_soupx_filtered/adjusted_counts/counts_mat.mtx",
        counts_rows = "results/{sample}/rna/seurat_soupx_filtered/adjusted_counts/counts_rows.txt",
        counts_cols = "results/{sample}/rna/seurat_soupx_filtered/adjusted_counts/counts_cols.txt",
    params:
        seed = config["seurat_seed"],
    log:
        console =  "logs/{sample}/rna/seurat_write_rna_soupx/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_write_rna_soupx.R"
