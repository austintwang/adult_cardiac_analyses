rule export_metadata: #
    """
    Export metadata
    """
    input:
        metadata_rna = lambda w: [f"results/{sample}/rna/seurat_build_rna_strict/metadata.tsv" for sample in groups[w.group]],
        metadata_atac = lambda w: [f"results/{sample}/atac/atac_qc_extended.tsv.gz" for sample in groups[w.group]],
        final_data = "results_groups/{group}/rna/seurat_write_rna_all/metadata.tsv"
    output:
        metadata = "export/{group}/metadata.tsv.gz",
        datasets = "export/{group}/datasets.txt"
    conda:
        "../envs/fetch.yaml"
    script:
        "../scripts/export_metadata.py"

rule export_l1_labels:
    """
    Export level-1 cell types
    """
    input:
        "results_groups/{group}/rna/seurat_write_rna_all/metadata.tsv"
    output:
        "export/rna/labels/cell_types.tsv.gz",
    conda:
        "../envs/fetch.yaml"
    script:
        "../scripts/export_rna_labels.py"

rule export_l1_markers:
    """
    Export level-1 cell type markers
    """
    input:
        markers = "results_groups/{group}/rna/seurat_write_l1_markers/markers"",
        genes = lambda w: [f"results/{sample}/fetch/features.tsv" for sample in groups[w.group]]
    output:
        directory("export/{group}/markers")
    conda:
        "../envs/fetch.yaml"
    script:
        "../scripts/export_l1_markers.py"

rule export_embeddings: #
    """
    Export cell embeddings
    """
    input:
        rna_pca = "results_merged/{group}/rna/seurat_write_embeddings/rna_pca.tsv",
        rna_harmony = "results_merged/{group}/rna/seurat_write_embeddings/rna_harmony.tsv",
        rna_umap = "results_merged/{group}/rna/seurat_write_embeddings/rna_umap.tsv"
    output:
        rna_pca = "export/{group}/embeddings/rna_pca.tsv.gz",
        rna_harmony = "export/{group}/embeddings/rna_harmony.tsv.gz",
        rna_umap = "export/{group}/embeddings/rna_umap.tsv.gz"
    conda:
        "../envs/fetch.yaml"
    script:
        "../scripts/export_rna_embeddings.py"

rule export_figures: #
    """
    Export figures
    """
    input:
        rna_umap_labels = "results_groups/{group}/rna/seurat_embed_all/umap_cell_types.pdf",
        rna_umap_samples = "results_groups/{group}/rna/seurat_embed_all/umap_dataset.pdf"
    output:
        scratch = directory("results_merged/rna/export_figures"),
        tarball = "export/{group}/figures.tar.gz"
    params:
        readme = workflow.source_path("../resources/figures_readme.txt")
    conda:
        "../envs/fetch.yaml"
    shell:
        "mkdir -p {output.scratch}; "
        "cp {input.rna_umap_labels} {output.scratch}/rna_umap_labels.pdf; "
        "cp {input.rna_umap_samples} {output.scratch}/rna_umap_samples.pdf; "
        "cp {params.readme} {output.scratch}/README.txt; "
        "tar -zcvf {output.tarball} {output.scratch}"