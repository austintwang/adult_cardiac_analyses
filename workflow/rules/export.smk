rule export_l1_metadata: #
    """
    Export metadata
    """
    input:
        metadata_rna = lambda w: [f"results/{sample}/rna/seurat_build_rna_strict/metadata.tsv" for sample in groups[w.group]],
        metadata_atac = lambda w: [f"results/{sample}/atac/atac_qc_extended.tsv.gz" for sample in groups[w.group]],
        dataset_data = lambda w: [f"results/{sample}/fetch/dataset_info.json" for sample in groups[w.group]],
        final_data = "results_groups/{group}/rna/seurat_write_embeddings/metadata.tsv",
    output:
        metadata = "export_l1/{group}/metadata.tsv.gz",
        datasets = "export_l1/{group}/datasets.txt"
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
        "export_l1/{group}/labels/cell_types_l1.tsv.gz",
    conda:
        "../envs/fetch.yaml"
    script:
        "../scripts/export_l1_labels.py"

rule export_l1_markers:
    """
    Export level-1 cell type markers
    """
    input:
        markers = "results_groups/{group}/rna/seurat_write_l1_markers/markers",
        genes = lambda w: [f"results/{sample}/fetch/features.tsv" for sample in groups[w.group]]
    output:
        directory("export_l1/{group}/markers")
    conda:
        "../envs/fetch.yaml"
    script:
        "../scripts/export_l1_markers.py"

rule export_l1_embeddings: #
    """
    Export cell embeddings
    """
    input:
        rna_pca = "results_groups/{group}/rna/seurat_write_embeddings/rna_pca.tsv",
        rna_harmony = "results_groups/{group}/rna/seurat_write_embeddings/rna_harmony.tsv",
        rna_umap = "results_groups/{group}/rna/seurat_write_embeddings/rna_umap.tsv",
        atac_lsi = "results_groups/{group}/rna/seurat_write_embeddings/atac_lsi.tsv",
        atac_harmony = "results_groups/{group}/rna/seurat_write_embeddings/atac_harmony.tsv",
        atac_umap = "results_groups/{group}/rna/seurat_write_embeddings/atac_umap.tsv"
    output:
        rna_pca = "export_l1/{group}/embeddings/rna_pca.tsv.gz",
        rna_harmony = "export_l1/{group}/embeddings/rna_harmony.tsv.gz",
        rna_umap = "export_l1/{group}/embeddings/rna_umap.tsv.gz",
        atac_lsi = "export_l1/{group}/embeddings/atac_lsi.tsv.gz",
        atac_harmony = "export_l1/{group}/embeddings/atac_harmony.tsv.gz",
        atac_umap = "export_l1/{group}/embeddings/atac_umap.tsv.gz"
    conda:
        "../envs/fetch.yaml"
    script:
        "../scripts/export_embeddings.py"

rule export_l1_figures: #
    """
    Export figures
    """
    input:
        umap_rna_dataset = "results_groups/{group}/rna/seurat_embed_all/umap_rna_dataset.pdf",
        umap_rna_cell_types = "results_groups/{group}/rna/seurat_embed_all/umap_rna_cell_types.pdf",
        umap_atac_dataset = "results_groups/{group}/rna/seurat_embed_all/umap_rna_dataset.pdf",
        umap_atac_cell_types = "results_groups/{group}/rna/seurat_embed_all/umap_rna_cell_types.pdf",
    output:
        scratch = directory("results_groups/{group}/rna/export_figures"),
        tarball = "export_l1/{group}/figures.tar.gz"
    params:
        readme = workflow.source_path("../resources/figures_readme.txt")
    conda:
        "../envs/fetch.yaml"
    shell:
        "mkdir -p {output.scratch}; "
        "cp {input.umap_rna_cell_types} {output.scratch}/rna_umap_labels.pdf; "
        "cp {input.umap_rna_dataset} {output.scratch}/rna_umap_samples.pdf; "
        "cp {input.umap_atac_cell_types} {output.scratch}/atac_umap_labels.pdf; "
        "cp {input.umap_atac_dataset} {output.scratch}/atac_umap_samples.pdf; "
        "cp {params.readme} {output.scratch}/README.txt; "
        "tar -zcvf {output.tarball} {output.scratch}"