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
        project_in = "results_groups/{group}/rna/seurat_import_atac_embeddings/proj.rds"
    output:
        metadata = "results_groups/{group}/rna/seurat_write_embeddings/metadata.tsv",
        rna_pca = "results_groups/{group}/rna/seurat_write_embeddings/rna_pca.tsv",
        rna_harmony = "results_groups/{group}/rna/seurat_write_embeddings/rna_harmony.tsv",
        rna_umap = "results_groups/{group}/rna/seurat_write_embeddings/rna_umap.tsv",
        atac_lsi = "results_groups/{group}/rna/seurat_write_embeddings/atac_lsi.tsv",
        atac_harmony = "results_groups/{group}/rna/seurat_write_embeddings/atac_harmony.tsv",
        atac_umap = "results_groups/{group}/rna/seurat_write_embeddings/atac_umap.tsv",
        umap_rna_dataset = "results_groups/{group}/rna/seurat_embed_all/umap_rna_dataset.pdf",
        umap_rna_cell_types = "results_groups/{group}/rna/seurat_embed_all/umap_rna_cell_types.pdf",
        umap_atac_dataset = "results_groups/{group}/rna/seurat_embed_all/umap_atac_dataset.pdf",
        umap_atac_cell_types = "results_groups/{group}/rna/seurat_embed_all/umap_atac_cell_types.pdf",
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/merged/{group}/rna/seurat_write_embeddings/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_write_embeddings.R"