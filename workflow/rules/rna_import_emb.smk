def get_seq_emb_cells(w):
    region, status, sex = w.group.split("-")
    region = region.upper()
    sex = sex.upper()
    return os.path.join(config["seq_emb_dir"], f"{status}_{region}_{sex}_cell_labels.txt")

def get_seq_emb_scbasset(w):
    region, status, sex = w.group.split("-")
    region = region.upper()
    sex = sex.upper()
    return os.path.join(config["seq_emb_dir"], f"{status}_{region}_{sex}_embeddings_scBasset_preharmony.mtx")

def get_seq_emb_cellspace(w):
    region, status, sex = w.group.split("-")
    region = region.upper()
    sex = sex.upper()
    return os.path.join(config["seq_emb_dir"], f"{status}_{region}_{sex}_embeddings_CellSpace_preharmony.mtx")


rule seurat_add_seq_emb:
    """
    Add ATAC sequence embeddings from Salil and harmonize
    """
    input:
        project_in = "results_groups/{group}/rna/seurat_subclusters_supergroups_to_groups/proj.rds",
        cells = get_seq_emb_cells,
        scbasset = get_seq_emb_scbasset,
        cellspace = get_seq_emb_cellspace
    output:
        project_out = "results_groups/{group}/rna/seurat_add_seq_emb/proj.rds",
        rna_test_harmony_mat = "results_groups/{group}/rna/seurat_add_seq_emb/rna_test_harmony_mat.mtx",
        rna_test_harmony_rows = "results_groups/{group}/rna/seurat_add_seq_emb/rna_test_harmony_rows.txt",
        rna_test_harmony_cols = "results_groups/{group}/rna/seurat_add_seq_emb/rna_test_harmony_cols.txt",
        scbasset_harmony_mat = "results_groups/{group}/rna/seurat_add_seq_emb/scbasset_harmony_mat.mtx",
        scbasset_harmony_rows = "results_groups/{group}/rna/seurat_add_seq_emb/scbasset_harmony_rows.txt",
        scbasset_harmony_cols = "results_groups/{group}/rna/seurat_add_seq_emb/scbasset_harmony_cols.txt",
        cellspace_harmony_mat = "results_groups/{group}/rna/seurat_add_seq_emb/cellspace_harmony_mat.mtx",
        cellspace_harmony_rows = "results_groups/{group}/rna/seurat_add_seq_emb/cellspace_harmony_rows.txt",
        cellspace_harmony_cols = "results_groups/{group}/rna/seurat_add_seq_emb/cellspace_harmony_cols.txt",
        umap_scbasset_l1 = "results_groups/{group}/rna/seurat_add_seq_emb/umap_scbasset_l1.pdf",
        umap_cellspace_l1 = "results_groups/{group}/rna/seurat_add_seq_emb/umap_cellspace_l1.pdf",
        umap_scbasset_l2 = "results_groups/{group}/rna/seurat_add_seq_emb/umap_scbasset_l2.pdf",
        umap_cellspace_l2 = "results_groups/{group}/rna/seurat_add_seq_emb/umap_cellspace_l2.pdf",
        metadata = "results_groups/{group}/rna/seurat_add_seq_emb/metadata.tsv"
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/merged/{group}/rna/seurat_add_seq_emb/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_import_seq_embeddings.R"

rule seurat_write_rna_all:
    """
    Write RNA all counts embeddings
    """
    input:
        project_in = "results_groups/{group}/rna/seurat_embed_all/proj.rds"
        # project_in = "results_groups/{group}/rna/seurat_clusters_to_groups/proj.rds"
    output:
        rna_harmony_mat = "results_groups/{group}/rna/seurat_write_rna_all/rna_harmony_mat.mtx",
        rna_harmony_rows = "results_groups/{group}/rna/seurat_write_rna_all/rna_harmony_rows.txt",
        rna_harmony_cols = "results_groups/{group}/rna/seurat_write_rna_all/rna_harmony_cols.txt",
        rna_pca_mat = "results_groups/{group}/rna/seurat_write_rna_all/rna_pca_mat.mtx",
        rna_pca_rows = "results_groups/{group}/rna/seurat_write_rna_all/rna_pca_rows.txt",
        rna_pca_cols = "results_groups/{group}/rna/seurat_write_rna_all/rna_pca_cols.txt",
        metadata = "results_groups/{group}/rna/seurat_write_rna_all/metadata.tsv"
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/merged/{group}/rna/seurat_write_rna_all/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_write_rna_all.R"


def get_atac_emb_cells(w):
    region, status, sex = w.group.split("-")
    region = region.upper()
    sex = sex.upper()
    return os.path.join(config["atac_emb_dir"], f"{status}_{region}_{sex}_rows.txt")

def get_atac_emb(w):
    region, status, sex = w.group.split("-")
    region = region.upper()
    sex = sex.upper()
    return os.path.join(config["atac_emb_dir"], f"{status}_{region}_{sex}_preh_mat.mtx")

rule seurat_import_atac_embeddings:
    """
    Add ATAC embeddings from Salil and harmonize
    """
    input:
        project_in = "results_groups/{group}/rna/seurat_embed_all/proj.rds",
        cells = get_atac_emb_cells,
        atac = get_atac_emb
    output:
        project_out = "results_groups/{group}/rna/seurat_import_atac_embeddings/proj.rds",
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/merged/{group}/rna/seurat_import_atac_embeddings/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_import_atac_embeddings.R"
