import gzip

HEADER_RNA_PCA = """
# Non-batch-corrected RNA PCA cell embeddings
# cell_id: Cell ID used in integrated analysis
# PC_1, PC_2, … PC_50: Columns of a 50-dimensional PCA embedding vector
"""

HEADER_RNA_HARMONY = """
# Unified RNA cell embeddings integrated using Harmony
# cell_id: Cell ID used in integrated analysis
# harmony_1, harmony_2, … harmony_50: Columns of a 50-dimensional Harmony embedding vector
"""

HEADER_RNA_UMAP = """
# RNA UMAP coordinates of each cell
# cell_id: Cell ID used in integrated analysis
# UMAP_1, UMAP_2: UMAP x and y coordinates, respectively
"""

HEADER_ATAC_LSI = """
# Non-batch-corrected ATAC LSI cell embeddings
# cell_id: Cell ID used in integrated analysis
# LSI_1, LSI_2, … LSI_98: Columns of a 98-dimensional LSI embedding vector
"""

HEADER_ATAC_HARMONY = """
# Unified ATAC cell embeddings integrated using Harmony
# cell_id: Cell ID used in integrated analysis
# harmony_1, harmony_2, … harmony_98: Columns of a 98-dimensional Harmony embedding vector
"""

HEADER_ATAC_UMAP = """
# ATAC UMAP coordinates of each cell
# cell_id: Cell ID used in integrated analysis
# UMAP_1, UMAP_2: UMAP x and y coordinates, respectively
"""

def export(in_path, out_path, header):
    with open(in_path) as f, gzip.open(out_path, "wt") as fo:
        fo.write(header)
        h = f.readline()
        fo.write("cell_id" + h)
        for line in f:
            fo.write(line)

def main(rna_pca_in_path, rna_harmony_in_path, rna_umap_in_path, atac_lsi_in_path, atac_harmony_in_path, atac_umap_in_path, rna_pca_out_path, rna_harmony_out_path, rna_umap_out_path, atac_lsi_out_path, atac_harmony_out_path, atac_umap_out_path):
    export(rna_pca_in_path, rna_pca_out_path, HEADER_RNA_PCA)
    export(rna_harmony_in_path, rna_harmony_out_path, HEADER_RNA_HARMONY)
    export(rna_umap_in_path, rna_umap_out_path, HEADER_RNA_UMAP)

    export(atac_lsi_in_path, atac_lsi_out_path, HEADER_RNA_PCA)
    export(atac_harmony_in_path, atac_harmony_out_path, HEADER_RNA_HARMONY)
    export(atac_umap_in_path, atac_umap_out_path, HEADER_RNA_UMAP)

rna_pca_in_path = snakemake.input["rna_pca"]
rna_harmony_in_path = snakemake.input["rna_harmony"]
rna_umap_in_path = snakemake.input["rna_umap"]
atac_lsi_in_path = snakemake.input["atac_lsi"]
atac_harmony_in_path = snakemake.input["atac_harmony"]
atac_umap_in_path = snakemake.input["atac_umap"]

rna_pca_out_path = snakemake.output["rna_pca"]
rna_harmony_out_path = snakemake.output["rna_harmony"]
rna_umap_out_path = snakemake.output["rna_umap"]
atac_lsi_out_path = snakemake.output["atac_lsi"]
atac_harmony_out_path = snakemake.output["atac_harmony"]
atac_umap_out_path = snakemake.output["atac_umap"]

main(rna_pca_in_path, rna_harmony_in_path, rna_umap_in_path, atac_lsi_in_path, atac_harmony_in_path, atac_umap_in_path, rna_pca_out_path, rna_harmony_out_path, rna_umap_out_path, atac_lsi_out_path, atac_harmony_out_path, atac_umap_out_path)