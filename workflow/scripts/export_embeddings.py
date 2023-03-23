import gzip

HEADER_PCA_HARMONY = """
# Non-batch-corrected PCA cell embeddings
# cell_id: Cell ID used in integrated analysis
# harmony_1, harmony_2, … harmony_50: Columns of a 50-dimensional Harmony embedding vector
"""

HEADER_RNA_HARMONY = """
# Unified cell embeddings integrated using Harmony
# cell_id: Cell ID used in integrated analysis
# PC_1, PC_2, … PC_50: Columns of a 50-dimensional PCA embedding vector
"""

HEADER_RNA_UMAP = """
# UMAP coordinates of each cell
# cell_id: Cell ID used in integrated analysis
# UMAP_1, UMAP_2: UMAP x and y coordinates, respectively
"""

def main(rna_pca_in_path, rna_harmony_in_path, rna_umap_in_path, rna_pca_out_path, rna_harmony_out_path, rna_umap_out_path):
    with open(rna_pca_in_path) as f, gzip.open(rna_pca_out_path, "wt") as fo:
        fo.write(HEADER_PCA_HARMONY)
        h = f.readline()
        fo.write("cell_id" + h)
        for line in f:
            fo.write(line)

    with open(rna_harmony_in_path) as f, gzip.open(rna_harmony_out_path, "wt") as fo:
        fo.write(HEADER_PCA_HARMONY)
        h = f.readline()
        fo.write("cell_id" + h)
        for line in f:
            fo.write(line)


    with open(rna_umap_in_path) as f, gzip.open(rna_umap_out_path, "wt") as fo:
        fo.write(HEADER_RNA_UMAP)
        h = f.readline()
        fo.write("cell_id" + h)
        for line in f:
            fo.write(line)

rna_pca_in_path = snakemake.input["rna_pca"]
rna_harmony_in_path = snakemake.input["rna_harmony"]
rna_umap_in_path = snakemake.input["rna_umap"]

rna_pca_out_path = snakemake.output["rna_pca"]
rna_harmony_out_path = snakemake.output["rna_harmony"]
rna_umap_out_path = snakemake.output["rna_umap"]

main(rna_pca_in_path, rna_harmony_in_path, rna_umap_in_path, rna_pca_out_path, rna_harmony_out_path, rna_umap_out_path)