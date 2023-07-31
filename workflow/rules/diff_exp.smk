rule pseudobulk_samples:
    """
    Build RNA pseudobulks for differential expression analyses
    """
    input:
        project_in = "results_groups/{group}/rna/seurat_embed_all/proj.rds"
    output:
        out_mats = [f"results_groups/{{group}}/rna/pseudobulk_samples/label/{label}/mat.tsv" for label in config["l1_labels"]],
        out_metadata = "results_groups/{group}/rna/pseudobulk_samples/metadata.tsv",
    params:
        seed = config["seurat_seed"],
        cell_types = config["l1_labels"],
    log:
        console = "logs/merged/{group}/rna/pseudobulk_samples/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/pseudobulk_samples.R"

rule collate_pseudobulks:
    """
    Collate pseudobulks across groups
    """
    input:
        mats = expand("results_groups/{group}/rna/pseudobulk_samples/label/{label}/mat.tsv", group=group_names, allow_missing=True),
        metadata = expand("results_groups/{group}/rna/pseudobulk_samples/metadata.tsv", group=group_names)
    output:
        mat = "results_diff_exp/label/{label}/collate_pseudobulks/mat.tsv",
        metadata = "results_diff_exp/label/{label}/collate_pseudobulks/metadata.tsv",
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/diff_exp/label/{label}/collate_pseudobulks/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/collate_pseudobulks.R"

rule run_deseq:
    """
    Run DESeq2 for differential expression analyses
    """
    input:
        mat = "results_diff_exp/label/{label}/collate_pseudobulks/mat.tsv",
        metadata = "results_diff_exp/label/{label}/collate_pseudobulks/metadata.tsv",
    output:
        dispersion_plot = "results_diff_exp/label/{label}/run_deseq/dispersion_estimates.pdf",
        coefficients = directory("results_diff_exp/label/{label}/run_deseq/coefficients")
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/diff_exp/label/{label}/run_deseq/console.log"
    conda:
        "../envs/deseq.yaml"
    script:
        "../scripts/run_deseq.R"

rule run_deseq_sva:
    """
    Run DESeq2 with SVA for differential expression analyses
    """
    input:
        mat = "results_diff_exp/label/{label}/collate_pseudobulks/mat.tsv",
        metadata = "results_diff_exp/label/{label}/collate_pseudobulks/metadata.tsv",
    output:
        dispersion_plot = "results_diff_exp/label/{label}/run_deseq_sva/dispersion_estimates.pdf",
        sv_plot = "results_diff_exp/label/{label}/run_deseq_sva/sv.pdf",
        coefficients = directory("results_diff_exp/label/{label}/run_deseq_sva/coefficients")
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/diff_exp/label/{label}/run_deseq_sva/console.log"
    conda:
        "../envs/deseq.yaml"
    script:
        "../scripts/run_deseq_sva.R"