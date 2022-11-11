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

rule seurat_build_rna_strict:
    """
    Build Seurat project (stringent cutoffs)
    """
    input:
        mat = "results/{sample}/fetch/matrix.mtx",
        features = "results/{sample}/fetch/features.tsv",
        cells = "results/{sample}/fetch/barcodes.tsv",
        metadata = "results/{sample}/atac/atac_qc.tsv"
    output:
        project_out = "results/{sample}/rna/seurat_build_rna_strict/proj.rds",
        metadata = "results/{sample}/rna/seurat_build_rna_strict/metadata.tsv",
        metadata_filtered = "results/{sample}/rna/seurat_build_rna_strict/metadata_filtered.tsv",
        qc_violin = "results/{sample}/rna/seurat_build_rna_strict/qc_violin.pdf",
        qc_scatter = "results/{sample}/rna/seurat_build_rna_strict/qc_scatter.pdf",
    params:
        sample_name = lambda w: w.sample,
        seed = config["seurat_seed"],
        min_frags = config["archr_min_frags"],
        min_tss_enr = config["archr_min_tss_enr"],
        min_count_rna = config["seurat_min_count"],
        max_pct_mito_rna = config["seurat_max_pct_mito"]
    log:
        console = "logs/{sample}/rna/seurat_build_rna_strict/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_build.R"

rule seurat_doublets_no_soupx:
    """
    Filter doublets (no ambient RNA removal)
    """
    input:
        project_in = "results/{sample}/rna/seurat_build_rna_strict/proj.rds",
        doubletfinder_library_dir = "resources/doubletfinder_lib"
    output:
        project_out_all = "results/{sample}/rna/seurat_doublets_no_soupx/proj_all.rds",
        project_out_filtered = "results/{sample}/rna/seurat_doublets_no_soupx/proj_filtered.rds",
        metadata = "results/{sample}/rna/seurat_doublets_no_soupx/metadata.tsv",
        umap = "results/{sample}/rna/seurat_doublets_no_soupx/umap.pdf",
        umap_filtered = "results/{sample}/rna/seurat_doublets_no_soupx/umap_filtered.pdf"
    params:
        seed = config["seurat_seed"],
        doublet_rate = config["doublet_formation_rate"],
        amulet_fdr = config["amulet_fdr"]
    log:
        console = "logs/{sample}/rna/seurat_doublets_no_soupx/console.log"
    threads:
        config["max_threads_per_rule"]
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_doublets_rna.R"

rule seurat_merge_doublet_plots_2:
    """
    Merge doublet plot pdf's (no ambient RNA removal)
    """
    input:
        expand("results/{sample}/rna/seurat_doublets_no_soupx/umap.pdf", sample=samples)
    output:
        "qc_all/seurat_doublet_no_soupx_umaps.pdf"
    conda:
        "../envs/fetch.yaml"
    shell:
        "pdfunite {input} {output}"

rule seurat_countsplit:
    """
    Create train and test pseudoreplicates
    """
    input:
        project_in = "results/{sample}/rna/seurat_doublets_no_soupx/proj_filtered.rds"
    output:
        project_out = "results/{sample}/rna/seurat_countsplit/proj.rds",
        umap = "results/{sample}/rna/seurat_countsplit/umap.pdf",
    params:
        seed = config["seurat_seed"],
        split_frac = config["pseudorep_split_frac"]
    log:
        console = "logs/{sample}/rna/seurat_countsplit/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_countsplit.R"

# rule seurat_integrate_rna:
#     """
#     Integrate RNA samples using Harmony
#     """
#     input:
#         projects_in = lambda w: [f"results/{i}/rna/seurat_doublets_rna/proj_filtered.rds" for i in groups[w.group]]
#     output:
#         project_out = "results_merged/{group}/rna/seurat_integrate_rna/proj.rds",
#         umap_dataset_pre_harmony = "results_merged/{group}/rna/seurat_integrate_rna/umap_dataset_pre_harmony.pdf",
#         umap_mixing_pre_harmony = "results_merged/{group}/rna/seurat_integrate_rna/umap_mixing_pre_harmony.pdf",
#         umap_dataset_harmony = "results_merged/{group}/rna/seurat_integrate_rna/umap_dataset_harmony.pdf",
#         umap_mixing_harmony = "results_merged/{group}/rna/seurat_integrate_rna/umap_mixing_harmony.pdf",
#     params:
#         seed = config["seurat_seed"],
#         samples = lambda w: groups[w.group],
#     log:
#         console = "logs/merged/{group}/rna/seurat_integrate_rna/console.log"
#     conda:
#         "../envs/seurat.yaml"
#     script:
#         "../scripts/seurat_integrate_rna.R"

# rule seurat_integrate_qc_extras:
#     """
#     Additional integration QC plots
#     """
#     input:
#         project_in = "results_merged/{group}/rna/seurat_integrate_rna/proj.rds"
#     output:
#         umap_counts = "results_merged/{group}/rna/seurat_integrate_rna/umap_counts.pdf",
#         umap_doubletfinder = "results_merged/{group}/rna/seurat_integrate_rna/umap_doubletfinder.pdf",
#         umap_amulet = "results_merged/{group}/rna/seurat_integrate_rna/umap_amulet.pdf"
#     params:
#         seed = config["seurat_seed"],
#         samples = lambda w: groups[w.group],
#     log:
#         console = "logs/merged/{group}/rna/seurat_integrate_rna_extras/console.log"
#     conda:
#         "../envs/seurat.yaml"
#     script:
#         "../scripts/seurat_integration_qc.R"

# rule seurat_merge_integration_plots:
#     """
#     Merge integration QC pdf's
#     """
#     input:
#         "results_merged/{group}/rna/seurat_integrate_rna/umap_dataset_pre_harmony.pdf",
#         "results_merged/{group}/rna/seurat_integrate_rna/umap_dataset_harmony.pdf",
#         "results_merged/{group}/rna/seurat_integrate_rna/umap_mixing_pre_harmony.pdf",
#         "results_merged/{group}/rna/seurat_integrate_rna/umap_mixing_harmony.pdf",
#         "results_merged/{group}/rna/seurat_integrate_rna/umap_counts.pdf",
#         "results_merged/{group}/rna/seurat_integrate_rna/umap_doubletfinder.pdf",
#         "results_merged/{group}/rna/seurat_integrate_rna/umap_amulet.pdf"
#     output:
#         "results_merged/{group}/rna/seurat_integrate_rna/integration_plots.pdf"
#     conda:
#         "../envs/fetch.yaml"
#     shell:
#         "pdfunite {input} {output}; "

# rule seurat_load_transfer_labels:
#     """
#     Add labels to integrated analyses
#     """
#     input:
#         project_integrated = "results_merged/{group}/rna/seurat_integrate_rna/proj.rds",
#         tables = expand(
#             "results_merged/{group}/rna/seurat_integrate_transfers/{reference}.tsv", 
#             reference = ["kramann", "ellinor", "teichmann", "azimuth"], 
#             allow_missing=True
#         )
#     output:
#         project_out = "results_merged/{group}/rna/seurat_load_transfer_labels/proj.rds",
#         umap_kramann_coarse = "results_merged/{group}/rna/seurat_load_transfer_labels/umap_kramann_coarse.pdf",
#         umap_kramann_fine = "results_merged/{group}/rna/seurat_load_transfer_labels/umap_kramann_fine.pdf",
#         umap_ellinor_coarse = "results_merged/{group}/rna/seurat_load_transfer_labels/umap_ellinor_coarse.pdf",
#         umap_ellinor_fine = "results_merged/{group}/rna/seurat_load_transfer_labels/umap_ellinor_fine.pdf",
#         umap_teichmann_coarse = "results_merged/{group}/rna/seurat_load_transfer_labels/umap_teichmann_coarse.pdf",
#         umap_teichmann_fine = "results_merged/{group}/rna/seurat_load_transfer_labels/umap_teichmann_fine.pdf",
#         umap_azimuth_coarse = "results_merged/{group}/rna/seurat_load_transfer_labels/umap_azimuth_coarse.pdf",
#         umap_azimuth_fine = "results_merged/{group}/rna/seurat_load_transfer_labels/umap_azimuth_fine.pdf",
#     params:
#         seed = config["seurat_seed"],
#     log:
#         console = "logs/merged/{group}/rna/seurat_load_transfer_labels/console.log"
#     conda:
#         "../envs/seurat.yaml"
#     script:
#         "../scripts/seurat_load_transfer_labels.R"

# rule seurat_merge_reference_label_plots:
#     """
#     Merge reference projection plot pdf's
#     """
#     input:
#         "results_merged/{group}/rna/seurat_load_transfer_labels/umap_kramann_coarse.pdf",
#         "results_merged/{group}/rna/seurat_load_transfer_labels/umap_kramann_fine.pdf",
#         "results_merged/{group}/rna/seurat_load_transfer_labels/umap_ellinor_coarse.pdf",
#         "results_merged/{group}/rna/seurat_load_transfer_labels/umap_ellinor_fine.pdf",
#         "results_merged/{group}/rna/seurat_load_transfer_labels/umap_teichmann_coarse.pdf",
#         "results_merged/{group}/rna/seurat_load_transfer_labels/umap_teichmann_fine.pdf",
#         "results_merged/{group}/rna/seurat_load_transfer_labels/umap_azimuth_coarse.pdf",
#         "results_merged/{group}/rna/seurat_load_transfer_labels/umap_azimuth_fine.pdf",
#     output:
#         "results_merged/{group}/rna/seurat_load_transfer_labels/reference_plots.pdf"
#     conda:
#         "../envs/fetch.yaml"
#     shell:
#         "pdfunite {input} {output}; "

# rule seurat_cluster_rna:
#     """
#     Seurat RNA clustering
#     """
#     input:
#         project_in = "results_merged/{group}/rna/seurat_load_transfer_labels/proj.rds",
#     output:
#         project_out = "results_merged/{group}/rna/seurat_cluster_rna/proj.rds",
#         umap = "results_merged/{group}/rna/seurat_cluster_rna/umap_clusters.pdf",
#     params:
#         seed = config["seurat_seed"],
#     log:
#         console = "logs/merged/{group}/rna/seurat_cluster_rna/console.log"
#     conda:
#         "../envs/seurat.yaml"
#     script:
#         "../scripts/seurat_cluster_rna.R"

# rule seurat_plot_cm:
#     """
#     Plot confusion matrices for clustering
#     """
#     input:
#         project_in = "results_merged/{group}/rna/seurat_cluster_rna/proj.rds",
#     output:
#         mat_kramann_coarse = "results_merged/{group}/rna/seurat_plot_cm/mat_kramann_coarse.pdf",
#         mat_kramann_fine = "results_merged/{group}/rna/seurat_plot_cm/mat_kramann_fine.pdf",
#         mat_ellinor_coarse = "results_merged/{group}/rna/seurat_plot_cm/mat_ellinor_coarse.pdf",
#         mat_ellinor_fine = "results_merged/{group}/rna/seurat_plot_cm/mat_ellinor_fine.pdf",
#         mat_teichmann_coarse = "results_merged/{group}/rna/seurat_plot_cm/mat_teichmann_coarse.pdf",
#         mat_teichmann_fine = "results_merged/{group}/rna/seurat_plot_cm/mat_teichmann_fine.pdf",
#         mat_azimuth_coarse = "results_merged/{group}/rna/seurat_plot_cm/mat_azimuth_coarse.pdf",
#         mat_azimuth_fine = "results_merged/{group}/rna/seurat_plot_cm/mat_azimuth_fine.pdf",
#     params:
#         seed = config["seurat_seed"],
#     log:
#         console = "logs/merged/{group}/rna/seurat_plot_cm/console.log"
#     conda:
#         "../envs/seurat.yaml"
#     script:
#         "../scripts/seurat_plot_cm.R"

# rule seurat_merge_cm_plots:
#     """
#     Merge confusion matrix plot pdf's
#     """
#     input:
#         "results_merged/{group}/rna/seurat_plot_cm/mat_kramann_coarse.pdf",
#         "results_merged/{group}/rna/seurat_plot_cm/mat_kramann_fine.pdf",
#         "results_merged/{group}/rna/seurat_plot_cm/mat_ellinor_coarse.pdf",
#         "results_merged/{group}/rna/seurat_plot_cm/mat_ellinor_fine.pdf",
#         "results_merged/{group}/rna/seurat_plot_cm/mat_teichmann_coarse.pdf",
#         "results_merged/{group}/rna/seurat_plot_cm/mat_teichmann_fine.pdf",
#         "results_merged/{group}/rna/seurat_plot_cm/mat_azimuth_coarse.pdf",
#         "results_merged/{group}/rna/seurat_plot_cm/mat_azimuth_fine.pdf",
#     output:
#         "results_merged/{group}/rna/seurat_plot_cm/mats.pdf"
#     conda:
#         "../envs/fetch.yaml"
#     shell:
#         "pdfunite {input} {output}; "

# rule seurat_integrate_l2:
#     """
#     Integrate all RNA samples using Harmony
#     """
#     input:
#         projects_in = lambda w: [f"results_merged/{i}/rna/seurat_cluster_rna/proj.rds" for i in group_names]
#     output:
#         project_out = "results_merged/all/rna/seurat_integrate_l2/proj.rds",
#         umap_dataset_pre_harmony = "results_merged/all/rna/seurat_integrate_l2/umap_dataset_pre_harmony.pdf",
#         umap_mixing_pre_harmony = "results_merged/all/rna/seurat_integrate_l2/umap_mixing_pre_harmony.pdf",
#         umap_dataset_harmony = "results_merged/all/rna/seurat_integrate_l2/umap_dataset_harmony.pdf",
#         umap_mixing_harmony = "results_merged/all/rna/seurat_integrate_l2/umap_mixing_harmony.pdf",
#         umap_group_harmony = "results_merged/all/rna/seurat_integrate_l2/umap_group_harmony.pdf",
#         umap_region_harmony = "results_merged/all/rna/seurat_integrate_l2/umap_region_harmony.pdf",
#         umap_status_harmony = "results_merged/all/rna/seurat_integrate_l2/umap_status_harmony.pdf",
#     params:
#         seed = config["seurat_seed"],
#         groups = group_names
#     log:
#         console = "logs/merged/all/rna/seurat_integrate_l2/console.log"
#     conda:
#         "../envs/seurat.yaml"
#     script:
#         "../scripts/seurat_integrate_l2.R"

# rule seurat_merge_integration_plots_l2:
#     """
#     Merge integration QC pdf's
#     """
#     input:
#         "results_merged/all/rna/seurat_integrate_l2/umap_dataset_pre_harmony.pdf",
#         "results_merged/all/rna/seurat_integrate_l2/umap_mixing_pre_harmony.pdf",
#         "results_merged/all/rna/seurat_integrate_l2/umap_dataset_harmony.pdf",
#         "results_merged/all/rna/seurat_integrate_l2/umap_mixing_harmony.pdf",
#         "results_merged/all/rna/seurat_integrate_l2/umap_group_harmony.pdf",
#         "results_merged/all/rna/seurat_integrate_l2/umap_region_harmony.pdf",
#         "results_merged/all/rna/seurat_integrate_l2/umap_status_harmony.pdf",
#     output:
#         "results_merged/all/rna/seurat_integrate_l2/integration_plots.pdf"
#     conda:
#         "../envs/fetch.yaml"
#     shell:
#         "pdfunite {input} {output}; "

# rule seurat_cluster_l2:
#     """
#     Seurat RNA clustering
#     """
#     input:
#         project_in = "results_merged/all/rna/seurat_integrate_l2/proj.rds"
#     output:
#         project_out = "results_merged/all/rna/seurat_cluster_rna/proj.rds",
#         umap = "results_merged/all/rna/seurat_cluster_rna/umap_clusters.pdf",
#         umap_kramann_coarse = "results_merged/all/rna/seurat_cluster_rna/umap_kramann_coarse.pdf",
#         umap_kramann_fine = "results_merged/all/rna/seurat_cluster_rna/umap_kramann_fine.pdf",
#         umap_ellinor_coarse = "results_merged/all/rna/seurat_cluster_rna/umap_ellinor_coarse.pdf",
#         umap_ellinor_fine = "results_merged/all/rna/seurat_cluster_rna/umap_ellinor_fine.pdf",
#         umap_teichmann_coarse = "results_merged/all/rna/seurat_cluster_rna/umap_teichmann_coarse.pdf",
#         umap_teichmann_fine = "results_merged/all/rna/seurat_cluster_rna/umap_teichmann_fine.pdf",
#         umap_azimuth_coarse = "results_merged/all/rna/seurat_cluster_rna/umap_azimuth_coarse.pdf",
#         umap_azimuth_fine = "results_merged/all/rna/seurat_cluster_rna/umap_azimuth_fine.pdf",
#     params:
#         seed = config["seurat_seed"],
#     log:
#         console = "logs/merged/all/rna/seurat_cluster_rna/console.log"
#     conda:
#         "../envs/seurat.yaml"
#     script:
#         "../scripts/seurat_cluster_l2.R"

# rule seurat_merge_label_plots_l2:
#     """
#     Merge reference projection plot pdf's
#     """
#     input:
#         "results_merged/all/rna/seurat_cluster_rna/umap_clusters.pdf",
#         "results_merged/all/rna/seurat_cluster_rna/umap_kramann_coarse.pdf",
#         "results_merged/all/rna/seurat_cluster_rna/umap_kramann_fine.pdf",
#         "results_merged/all/rna/seurat_cluster_rna/umap_ellinor_coarse.pdf",
#         "results_merged/all/rna/seurat_cluster_rna/umap_ellinor_fine.pdf",
#         "results_merged/all/rna/seurat_cluster_rna/umap_teichmann_coarse.pdf",
#         "results_merged/all/rna/seurat_cluster_rna/umap_teichmann_fine.pdf",
#         "results_merged/all/rna/seurat_cluster_rna/umap_azimuth_coarse.pdf",
#         "results_merged/all/rna/seurat_cluster_rna/umap_azimuth_fine.pdf",
#     output:
#         "results_merged/all/rna/seurat_cluster_rna/label_plots.pdf"
#     conda:
#         "../envs/fetch.yaml"
#     shell:
#         "pdfunite {input} {output}; "

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
        umap_samples = "results_merged/rna/seurat_integrate_rna/umap_dataset.pdf"
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