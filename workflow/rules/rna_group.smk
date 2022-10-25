

rule seurat_integrate_rna:
    """
    Integrate RNA samples using Harmony
    """
    input:
        projects_in = lambda w: [f"results/{sample}/rna/seurat_countsplit/proj.rds" for sample in groups[w.group]]
    output:
        project_out = "results_groups/{group}/rna/seurat_integrate_rna/proj.rds",
        umap_dataset_pre_harmony = "results_groups/{group}/rna/seurat_integrate_rna/umap_dataset_pre_harmony.pdf",
        umap_mixing_pre_harmony = "results_groups/{group}/rna/seurat_integrate_rna/umap_mixing_pre_harmony.pdf",
        umap_dataset_harmony = "results_groups/{group}/rna/seurat_integrate_rna/umap_dataset_harmony.pdf",
        umap_mixing_harmony = "results_groups/{group}/rna/seurat_integrate_rna/umap_mixing_harmony.pdf",
    params:
        seed = config["seurat_seed"],
        samples = lambda w: groups[w.group],
    log:
        console = "logs/merged/{group}/rna/seurat_integrate_rna/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_integrate_rna.R"

rule seurat_write_integration_metadata:
    """
    Write integration metadata
    """
    input:
        project_in = "results_groups/{group}/rna/seurat_integrate_rna/proj.rds"
    output:
        metadata = "results_groups/{group}/rna/seurat_integrate_rna/metadata.tsv"
    params:
        seed = config["seurat_seed"]
    log:
        console = "logs/merged/{group}/rna/seurat_integrate_rna/console_metadata.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_write_metadata.R"

rule seurat_integrate_qc_extras:
    """
    Additional integration QC plots
    """
    input:
        project_in = "results_groups/{group}/rna/seurat_integrate_rna/proj.rds"
    output:
        umap_counts = "results_groups/{group}/rna/seurat_integrate_rna/umap_counts.pdf",
        umap_doubletfinder = "results_groups/{group}/rna/seurat_integrate_rna/umap_doubletfinder.pdf",
        umap_amulet = "results_groups/{group}/rna/seurat_integrate_rna/umap_amulet.pdf"
    params:
        seed = config["seurat_seed"],
        samples = lambda w: groups[w.group],
    log:
        console = "logs/merged/{group}/rna/seurat_integrate_rna_extras/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_integration_qc.R"

rule seurat_merge_integration_plots:
    """
    Merge integration QC pdf's
    """
    input:
        "results_groups/{group}/rna/seurat_integrate_rna/umap_dataset_pre_harmony.pdf",
        "results_groups/{group}/rna/seurat_integrate_rna/umap_dataset_harmony.pdf",
        "results_groups/{group}/rna/seurat_integrate_rna/umap_mixing_pre_harmony.pdf",
        "results_groups/{group}/rna/seurat_integrate_rna/umap_mixing_harmony.pdf",
        "results_groups/{group}/rna/seurat_integrate_rna/umap_counts.pdf",
        "results_groups/{group}/rna/seurat_integrate_rna/umap_doubletfinder.pdf",
        "results_groups/{group}/rna/seurat_integrate_rna/umap_amulet.pdf"
    output:
        "results_groups/{group}/rna/seurat_integrate_rna/integration_plots.pdf"
    conda:
        "../envs/fetch.yaml"
    shell:
        "pdfunite {input} {output}; "

rule seurat_load_transfer_labels:
    """
    Add labels to integrated analyses
    """
    input:
        project_integrated = "results_groups/{group}/rna/seurat_integrate_rna/proj.rds",
        tables = expand(
            "results_groups/{group}/rna/seurat_integrate_transfers/{reference}.tsv", 
            reference = ["kramann", "ellinor", "teichmann", "azimuth"], 
            allow_missing=True
        )
    output:
        project_out = "results_groups/{group}/rna/seurat_load_transfer_labels/proj.rds",
        umap_kramann_coarse = "results_groups/{group}/rna/seurat_load_transfer_labels/umap_kramann_coarse.pdf",
        umap_kramann_fine = "results_groups/{group}/rna/seurat_load_transfer_labels/umap_kramann_fine.pdf",
        umap_ellinor_coarse = "results_groups/{group}/rna/seurat_load_transfer_labels/umap_ellinor_coarse.pdf",
        umap_ellinor_fine = "results_groups/{group}/rna/seurat_load_transfer_labels/umap_ellinor_fine.pdf",
        umap_teichmann_coarse = "results_groups/{group}/rna/seurat_load_transfer_labels/umap_teichmann_coarse.pdf",
        umap_teichmann_fine = "results_groups/{group}/rna/seurat_load_transfer_labels/umap_teichmann_fine.pdf",
        umap_azimuth_coarse = "results_groups/{group}/rna/seurat_load_transfer_labels/umap_azimuth_coarse.pdf",
        umap_azimuth_fine = "results_groups/{group}/rna/seurat_load_transfer_labels/umap_azimuth_fine.pdf",
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/merged/{group}/rna/seurat_load_transfer_labels/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_load_transfer_labels.R"

rule seurat_merge_reference_label_plots:
    """
    Merge reference projection plot pdf's
    """
    input:
        "results_groups/{group}/rna/seurat_load_transfer_labels/umap_kramann_coarse.pdf",
        "results_groups/{group}/rna/seurat_load_transfer_labels/umap_kramann_fine.pdf",
        "results_groups/{group}/rna/seurat_load_transfer_labels/umap_ellinor_coarse.pdf",
        "results_groups/{group}/rna/seurat_load_transfer_labels/umap_ellinor_fine.pdf",
        "results_groups/{group}/rna/seurat_load_transfer_labels/umap_teichmann_coarse.pdf",
        "results_groups/{group}/rna/seurat_load_transfer_labels/umap_teichmann_fine.pdf",
        "results_groups/{group}/rna/seurat_load_transfer_labels/umap_azimuth_coarse.pdf",
        "results_groups/{group}/rna/seurat_load_transfer_labels/umap_azimuth_fine.pdf",
    output:
        "results_groups/{group}/rna/seurat_load_transfer_labels/reference_plots.pdf"
    conda:
        "../envs/fetch.yaml"
    shell:
        "pdfunite {input} {output}; "

rule seurat_cluster_rna:
    """
    Seurat RNA clustering
    """
    input:
        project_in = "results_groups/{group}/rna/seurat_load_transfer_labels/proj.rds",
    output:
        project_out = "results_groups/{group}/rna/seurat_cluster_rna/proj.rds",
        umap = "results_groups/{group}/rna/seurat_cluster_rna/umap_clusters.pdf",
    params:
        seed = config["seurat_seed"],
        resolution = config["rna_cluster_resolution"]
    log:
        console = "logs/merged/{group}/rna/seurat_cluster_rna/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_cluster_rna.R"

rule seurat_write_cluster_metadata:
    """
    Seurat RNA clustering
    """
    input:
        project_in = "results_groups/{group}/rna/seurat_cluster_rna/proj.rds"
    output:
        metadata = "results_groups/{group}/rna/seurat_cluster_rna/metadata.tsv"
    params:
        seed = config["seurat_seed"]
    log:
        console = "logs/merged/{group}/rna/seurat_cluster_rna/console_metadata.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_write_metadata.R"

rule seurat_plot_cm:
    """
    Plot confusion matrices for clustering
    """
    input:
        project_in = "results_groups/{group}/rna/seurat_cluster_rna/proj.rds",
    output:
        mat_kramann_coarse = "results_groups/{group}/rna/seurat_plot_cm/mat_kramann_coarse.pdf",
        mat_kramann_fine = "results_groups/{group}/rna/seurat_plot_cm/mat_kramann_fine.pdf",
        mat_ellinor_coarse = "results_groups/{group}/rna/seurat_plot_cm/mat_ellinor_coarse.pdf",
        mat_ellinor_fine = "results_groups/{group}/rna/seurat_plot_cm/mat_ellinor_fine.pdf",
        mat_teichmann_coarse = "results_groups/{group}/rna/seurat_plot_cm/mat_teichmann_coarse.pdf",
        mat_teichmann_fine = "results_groups/{group}/rna/seurat_plot_cm/mat_teichmann_fine.pdf",
        mat_azimuth_coarse = "results_groups/{group}/rna/seurat_plot_cm/mat_azimuth_coarse.pdf",
        mat_azimuth_fine = "results_groups/{group}/rna/seurat_plot_cm/mat_azimuth_fine.pdf",
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/merged/{group}/rna/seurat_plot_cm/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_plot_cm.R"

rule seurat_merge_cm_plots:
    """
    Merge confusion matrix plot pdf's
    """
    input:
        "results_groups/{group}/rna/seurat_plot_cm/mat_kramann_coarse.pdf",
        "results_groups/{group}/rna/seurat_plot_cm/mat_kramann_fine.pdf",
        "results_groups/{group}/rna/seurat_plot_cm/mat_ellinor_coarse.pdf",
        "results_groups/{group}/rna/seurat_plot_cm/mat_ellinor_fine.pdf",
        "results_groups/{group}/rna/seurat_plot_cm/mat_teichmann_coarse.pdf",
        "results_groups/{group}/rna/seurat_plot_cm/mat_teichmann_fine.pdf",
        "results_groups/{group}/rna/seurat_plot_cm/mat_azimuth_coarse.pdf",
        "results_groups/{group}/rna/seurat_plot_cm/mat_azimuth_fine.pdf",
    output:
        "results_groups/{group}/rna/seurat_plot_cm/mats.pdf"
    conda:
        "../envs/fetch.yaml"
    shell:
        "pdfunite {input} {output}; "

rule seurat_add_batch_data:
    """
    Add experimental batch data
    """
    input:
        project_in = "results_groups/{group}/rna/seurat_cluster_rna/proj.rds",
    output:
        project_out = "results_groups/{group}/rna/seurat_add_batch_data/proj.rds",
        umap = "results_groups/{group}/rna/seurat_add_batch_data/umap_batch.pdf",
        metadata = "results_groups/{group}/rna/seurat_add_batch_data/metadata.tsv"
    params:
        seed = config["seurat_seed"],
        samples = lambda w: groups[w.group],
        batches = lambda w: [batch_data[i] for i in groups[w.group]]
    log:
        console = "logs/merged/{group}/rna/seurat_add_batch_data/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_add_batch_data.R"

rule count_bc_multiplexing:
    """
    Calculate barcode multiplexing stats
    """
    input:
        data = expand("results_groups/{group}/rna/seurat_add_batch_data/metadata.tsv", group=groups.keys())
    output:
        data = "qc_global/count_bc_multiplexing/counts.tsv"
    log:
        console = "logs/qc_global/count_bc_multiplexing/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/count_bc_multiplexing.R"

rule seurat_add_multiplexing_data:
    """
    Add barcode multiplexing data
    """
    input:
        project_in = "results_groups/{group}/rna/seurat_add_batch_data/proj.rds",
        data = "qc_global/count_bc_multiplexing/counts.tsv"
    output:
        project_out = "results_groups/{group}/rna/seurat_add_multiplexing_data/proj.rds",
        umap = "results_groups/{group}/rna/seurat_add_multiplexing_data/umap.pdf",
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/merged/{group}/rna/seurat_add_multiplexing_data/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_add_multiplexing_data.R"