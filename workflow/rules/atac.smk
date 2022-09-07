

rule archr_install_bsgenome:
    """
    Install custom bsgenome object for ArchR
    """
    input:
        bsgenome = "resources/bsgenome.tar.gz",
    output:
        bsgenome_library_dir = directory("resources/bsgenome_lib")
    log:
        console = "logs/resources/archr_install_bsgenome/console.log",
    conda:
        "../envs/archr.yaml"
    script:
        "../scripts/archr_install_bsgenome.R"

rule archr_build:
    """
    Build ArchR project
    """
    input:
        frag = "results/{sample}/fetch/fragments.tsv.gz", 
        frag_ind = "results/{sample}/fetch/fragments.tsv.gz.tbi",
        blacklist = "resources/blacklist.bed",
        bsgenome_library_dir = "resources/bsgenome_lib",
        gene_anno = "resources/gene_annotation.rda"
    output:
        qc_dir = directory("results/{sample}/atac/archr_qc"),
        project_dir = directory("results/{sample}/atac/archr_init"),
        arrow_dir = directory("results/{sample}/atac/archr_arrows"),
        metadata = "results/{sample}/atac/archr_metadata.tsv"
    params:
        sample_name = lambda w: w.sample,
        seed = config["archr_seed"],
        min_frags = config["archr_min_frags"],
        min_tss_enr = config["archr_min_tss_enr"],
        bsgenome = config["bsgenome_name"],
        gene_anno = config["gene_anno_name"]
    log:
        console = "logs/{sample}/atac/archr_build/console.log",
        arrow_create = "logs/{sample}/atac/archr_build/arrow_create.log",
        save = "logs/{sample}/atac/archr_build/save.log"
    threads:
        max_threads
    conda:
        "../envs/archr.yaml"
    script:
        "../scripts/archr_build.R"

rule archr_write_qc:
    """
    ArchR write sample QC data
    """
    input:
        qc_dir = "results/{sample}/atac/archr_qc"
    output:
        out_dir = directory("results/{sample}/atac/archr_qc_parsed")
    log:
        console = "logs/{sample}/atac/archr_write_qc/console.log"
    conda:
        "../envs/archr.yaml"
    script:
        "../scripts/archr_write_qc.R"

rule extract_barcodes:
    """
    Extract barcode data from fragments
    """
    input:
        "results/{sample}/fetch/fragments.tsv.gz"
    output:
        "results/{sample}/atac/amulet_barcode_data.csv"
    conda:
        "../envs/fetch.yaml"
    shell:
        "echo 'barcode,is__cell_barcode' > {output} && " 
        "zcat {input} | awk '{{ print $4 }}' - | sort -u | sed $'s/$/,1/' >> {output}"

rule write_amulet_barcodes:
    """
    Filter amulet input barcodes
    """
    input:
        metadata = "results/{sample}/atac/archr_qc_parsed",
        amulet_barcodes = "results/{sample}/atac/amulet_barcode_data.csv"
    output:
        "results/{sample}/atac/amulet_barcode_data_final.csv"
    params:
        sample_name = lambda w: w.sample,
        min_frags = config["archr_min_frags"],
        min_tss_enr = config["archr_min_tss_enr"]
    conda:
        "../envs/fetch.yaml"
    script:
        "../scripts/write_amulet_barcodes.py"
    
rule run_amulet:
    """
    Run AMULET batch correction tool
    """
    input:
        frag = "results/{sample}/fetch/fragments.tsv.gz", 
        barcodes = "results/{sample}/atac/amulet_barcode_data_final.csv",
        chroms = "resources/GRCh38.chroms.tsv",
        blacklist = "resources/blacklist.bed",
        amulet_dir = "resources/amulet",
        amulet_lowmem_dir = "resources/amulet_lowmem"
    output:
        directory("results/{sample}/atac/amulet")
    conda:
        "../envs/amulet.yaml"
    shell:
        "mkdir -p {output}; "
        "{input.amulet_dir}/AMULET.sh {input.frag} {input.barcodes} {input.chroms} {input.blacklist} {output} {input.amulet_dir} || "
        "{input.amulet_lowmem_dir}/AMULET.sh {input.frag} {input.barcodes} {input.chroms} {input.blacklist} {output} {input.amulet_lowmem_dir}"

rule write_atac_qc:
    """
    Write ATAC QC data
    """
    input:
        metadata = "results/{sample}/atac/archr_qc_parsed",
        final_data = "results/{sample}/atac/archr_metadata.tsv",
        amulet = "results/{sample}/atac/amulet",
        bc_atac = "resources/whitelist_atac.txt",
        bc_rna = "resources/whitelist_rna.txt",
        barcodes = "results/{sample}/atac/amulet_barcode_data.csv"
    output:
        "results/{sample}/atac/atac_qc.tsv"
    params:
        sample_name = lambda w: w.sample,
    conda:
        "../envs/fetch.yaml"
    script:
        "../scripts/write_atac_qc.py"