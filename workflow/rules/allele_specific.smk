rule build_rasqual:
    output:
        build_dir = directory("resources/rasqual")
    params:
        rasqual_url = config["rasqual_url"]
    conda:
        "../envs/allele_specific.yaml"
    shell:
        "mkdir -p {output.build_dir}; "
        "git clone {params.rasqual_url} {output.build_dir}; "
        "cd {output.build_dir}/src; "
        "export CFLAGS=\"-I$CONDA_PREFIX/include/INCLUDE -I$CONDA_PREFIX/include -fcommon\"; "
        "export LDFLAGS=\"-L$CONDA_PREFIX/lib\"; "
        "make; "
        "make install; "
        "cd ASVCF; "
        "make"


rule merge_vcfs:
    input:
        lambda w: [config["vcf_path"].format(dir=d) for d, _ in groups_unisex_wgs[w.group]]
    output:
        "results_unisex/{group}/allele_specific/genotypes/merged.vcf.gz"
    conda:
        "../envs/allele_specific.yaml"
    threads:
        8
    shell:
        "bcftools merge --force-samples -0 -O z -o {output} --threads {threads} {input}"


rule index_vcf:
    input:
        "results_unisex/{group}/allele_specific/genotypes/merged.vcf.gz"
    output:
        "results_unisex/{group}/allele_specific/genotypes/merged.vcf.gz.tbi"
    conda:
        "../envs/allele_specific.yaml"
    shell:
        "bcftools index --tbi {input}"


rule split_bam_atac:
    input:
        metadata = lambda w: f"export_l1/{sample_to_group[w.sample]}/metadata.tsv.gz",
        labels = lambda w: f"export_l1/{sample_to_group[w.sample]}/labels/cell_types_l1.tsv.gz",
        dataset_info = "results/{sample}/fetch/dataset_info.json",
        bam = "results/{sample}/fetch/atac.bam"
    output:
        bam = "results/{sample}/split_bam/label/{label}/atac.bam"
    conda:
        "../envs/allele_specific.yaml"
    script:
        "../scripts/allele_specific/split_bam_atac.py"


rule index_bam_atac:
    input:
        "results/{sample}/split_bam/label/{label}/atac.bam"
    output:
        "results/{sample}/split_bam/label/{label}/atac.bam.bai"
    conda:
        "../envs/allele_specific.yaml"
    shell:
        "samtools index {input}"


rule rasqual_bam_list_atac:
    input:
        bams = lambda w: [f"results/{s}/split_bam/label/{w.label}/atac.bam" for _, s in groups_unisex_wgs[w.group]],
        inds = lambda w: [f"results/{s}/split_bam/label/{w.label}/atac.bam.bai" for _, s in groups_unisex_wgs[w.group]]
    output:
        "results_unisex/{group}/allele_specific/label/{label}/atac/bam_list.txt"
    conda:
        "../envs/allele_specific.yaml"
    shell:
        "echo {input.bams} | tr ' ' '\\n' > {output}"


rule vcf_as_counts_atac:
    input:
        vcf = "results_unisex/{group}/allele_specific/genotypes/merged.vcf.gz",
        index = "results_unisex/{group}/allele_specific/genotypes/merged.vcf.gz.tbi",
        bam_list = "results_unisex/{group}/allele_specific/label/{label}/atac/bam_list.txt",
        rasqual_dir = "resources/rasqual"
    output:
        "results_unisex/{group}/allele_specific/label/{label}/atac/as_counts.vcf.gz",
    conda:
        "../envs/allele_specific.yaml"
    shell:
        "export RASQUALDIR={input.rasqual_dir}; "
        "bash {input.rasqual_dir}/src/ASVCF/createASVCF.sh paired_end {input.bam_list} {input.vcf} {output} atac"


rule index_vcf_as_atac:
    input:
        "results_unisex/{group}/allele_specific/label/{label}/atac/as_counts.vcf.gz",
    output:
        "results_unisex/{group}/allele_specific/label/{label}/atac/as_counts.vcf.gz.tbi",
    conda:
        "../envs/allele_specific.yaml"
    shell:
        "bcftools index --tbi {input}"


rule organize_total_counts_atac:
    input:
        counts = config["total_counts_atac_path"]
    output:
        counts = "results_unisex/{group}/allele_specific/label/{label}/atac/total_counts.tsv",
    params:
        samples = lambda w: [s for _, s in groups_unisex_wgs[w.group]]
    conda:
        "../envs/allele_specific.yaml"
    script:
        "../scripts/allele_specific/organize_total_counts.py"


rule rasqual_offsets_atac:
    input:
        ytxt = "results_unisex/{group}/allele_specific/label/{label}/atac/total_counts.tsv"
    output:
        ktxt = "results_unisex/{group}/allele_specific/label/{label}/atac/offsets.txt",
    conda:
        "../envs/allele_specific.yaml"
    script:
        "../scripts/allele_specific/rasqual_makeOffset.R"


rule rasqual_txt_to_bin_atac:
    input:
        ytxt = "results_unisex/{group}/allele_specific/label/{label}/atac/total_counts.tsv",
        ktxt = "results_unisex/{group}/allele_specific/label/{label}/atac/offsets.txt",
    output:
        ybin = "results_unisex/{group}/allele_specific/label/{label}/atac/total_counts.bin",
        kbin = "results_unisex/{group}/allele_specific/label/{label}/atac/offsets.bin",
    conda:
        "../envs/allele_specific.yaml"
    script:
        "../scripts/allele_specific/rasqual_txt2bin.R"


rule run_rasqual_atac:
    input:
        rasqual_dir = "resources/rasqual",
        vcf = "results_unisex/{group}/allele_specific/label/{label}/atac/as_counts.vcf.gz",
        index = "results_unisex/{group}/allele_specific/label/{label}/atac/as_counts.vcf.gz.tbi",
        total_counts = "results_unisex/{group}/allele_specific/label/{label}/atac/total_counts.bin",
        offsets = "results_unisex/{group}/allele_specific/label/{label}/atac/offsets.bin",
        peaks = get_peaks,
    output:
        "results_unisex/{group}/allele_specific/label/{label}/atac/rasqual_output.txt",
    params:
        samples = lambda w: [s for _, s in groups_unisex_wgs[w.group]]
    conda:
        "../envs/allele_specific.yaml"
    script:
        "../scripts/allele_specific/run_rasqual.py"