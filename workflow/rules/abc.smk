rule download_abc_genes:
    """
    Download ABC genes reference
    """
    output:
        "resources/abc/CollapsedGeneBounds.hg38.bed"
    params:
        url = config["abc_genes_url"]
    conda:
        "../envs/fetch.yaml"
    shell:
        "curl --no-progress-meter -L -f '{params.url}' > {output}"


rule sort_abc_genes:
    input:
        "resources/abc/CollapsedGeneBounds.hg38.bed"
    output:
        "resources/abc/CollapsedGeneBounds_sorted.hg38.bed"
    conda:
        "../envs/abc.yaml"
    shell:
        "LC_COLLATE='C' sort -k1,1 -k2,2n {input} > {output}"


rule sort_chrom_sizes:
    input:
        "resources/GRCh38.chrom.sizes.tsv"
    output:
        "resources/abc/GRCh38.chrom.sizes.sorted.tsv"
    conda:
        "../envs/abc.yaml"
    shell:
        "LC_COLLATE='C' sort -k1,1 -k2,2n {input} > {output}"


rule generate_chrom_sizes_bed_file:
    input:
        chrom_sizes = "resources/GRCh38.chrom.sizes.sorted.tsv"
    output:
        chrom_sizes_bed = "resources/abc/GRCh38.chrom.sizes.sorted.tsv.bed"

    shell:
        "awk 'BEGIN {{OFS=\"\\t\"}} {{if (NF > 0) print $1,\"0\",$2 ; else print $0}}' {input.chrom_sizes} > {output.chrom_sizes_bed}"


rule download_hic:
    """
    Download HiC data
    """
    output:
        "resources/abc/hic/{experiment}.hic"
    params:
        url = lambda w: config["hic_url"].format(file=hic_files[w.experiment])
    conda:
        "../envs/fetch.yaml"
    shell:
        "curl --no-progress-meter -L -f '{params.url}' > {output}"


status_map = {
    "healthy": "H",
    "ischemic": "I",
    "nonischemic": "NI",
    "inflammatory": "In",
}


def get_pseudobulk_frags(w):
    label = w.label
    region, status = w.group.split("-")
    region = region.upper()
    status = status_map[status]
    return os.path.join(config["pseudobulk_frags_dir"], f"{label}_{region}_{status}_sorted.tsv")

rule frags_to_ta:
    input:
        get_pseudobulk_frags
    output:
        "results_unisex/{group}/abc/{label}/fragments.tagAlign.gz"
    conda:
        "../envs/abc.yaml"
    script:
        "../scripts/abc/fragments_to_ta.py"


rule bgzip_ta:
    input:
        "results_unisex/{group}/abc/{label}/fragments.tagAlign.gz"
    output:
        "results_unisex/{group}/abc/{label}/fragments_bgzip.tagAlign.gz"
    params:
        tmpdir = config["scratch_dir"]
    conda:
        "../envs/abc.yaml"
    shell:
        "zcat {input} "
        "| LC_COLLATE='C' sort -T {params.tmpdir} -k1,1 -k2,2n"
        "| bgzip > {output}"


def get_peaks(w):
    label = w.label
    region, status = w.group.split("-")
    region = region.upper()
    status = status_map[status]
    return os.path.join(config["peaks_dir"], f"{label}_{region}_{status}_peaks_overlap_filtered.narrowPeak")


rule abc_neighborhoods:
    input:
        tagalign = "results_unisex/{group}/abc/{label}/fragments_bgzip.tagAlign.gz",
        peaks = get_peaks,
        chrom_sizes = "resources/abc/GRCh38.chrom.sizes.sorted.tsv",
        chrom_sizes_bed = "resources/abc/GRCh38.chrom.sizes.tsv.bed",
        genes = "resources/abc/CollapsedGeneBounds_sorted.hg38.bed",
    output:
        genes = "results_unisex/{group}/abc/{label}/neighborhoods/GeneList.txt",
        enhancers = "results_unisex/{group}/abc/{label}/neighborhoods/EnhancerList.txt",
    conda:
        "../envs/abc.yaml"
    script:
        "../scripts/abc/neighborhoods.py"


rule abc_predict:
    input:
        genes = "results_unisex/{group}/abc/{label}/neighborhoods/GeneList.txt",
        enhancers = "results_unisex/{group}/abc/{label}/neighborhoods/EnhancerList.txt",
        chrom_sizes = "resources/GRCh38.chrom.sizes.tsv",
        hic = "resources/abc/hic/{experiment}.hic",
    output:
        predictions = "results_unisex/{group}/abc/{label}/predictions/{experiment}/EnhancerPredictionsAllPutative.tsv.gz",
        variants = "results_unisex/{group}/abc/{label}/predictions/{experiment}/EnhancerPredictionsAllPutative.ForVariantOverlap.shrunk150bp.tsv.gz",
    conda:
        "../envs/abc.yaml"
    script:
        "../scripts/abc/predict.py"
