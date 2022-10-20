Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(parallel)

library(ArchR)

addArchRVerbose(verbose = FALSE)
addArchRChrPrefix(chrPrefix = FALSE)

# Disable HDF5 file locking
# Workaround for HDF5 I/O issues on NFS
# https://github.com/GreenleafLab/ArchR/issues/248#issuecomment-789453997
Sys.setenv("HDF5_USE_FILE_LOCKING" = "FALSE")
Sys.setenv("RHDF5_USE_FILE_LOCKING" = "FALSE")

build_archr_project <- function(params, input_paths, output_paths, threads, log_paths) {
    arrow_sample_names <- params[["sample_names"]]
    bsgenome_name <- params[["bsgenome"]]
    gene_anno_name <- params[["gene_anno"]]
    seed <- params[["seed"]]
    min_frags <- params[["min_frags"]]
    min_tss_enr <- params[["min_tss_enr"]]

    set.seed(seed)

    seurat_data <- read.table(file = input_paths[["seurat_data"]], sep = '\t', header = TRUE)
    bc_names <- paste0(seurat_data[["dataset"]], "#", seurat_data[["barcode_atac"]])

    blacklist_path <- input_paths[["blacklist"]]
    regions <- read.table(blacklist_path, sep = '\t', header = FALSE)
    colnames(regions) <- c('chr','start','end')
    blacklist <- GRanges(regions)

    # bsgenome_path <- input_paths[["bsgenome"]]
    # install.packages(bsgenome_path, repos = NULL, type = "source")
    library(bsgenome_name, character.only = TRUE, lib.loc=input_paths[["bsgenome_library_dir"]])
    bsgenome <- get(bsgenome_name)

    chromSizes <- GRanges(names(seqlengths(bsgenome)), IRanges(1, seqlengths(bsgenome)))
    chromSizes <- filterChrGR(chromSizes, remove = c("chrM"))
    seqlengths(chromSizes) <- end(chromSizes)

    genome_annotation <- createGenomeAnnotation(
        genome = bsgenome,
        chromSizes = chromSizes,
        blacklist = blacklist
    )

    gene_anno_path <- input_paths[["gene_anno"]]
    load(gene_anno_path)
    gene_annotation <- get(gene_anno_name)
    
    addArchRThreads(threads = threads)

    # addArchRGenome("hg38")

    frag_paths <- input_paths[["frags"]]

    arrow_output_dir <- output_paths[["arrow_dir"]]
    arrow_output_names <- file.path(arrow_output_dir, unlist(arrow_sample_names))
    # print(arrow_output_dir) ####
    dir.create(arrow_output_dir, recursive = TRUE)
    arrows <- createArrowFiles(
        inputFiles = unlist(frag_paths),
        sampleNames = unlist(arrow_sample_names),
        outputNames = arrow_output_names,
        geneAnnotation = gene_annotation,
        genomeAnnotation = genome_annotation,
        offsetPlus = 0,
        offsetMinus = 0,
        minFrags = min_frags,
        minTSS = min_tss_enr,
        addTileMat = TRUE,
        addGeneScoreMat = TRUE,
        force = TRUE,
        subThreading = FALSE, # required for no file locking
        logFile = log_paths[["arrow_create"]],
        QCDir = output_paths[["qc_dir"]]
    )

    # Create project
    dir.create(output_paths[["project_dir"]])
    proj <- ArchRProject(
        ArrowFiles = arrows, 
        outputDirectory = output_paths[["project_dir"]],
        copyArrows = FALSE,
        geneAnnotation = gene_annotation,
        genomeAnnotation = genome_annotation,
    )
    proj <- subsetCells(ArchRProj = proj, cellNames = bc_names)

    saveArchRProject(
        ArchRProj = proj,
        overwrite = TRUE,
        load = FALSE,
        logFile = log_paths[["save"]],
    )

}

build_archr_project(snakemake@params, snakemake@input, snakemake@output, snakemake@threads, snakemake@log)

sink(type = "message")
sink()