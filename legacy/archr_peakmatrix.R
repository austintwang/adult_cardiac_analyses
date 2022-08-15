Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(parallel)

library(ArchR)
library(rtracklayer)

addArchRVerbose(verbose = FALSE)
addArchRChrPrefix(chrPrefix = FALSE)

# Disable HDF5 file locking
# Workaround for HDF5 I/O issues on NFS
# https://github.com/GreenleafLab/ArchR/issues/248#issuecomment-789453997
Sys.setenv("HDF5_USE_FILE_LOCKING" = "FALSE")
Sys.setenv("RHDF5_USE_FILE_LOCKING" = "FALSE")

##########

params = snakemake@params 
input_paths = snakemake@input
output_paths = snakemake@output
threads = snakemake@threads
log_paths = snakemake@log

seed <- params[["seed"]]
set.seed(seed)

addArchRThreads(threads = threads)

##########

# Load and move project
proj <- loadArchRProject(path = input_paths[["project_in"]], force = FALSE, showLogo = FALSE)
proj <- saveArchRProject(
    ArchRProj = proj,
    outputDirectory = output_paths[["project_out"]],
    overwrite = TRUE,
    load = TRUE,
    logFile = log_paths[["move"]],
)

##########

proj <- addGroupCoverages(
  ArchRProj = proj,
  groupBy = "Clusters_ATAC",
  logFile = log_paths[["pseudobulks"]]
)

# Call peaks
proj <- addReproduciblePeakSet(
    ArchRProj = proj, 
    groupBy = "Clusters_ATAC", 
    logFile = log_paths[["call_peaks"]]
)
proj <- addPeakMatrix(proj, logFile = log_paths[["add_peak_mat"]])

##########

se <- getMatrixFromProject(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  logFile = log_paths[["save_peak_mat"]]
)
print(rowRanges(se)) ####
mat <- Matrix(assay(se), sparse = TRUE)
barcodes <- colData(se)
peaks <- rowRanges(se)

writeMM(obj = mat, file=output_paths[["mat_out"]])
write.table(barcodes, file=output_paths[["barcodes_out"]], quote=FALSE, sep='\t', col.names = NA)
export.bed(peaks, con=output_paths[["peaks_out"]])


saveArchRProject(
    ArchRProj = proj,
    outputDirectory = output_paths[["project_out"]],
    overwrite = TRUE,
    load = FALSE,
    logFile = log_paths[["save"]],
)


sink(type = "message")
sink()