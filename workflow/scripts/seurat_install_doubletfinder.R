console_log <- file(snakemake@log[[1]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

params = snakemake@params 
input_paths = snakemake@input
output_paths = snakemake@output
threads = snakemake@threads
log_paths = snakemake@log

dir.create(output_paths[["doubletfinder_library_dir"]], recursive = TRUE)
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', lib = output_paths[["doubletfinder_library_dir"]], force = TRUE)

sink(type = "message")
sink()