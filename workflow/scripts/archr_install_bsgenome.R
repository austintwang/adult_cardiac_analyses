console_log <- file(snakemake@log[[1]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

params = snakemake@params 
input_paths = snakemake@input
output_paths = snakemake@output
threads = snakemake@threads
log_paths = snakemake@log

dir.create(out_dir)

bsgenome_path <- input_paths[["bsgenome"]]
dir.create(output_paths[["bsgenome_library_dir"]], recursive = TRUE)
install.packages(bsgenome_path, repos = NULL, type = "source", lib = output_paths[["bsgenome_library_dir"]])

sink(type = "message")
sink()