
def main(in_path, out_path):
    with open(in_path) as fi, open(out_path, "w") as fo:
        next(fi)
        for line in fi:
            archr_bc = line.rstrip("\n").split("\t")[0]
            bc = archr_bc.split("#")[1]
            fo.write(f"{bc}\n")

in_path, = snakemake.input
out_path, = snakemake.output

main(in_path, out_path)