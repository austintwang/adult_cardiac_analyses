
def main(in_path, out_path):
    with open(in_path) as fi, open(out_path, "w") as fo:
        for line in fi:
            chrom, start, end = line.rstrip("\n").split("\t")[:3]
            fo.write(f"{chrom}:{start}-{end}\n")

in_path, = snakemake.input
out_path, = snakemake.output

main(in_path, out_path)