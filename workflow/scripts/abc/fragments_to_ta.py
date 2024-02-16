import gzip

def main(fragment_file, ta_file):
    with open(fragment_file, 'rt') as f, gzip.open(ta_file, 'wt') as g:
        for line in f:
            chrom, start, end, bc, _ = line.rstrip().split('\t', maxsplit=4)
            start = int(start)
            end = int(end)
            
            g.write(f'{chrom}\t{start}\t{start+1}\t{bc}\t1000\t+\n')
            g.write(f'{chrom}\t{end-1}\t{end}\t{bc}\t1000\t-\n')


frag, = snakemake.input
ta, = snakemake.output

main(frag, ta)