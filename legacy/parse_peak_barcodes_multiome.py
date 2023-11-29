
REV_COMP = str.maketrans("ATGC", "TACG")
def reverse_complement(seq):
    return str.translate(seq, REV_COMP)[::-1]

def load_whitelists(wl_atac_path, wl_rna_path):
    bc_map = {}
    bc_map_rc = {}
    with open(wl_atac_path) as a, open(wl_rna_path) as r:
        for line_a, line_r in zip(a, r):
            bc_a = line_a.rstrip("\n")
            bc_r = line_r.rstrip("\n")
            bc_a_rc = reverse_complement(bc_a)

            bc_map[bc_a] = bc_r
            bc_map_rc[bc_a_rc] = bc_r

    return bc_map, bc_map_rc

def main(in_path, out_path, wl_atac_path, wl_rna_path):
    bc_map, bc_map_rc = load_whitelists(wl_atac_path, wl_rna_path)

    bcs = []
    with open(in_path) as fi:
        next(fi)
        for line in fi:
            archr_bc = line.rstrip("\n").split("\t")[0]
            bc = archr_bc.split("#")[1]
            bcs.append(bc)

    try:
        bcs_out = [bc_map[i] for i in bcs]
    except KeyError:
        bcs_out = [bc_map_rc[i] for i in bcs]
    
    bcs_out.append("")
    with open(out_path, "w") as f:
        f.write("\n".join(bcs_out))

in_path = snakemake.input["archr_data"]
wl_atac_path = snakemake.input["wl_atac"]
wl_rna_path = snakemake.input["wl_rna"]

out_path, = snakemake.output

main(in_path, out_path, wl_atac_path, wl_rna_path)