import os

def load_passing_bcs(metadata_path, min_frags, min_tss_enr):
    passes = set()
    with open(metadata_path) as f:
        h = f.readline().rstrip('\n').split('\t')
        cell_id_ind = h.index("cellNames")
        frag_ind = h.index("nFrags")
        tss_ind = h.index("TSSEnrichment")

        for line in f:
            entries = line.rstrip('\n').split('\t')

            cell_id = entries[cell_id_ind]
            _, barcode = cell_id.split("#")
            
            frag_count = int(entries[frag_ind])
            tss_enr = float(entries[tss_ind])

            if (frag_count > min_frags) and (tss_enr > min_tss_enr):
                passes.add(barcode)

    return passes

def main(metadata_dir, in_path, out_path, min_frags, min_tss_enr, sample_name):
    metadata_path = os.path.join(metadata_dir, sample_name, "metadata.tsv")
    passes = load_passing_bcs(metadata_path, min_frags, min_tss_enr)

    with open(in_path) as fi, open(out_path, "w") as fo:
        header = fi.readline()
        fo.write(header)

        for line in fi:
            barcode, _ = line.rstrip('\n').split(',')
            pass_filter = barcode in passes
            fo.write(f"{barcode},{int(pass_filter)}\n")

metadata_dir = snakemake.input["metadata"]
in_path = snakemake.input["amulet_barcodes"]

out_path, = snakemake.output

sample_name = snakemake.params["sample_name"]
min_frags = snakemake.params["min_frags"]
min_tss_enr = snakemake.params["min_tss_enr"]

main(metadata_dir, in_path, out_path, min_frags, min_tss_enr, sample_name)