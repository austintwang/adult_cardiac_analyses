import os
import gzip

REV_COMP = str.maketrans("ATGC", "TACG")
def reverse_complement(seq):
    return str.translate(seq, REV_COMP)[::-1]

def load_barcode_map(bc_atac_path, bc_rna_path):
    bc_map = {}
    with open(bc_atac_path) as a, open(bc_rna_path) as r:
        for la, lr in zip(a, r):
            ba = la.rstrip("\n")
            br = lr.rstrip("\n")
            bc_map[ba] = br
            ba_rc = reverse_complement(ba)
            bc_map[ba_rc] = br

    return bc_map

def load_barcodes(barcodes_path):
    barcodes = []
    with open(barcodes_path) as f:
        h = f.readline().rstrip('\n').split(',')
        bc_ind = h.index("barcode")
        for line in f:
            entries = line.rstrip("\n").split(",")
            b = entries[bc_ind]
            barcodes.append(b)

    return barcodes

def load_records(barcodes, metadata_path, barcode_map):
    records = {}
    with open(metadata_path) as f:
        h = f.readline().rstrip('\n').split('\t')
        cell_id_ind = h.index("cellNames")
        frag_ind = h.index("nFrags")
        tss_ind = h.index("TSSEnrichment")

        for line in f:
            entries = line.rstrip('\n').split('\t')

            cell_id = entries[cell_id_ind]
            _, barcode = cell_id.split("#")
            barcode_rna = barcode_map[barcode]
            
            frag_count = int(entries[frag_ind])
            tss_enr = float(entries[tss_ind])

            record = {
                "barcode_rna": barcode_rna,
                "barcode_atac": barcode,
                "frag_count": frag_count,
                "tss_enr": tss_enr,
                "pass_atac_filter": False
            }
            records[barcode] = record

    for b in barcodes:
        if b not in records:
            barcode = b
            barcode_rna = barcode_map[barcode]
            record = {
                "barcode_rna": barcode_rna,
                "barcode_atac": barcode,
                "frag_count": "NA",
                "tss_enr": "NA",
                "pass_atac_filter": False
            }
            records[b] = record

    return records

def load_final_data(final_data_path):
    ids = []
    with open(final_data_path) as f:
        next(f)
        for line in f:
            cell_id = line.rstrip('\n').split('\t')[0]
            _, barcode = cell_id.split("#")
            ids.append(barcode)

    return ids

def load_amulet_data(barcodes, amulet_data_path):
    records = {}
    with open(amulet_data_path) as f:
        h = f.readline().rstrip('\n').split('\t')
        barcode_ind = h.index("barcode")
        p_ind = h.index("p-value")
        q_ind = h.index("q-value")

        for line in f:
            entries = line.rstrip('\n').split('\t')

            barcode = entries[barcode_ind]
            
            pval = float(entries[p_ind])
            qval = float(entries[q_ind])

            record = {
                "amulet_pval": pval,
                "amulet_qval": qval
            }
            records[barcode] = record

    for b in barcodes:
        if b not in records:
            record = {
                "amulet_pval": "NA",
                "amulet_qval": "NA"
            }
            records[b] = record

    return records


def main(sample_name, metadata_dir, final_data_path, amulet_data_dir, bc_atac_path, bc_rna_path, barcodes_path, out_path):
    metadata_path = os.path.join(metadata_dir, sample_name, "metadata.tsv")

    barcodes = load_barcodes(barcodes_path)

    bc_map = load_barcode_map(bc_atac_path, bc_rna_path)
    records = load_records(barcodes, metadata_path, bc_map)

    final_ids = load_final_data(final_data_path)
    for i in final_ids:
        records[i]["pass_filter"] = True

    amulet_data_path = os.path.join(amulet_data_dir, "MultipletProbabilities.txt")
    records_amulet = load_amulet_data(barcodes, amulet_data_path)
    for k, v in records_amulet.items():
        records[k].update(v)

    with gzip.open(out_path, "wt") as f:
        sample = next(iter(records.values()))
        cols = list(sample.keys())
        f.write("\t".join(cols) + "\n")

        for cell_id, r in records.items():
            line = "\t".join(r[c] for c in cols) + "\n"
            f.write(line)

metadata_dir = snakemake.input["metadata"]
final_data_path = snakemake.input["final_data"]
amulet_data_dir = snakemake.input["amulet"]
bc_atac_path = snakemake.input["bc_atac"]
bc_rna_path = snakemake.input["bc_rna"]
barcodes_path = snakemake.input["barcodes"]

out_path, = snakemake.output

sample_name = snakemake.params["sample_name"]

main(sample_name, metadata_dir, final_data_path, amulet_data_dir, bc_atac_path, bc_rna_path, barcodes_path, out_path)