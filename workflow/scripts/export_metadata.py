import gzip
import json

HEADER = """
# cell_id: Cell ID used in integrated analysis
# rna_dataset: ENCODE snRNA-Seq dataset ID
# rna_barcode: snRNA-Seq barcode
# atac_dataset: ENCODE snATAC-Seq dataset ID
# atac_barcode: snATAC-Seq barcode
# rna_umi_count: snRNA UMI's per cell
# atac_fragment_count: snATAC-Seq fragments per cell
# rna_frac_mito: Fraction of mitochondrial RNA reads
# rna_frac_ribo: Fraction of ribosomal RNA reads
# atac_tss_enrichment: snATAC-Seq transcription start site enrichment
# passed_filtering: A binary (0/1) value indicating whether the cell passed manual filtering
"""

COLUMNS = (
    "cell_id\trna_dataset\trna_barcode\tatac_dataset\tatac_barcode\t"
    "rna_umi_count\tatac_fragment_count\trna_frac_mito\trna_frac_ribo\tatac_tss_enrichment\tpassed_filtering\n"
)

def load_records_rna(metadata_path, dataset_data):
    records = {}
    with open(metadata_path) as f:
        h = f.readline().rstrip('\n').split('\t')

        dataset_ind = h.index("dataset")
        barcode_rna_ind = h.index("barcode_rna")
        barcode_atac_ind = h.index("barcode_atac")
        count_ind = h.index("nCount_RNA")
        mito_ind = h.index("percent.mt")
        ribo_ind = h.index("percent.ribo")

        for line in f:
            entries = line.rstrip('\n').split('\t')

            dataset = entries[dataset_ind]
            dataset_rna = dataset_data[dataset]["rna"]
            dataset_atac = dataset_data[dataset]["atac"]
            barcode_rna = entries[barcode_rna_ind]
            barcode_atac = entries[barcode_atac_ind]
            cell_id = f"{dataset}_{barcode_rna}"
            
            umi_count = int(entries[count_ind])
            frac_mito = float(entries[mito_ind]) / 100
            frac_ribo = float(entries[ribo_ind]) / 100

            record = {
                "dataset_rna": dataset_rna,
                "dataset_atac": dataset_atac,
                "barcode_rna": barcode_rna,
                "barcode_atac": barcode_atac,
                "umi_count": umi_count,
                "frac_mito": frac_mito,
                "frac_ribo": frac_ribo,
                "frag_count": "NA",
                "tss_enr": "NA",
                "pass_filter": False
            }
            records[cell_id] = record

    return records

def add_records_atac(metadata_path, dataset_data, records):

    with gzip.open(metadata_path, "rt") as f:
        h = f.readline().rstrip('\n').split('\t')

        dataset_ind = h.index("sample")
        barcode_rna_ind = h.index("barcode_rna")
        barcode_atac_ind = h.index("barcode_atac")
        frag_ind = h.index("frag_count")
        tss_ind = h.index("tss_enr")

        for line in f:
            entries = line.rstrip('\n').split('\t')

            dataset = entries[dataset_ind]
            dataset_rna = dataset_data[dataset]["rna"]
            dataset_atac = dataset_data[dataset]["atac"]
            barcode_rna = entries[barcode_rna_ind]
            barcode_atac = entries[barcode_atac_ind]
            cell_id = f"{dataset}_{barcode_rna}"
            
            frag_count = int(entries[frag_ind])
            tss_enr = float(entries[tss_ind])

            if cell_id in records:
                records[cell_id].update({
                    "frag_count": frag_count,
                    "tss_enr": tss_enr
                })
            else:
                record = {
                    "dataset_rna": dataset_rna,
                    "dataset_atac": dataset_atac,
                    "barcode_rna": barcode_rna,
                    "barcode_atac": barcode_atac,
                    "umi_count": "NA",
                    "frac_mito": "NA",
                    "frac_ribo": "NA",
                    "frag_count": frag_count,
                    "tss_enr": tss_enr,
                    "pass_filter": False
                }
                records[cell_id] = record


def load_final_data(final_data_path):
    ids = []
    with open(final_data_path) as f:
        next(f)
        for line in f:
            cell_id = line.rstrip('\n').split('\t')[0]
            ids.append(cell_id)

    return ids

def main(rna_metadata_paths, atac_metadata_paths, dataset_data_paths, final_data_path, out_path, out_datasets_path):
    dataset_data = {}
    for i in dataset_data_paths:
        with open(i) as f:
            entry = json.load(i)
        dataset_data |= entry

    with open(out_datasets_path, "w") as f:
        for v1 in dataset_data.values():
            for v2 in v1.values():
                f.write(f"{v2}\n")

    records = {}
    for i in rna_metadata_paths:
        records |= load_records_rna(i, dataset_data)
    for i in atac_metadata_paths:
        add_records_atac(i, dataset_data, records)

    final_ids = load_final_data(final_data_path)
    for i in final_ids:
        records[i]["pass_filter"] = True

    with gzip.open(out_path, "wt") as f:
        f.write(HEADER)
        f.write(COLUMNS)

        for cell_id, r in records.items():
            line = f"{cell_id}\t{r['dataset_rna']}\t{r['barcode_rna']}\t{r['dataset_atac']}\t{r['barcode_atac']}\t{r['umi_count']}\{r['frag_count']}\t{r['frac_mito']}\t{r['frac_ribo']}\t{r['tss_enr']}\t{int(r['pass_filter'])}\n"
            f.write(line)

rna_metadata_paths = snakemake.input["metadata_rna"]
atac_metadata_paths = snakemake.input["metadata_atac"]
dataset_data_paths = snakemake.input["dataset_data"]
final_data_path = snakemake.input["final_data"]

out_path = snakemake.output["metadata"]
out_datasets_path = snakemake.output["metadata"]


main(rna_metadata_paths, atac_metadata_paths, dataset_data_paths, final_data_path, out_path, out_datasets_path)