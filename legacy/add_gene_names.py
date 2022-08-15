import mygene

def query_genes(gene_list):
    mg = mygene.MyGeneInfo()
    query = mg.querymany(gene_list , scopes='ensembl.gene', fields='symbol', species='human')
    mapping = {i["query"]: i["symbol"] for i in query if "symbol" in i}
    names = [mapping.get(g, g) for g in gene_list]
    return names

def main(in_path, out_path):
    ids = []
    with open(in_path) as f:
        for line in f:
            gene_id = line.rstrip("\n").split()[0]
            ids.append(gene_id)

    names = query_genes(ids)
    with open(out_path, "w") as f:
        for i, n in zip(ids, names):
            f.write(f"{i}\t{n}\n")

in_path, = snakemake.input
out_path, = snakemake.output

main(in_path, out_path)