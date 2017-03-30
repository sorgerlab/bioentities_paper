import pandas
from indra.databases import uniprot_client

def get_signaling_genes():
    df = pandas.read_csv('signaling_proteins.tsv', sep='\t', index_col=None)
    gene_names = []
    for up_id in df['Identifier']:
        #terms = mn.split(' ')
        #if len(terms) == 2:
        #    gene_name = terms[1]
        #    gene_names.append(gene_name)
        if not uniprot_client.is_human(up_id):
            print("%s is not a human gene" % up_id)
            continue
        gene_name = uniprot_client.get_gene_name(up_id)
        if not gene_name:
            print("Could not get gene name for %s" % up_id)
            continue
        gene_names.append(gene_name)
    gene_names = sorted(list(set(gene_names)))
    with open('signaling_genes.txt', 'wt') as fh:
        for gn in gene_names:
            fh.write('%s\n' % gn)
    return gene_names

def get_relations_genes():
    df = pandas.read_csv('../relations.csv', index_col=None, header=None)
    gene_names = []
    for _, row in df.iterrows():
        if row[0] == 'HGNC':
            gene_names.append(row[1])
    gene_names = sorted(list(set(gene_names)))
    return gene_names

if __name__ == '__main__':
    signaling_genes = get_signaling_genes()
    relations_genes = get_relations_genes()
    covered_genes = set(signaling_genes).intersection(set(relations_genes))
    missing_genes = set(signaling_genes).difference(set(relations_genes))
    combined_genes = set(signaling_genes).union(set(relations_genes))

    with open('combined_genes.txt', 'wt') as f:
        for gene in combined_genes:
            f.write('%s\n' % gene)

