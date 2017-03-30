from indra.literature import pubmed_client
import pickle
import csv

if __name__ == '__main__':
    # Open the list of Reactome signaling genes
    with open('combined_genes.txt', 'rt') as f:
        genes = [line.strip() for line in f.readlines()]

    # Assemble a list of PMIDs curated in Entrez gene
    pmids_for_genes = {}
    for gene_ix, gene in enumerate(genes):
        print("%d of %d: Getting PMIDs for %s" % (gene_ix + 1, len(genes),
                                                  gene))
        try:
            pmids_for_genes[gene] = pubmed_client.get_ids_for_gene(gene)
        except:
            print("Skipping %s" % gene)
            continue
    # Save the dict mapping genes to publications
    with open('combined_genes_to_pmids.pkl', 'wb') as f:
        pickle.dump(pmids_for_genes, f, protocol=2)

    # Get the list of unique PMIDs
    pmids = set([pmid for pmid_list in pmids_for_genes.values()
                      for pmid in pmid_list])

    # Save the PMIDs to a file
    print("Saving PMIDs")
    with open('combined_pmids.txt', 'wt') as f:
        for pmid in pmids:
            f.write('%s\n' % pmid)
