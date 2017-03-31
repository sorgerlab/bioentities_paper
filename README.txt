- Divide the no-bioentities corpus in to a 90% (training) set and a 10%
  validation set.
- Process the 90% set, and look at the frequency of ungrounded entities
  and ungrounded relations.
- Plot frequency distribution of ungrounded entities
- Curate top XXX entities, categorize by type
- Repeat for grounded entities, identifying errors
- 

Methods (so far)

Step 1: Get Proteins->Genes->PMIDs

- Download Reactome signaling proteins list (see signaling_proteins_readme.txt)
  - output: signaling_proteins.tsv
- process_proteins.py
  - For each protein, get gene name.
  - output: signaling_genes.txt
  - Then get genes identified in bioentities/relations.csv (which included
    families from BEL and from the Ras family).
  - Combine the gene lists together
  - output: combined_genes.txt.
- get_pmids.py
  - For each gene in combined_genes.txt, get papers curated from Entrez gene.
  - Save dict mapping gene to list of papers
  - output: combined_genes_to_pmids.pkl
  - Make set of unique PMIDs, save
  - output: combined_pids.txt
  - (separately, make a subset of 2000 papers for pilot reading purposes:
     combined_pmids_subset.txt)

Step 2: Build REACH without Bioentities and Read PMIDs

- Build REACH, b4a284, without Bioentities (see email threads)
- Run on EC2, upload JSON files to key PMIDXXXXX/reach_no_bioentities

Step 3: Subsample the data into two files
- get_training_test_stmts.py
  - output: training_pmid_stmts.pkl
  - output: test_pmid_stmts.pkl

Step 4: stmt_entity_stats




