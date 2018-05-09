Scripts to reproduce results reported in FamPlex manuscript
===========================================================

This repository is organized into six folders corresponding to each stage of
constructing and evaluating the FamPlex resource.

Step 1: Get list of genes and corresponding PMIDs
-------------------------------------------------

- Download Reactome signaling proteins list
  (see step1_genes_pmids/signaling_proteins_readme.txt)
  - output: signaling_proteins.tsv
- Run process_proteins.py to get a gene list for literature search
  - For each protein, get gene name.
  - output: signaling_genes.txt
  - Get genes identified in famplex/relations.csv (which included
    families from BEL and from the Ras family).
  - Combine the two gene lists
  - output: combined_genes.txt.
- Run get_pmids.py to get a list of PMIDs from gene list
  - For each gene in combined_genes.txt, get papers curated from Entrez gene.
  - Save dict mapping gene to list of papers
  - output: combined_genes_to_pmids.pkl
  - Make set of unique PMIDs, save
  - output: combined_pids.txt

Step 2: Build REACH without FamPlex and read PMIDs
--------------------------------------------------

- Build REACH without FamPlex
- output: combined_genes_no_famplex_stmts.pkl (the pickle file contains a dictionary
  mapping each PMID to a list of Statements)


Step 3: Shuffle the PMIDs and subsample the data into training and test sets
----------------------------------------------------------------------------
- 80% of papers to training set, 20% of papers to test set
- get_training_test_stmts.py
  - output: training_pmid_stmts.pkl
  - output: test_pmid_stmts.pkl
  - output: training_pmids.txt
  - output: test_pmids.txt

Step 4: stmt_entity_stats
-------------------------
- Run get_agents_for_evaluation.py. Creates random sample of agents from
  training (or test) set for evaluation.
  - output: training_agents_sample.csv
  - output: training_agents.csv
  - output: training_ungrounded.csv
  - output: training_agent_distribution.pdf
  - output: training_ungrounded_distribution.pdf
- Open training_agents_sample_curated.csv in Excel to curate
  - Curate agents. Code:
    - P: protein
    - F: family
    - C: complex
    - X: complex of families
    - S: small molecule
    - B: biological process
    - U: unknown/other
    - M: microRNA
- For ease of curation, run generate_agent_links.py to create an HTML table
  to check groundings in different databases.
- Same for test set.

Step 5: Curate grounding map
----------------------------
  - texts_for_gene.py helps in finding lexical synonyms for unmapped families.
  - Fraction of most frequent agents that in
clude family or complex
    - Fraction of mentions
  - Fraction of most frequent ungrounded that include family or complex
    - Fraction of ungrounded mentions
  - Table with 3 columns Entity, Frequency, Category (F, C, X, or empty),
    Curator

Step 6: Evaluate FamPlex resource
---------------------------------
- Creates a plot of frequency of groundings to each entity
 
