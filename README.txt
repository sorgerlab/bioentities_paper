Methods (so far)

Step 1: Get Proteins->Genes->PMIDs

- Download Reactome signaling proteins list
  (see step1_genes_pmids/signaling_proteins_readme.txt)
  - output: signaling_proteins.tsv
- Run process_proteins.py
  - For each protein, get gene name.
  - output: signaling_genes.txt
- Get genes identified in bioentities/relations.csv (which included
    families from BEL and from the Ras family).
  - Combine the gene lists together
  - output: combined_genes.txt.
- get_pmids.py
  - For each gene in combined_genes.txt, get papers curated from Entrez gene.
  - Save dict mapping gene to list of papers
  - output: combined_genes_to_pmids.pkl
  - Make set of unique PMIDs, save
  - output: combined_pids.txt

Step 2: Build REACH without Bioentities and Read PMIDs

- Build REACH, b4a284, without Bioentities (see email threads)
- Run on EC2, upload JSON files to key PMIDXXXXX/reach_no_bioentities
  - output: combined_genes_no_bioentities_stmts.pkl (pmid->stmts file, Python3)

Step 3: Shuffle the PMIDs and subsample the data into two files
- 80% of papers to training set, 20% of papers to test set
- get_training_test_stmts.py
  - output: training_pmid_stmts.pkl
  - output: test_pmid_stmts.pkl

Step 4: stmt_entity_stats
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
  - JAB curated rows 2-101; BMG 102-251; JAB 252-301 (row 1 is header)
- For ease of curation, run generate_agent_links.py to create an HTML table
  to check groundings in different databases.
- Same for test set.

Step 5: Curate grounding map
  - texts_for_gene.py helps in finding lexical synonyms for unmapped families.
  - Fraction of most frequent agents that in
clude family or complex
    - Fraction of mentions
  - Fraction of most frequent ungrounded that include family or complex
    - Fraction of ungrounded mentions
  - Table with 3 columns Entity, Frequency, Category (F, C, X, or empty),
    Curator

Measure coverage
- test_entities_coverage.py:
  - Creates a plot of frequency of groundings to each entity
  -
