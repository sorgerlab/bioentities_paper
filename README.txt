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
- get_agents_for_evaluation.py
  - output: training_agents_sample.csv
  - output: training_agents.csv
  - output: training_ungrounded.csv
  - output: training_agent_distribution.pdf
  - output: training_ungrounded_distribution.pdf
- Open training_agents_sample_curated.csv in Excel to curate
  - Curate agents. Code:
    - S: small molecule
    - P: protein
    - B: biological process
    - F: family
    - C: complex
    - X: complex of families
    - M: microRNA
    - U: unknown
  - JAB curated rows 2-101; BMG 102-251; JAB 252-301 (row 1 is header)
- Estimate percentage of entitities that are family, complex, or combined
  - family_frequency.py
  - Analyzes training_agents_sample_curated.csv, plots percentage


  - Fraction of most frequent agents that in
clude family or complex
    - Fraction of mentions
  - Fraction of most frequent ungrounded that include family or complex
    - Fraction of ungrounded mentions
  - Table with 3 columns Entity, Frequency, Category (F, C, X, or empty),
    Curator


