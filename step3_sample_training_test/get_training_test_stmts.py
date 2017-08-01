import random
import pickle

# Load dict of PMIDs->statements
with open('../step2_read_no_bioentities/'
          'combined_genes_no_bioentities_stmts.pkl', 'rb') as f:
    print("Loading statements")
    pmid_stmts = pickle.load(f)

# Get the PMIDs and sort
pmids = sorted(list(pmid_stmts.keys()))
# Seed the random number generator for reproducibility
random.seed(1)
# Shuffle the PMIDs
random.shuffle(pmids)

# Partition the shuffled PMIDs into training/test sets
training_fraction = 0.8
partition = int(training_fraction * len(pmids))
training_pmids = pmids[0:partition]
test_pmids = pmids[partition:]

# Get the statements for the training/test sets
training_stmts = {pmid: pmid_stmts[pmid] for pmid in training_pmids}
test_stmts = {pmid: pmid_stmts[pmid] for pmid in test_pmids}

# Save the training/test sets to pickle files
with open('training_pmid_stmts.pkl', 'wb') as f:
    print("Saving training stmts")
    pickle.dump(training_stmts, f)
with open('test_pmid_stmts.pkl', 'wb') as f:
    print("Saving test stmts")
    pickle.dump(test_stmts, f)





