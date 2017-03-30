import random
import pickle

# Load dict of PMIDs->statements
with open('../step2_read_no_bioentities/combined_genes_no_bioentities_stmts.pkl',
          'rb') as f:
    print("Loading statements")
    pmid_stmts = pickle.load(f)

pmids = list(pmid_stmts.keys())

training_fraction = 0.9
partition = int(training_fraction * len(pmids))
training_pmids = pmids[0:partition]
test_pmids = pmids[partition:]

training_stmts = {pmid: pmid_stmts[pmid] for pmid in training_pmids}
test_stmts = {pmid: pmid_stmts[pmid] for pmid in test_pmids}

with open('training_pmid_stmts.pkl', 'wb') as f:
    print("Saving training stmts")
    pickle.dump(training_stmts, f)

with open('test_pmid_stmts.pkl', 'wb') as f:
    print("Saving test stmts")
    pickle.dump(test_stmts, f)





