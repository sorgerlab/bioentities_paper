import pickle
from indra.preassembler.grounding_mapper import *
from fuzzywuzzy import fuzz
#from process_proteins import get_relations_genes

# Take a given gene symbol

# Look up synonyms (in HGNC? Uniprot?)

# Search the list of entities/freqs for that synonym

def load_texts(stmts_filename):
    with open(stmts_filename, 'rb') as f:
        stmts_by_paper = pickle.load(f)
    stmts = [s for stmt_list in stmts_by_paper.values() for s in stmt_list]

    texts = agent_texts_with_grounding(stmts)
    return texts

def get_keyword_matches(kw, texts, match_type='partial'):
    hits = []
    for entry in texts:
        text = entry[0]
        if match_type == 'partial':
            ratio = fuzz.partial_ratio(kw.upper(), text.upper())
        else:
            ratio = fuzz.ratio(kw.upper(), text.upper())
        hits.append((entry, ratio))
    hits.sort(key=lambda x: x[1], reverse=True)
    return hits

if __name__ == '__main__':
    texts = load_texts('../step3_sample_training_test/training_pmid_stmts.pkl')
    #with open('../entities.csv', 'rt') as f:
    #    bioents = [line.strip() for line in f.readlines()]
    def getp(kw):
        return get_keyword_matches(kw, texts, match_type='partial')
    def getr(kw):
        return get_keyword_matches(kw, texts, match_type='ratio')

    # Workflow:
    #
    # 1. Finish reading the 270k papers with the latest REACH
    # 1. Update to get abstracts when text can't be obtained from Elsevier
    # 1. Update existing S3 REACH keys with new structure
    # 2. Recompile REACH per Mihai's instructions
    # 3. Re-read papers with no Bioentities version
    # 4. Get statements, partition randomly into 90% and 10% sets, by paper
    # 5. Plot top list of ungrounded entities
    # 6. Plot top list of entities, by frequency, noting misgrounded ones
    # 7. Estimate grounding error rate on all entities and all families

    # 8. Run curation procedure on 90% set:
    #   - For each bioentity (of the 400 or so), and also for all of the
    #     families (complexes?) defined in Reactome, do a search of the texts
    #   - Dump out a text file with the candidate texts, csv-separated with
    #     the BE identifier used to search
    #   - Add texts to grounding map for

    #
    #    How to store multiple REACH versions? Change to use key like
    #    PMID1233556/reach/1.3.3-asdf
    #    Each one would have the source material used in the head
    #    Would require updating all keys?? Not necessarily--

