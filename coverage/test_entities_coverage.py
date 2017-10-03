import csv
import pickle
from collections import defaultdict
from util import *

def get_coverage_stats(stmts):
    counts = defaultdict(int)
    for stmt in stmts:
        for agent in stmt.agent_list():
            if agent is not None:
                be_id = agent.db_refs.get('BE')
                if be_id:
                    counts[be_id] += 1
    return counts


def get_missing_entries(entries, counts):
    return set(entries) - set(counts.keys())


if __name__ == '__main__':
    fname = '../step3_sample_training_test/bioentities_test_stmts_mapped.pkl'
    with open(fname, 'rb') as fh:
        stmts = pickle.load(fh)

    entries = load_entity_list('../../../bioentities/entities.csv')

    counts = get_coverage_stats(stmts)
    missing_entries = get_missing_entries(entries, counts)
