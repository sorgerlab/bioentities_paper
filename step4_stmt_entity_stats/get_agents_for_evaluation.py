import os
import sys
import pickle
from indra.preassembler import grounding_mapper as gm
from indra.tools.reading import reading_results_stats as rrs

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: %s [training|test]' % os.path.basename(__file__))
        sys.exit()
    mode = sys.argv[1]
    if mode not in ('training', 'test'):
        print('Usage: %s [training|test]' % os.path.basename(__file__))
        sys.exit()

    filename = '../step3_sample_training_test/%s.pkl' % \
            ('training_pmid_stmts' if mode == 'training' else \
            'bioentities_test_stmts_mapped')
    with open(filename, 'rb') as f:
        stmts = pickle.load(f)

    if mode == 'training':
        print("No. of papers in %s set: %d" % len(stmts.keys()))

        # Sort the statements by PMID key
        sorted_pmids = sorted(stmts.keys())
        stmts_list = []
        for pmid in sorted_pmids:
            for s in stmts[pmid]:
                stmts_list.append(s)
        stmts = stmts_list

    # Get 10000 randomly selected agents from the full statement list
    out_file = 'training_agents_sample.csv' if mode == 'training' else \
        'test_agents_with_be_sample.csv'
    gm.raw_agent_texts(stmts, 10000, out_file)

    plot_prefix = 'training' if mode == 'training' else 'test_with_be'
    rrs.report_grounding(stmts, bin_interval=100,
                         plot_prefix=plot_prefix)

