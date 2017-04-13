
if __name__ == '__main__':
    import pickle
    from indra.preassembler import grounding_mapper as gm
    from indra.tools.reading import reading_results_stats as rrs

    filename = '../step3_sample_training_test/training_pmid_stmts.pkl'
    #filename = '../step3_sample_training_test/test_pmid_stmts.pkl'
    with open(filename, 'rb') as f:
        pmid_stmts = pickle.load(f)
    print("No. of papers in training set: %d" % len(pmid_stmts.keys()))

    # Sort the statements by PMID key
    sorted_pmids = sorted(pmid_stmts.keys())
    stmts = []
    for pmid in sorted_pmids:
        for s in pmid_stmts[pmid]:
            stmts.append(s)

    gm.raw_agent_texts(stmts, 10000, 'raw_agent_texts_before.csv')

    rrs.report_grounding(stmts, bin_interval=100,
                         plot_prefix='training')

