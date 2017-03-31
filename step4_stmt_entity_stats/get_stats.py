
if __name__ == '__main__':
    import pickle
    from indra.tools.reading import reading_results_stats as rrs

    filename = '../step3_sample_training_test/training_pmid_stmts.pkl'
    with open(filename, 'rb') as f:
        pmid_stmts = pickle.load(f)
    print("No. of papers in training set: %d" % len(pmid_stmts.keys()))

    stmts = [s for stmt_list in pmid_stmts.values() for s in stmt_list]
    rrs.report_grounding(stmts, list_length=1000, bin_interval=100,
                         plot_prefix='training')

