"""This script compares the frequency of entities that
appear in the test set with those in the training set, and
produces plots that show the comparison."""
import os
import pickle
import numpy
from collections import Counter
import matplotlib.pyplot as plt


def get_strings(stmts):
    """Return entitiy strings from a list of Statements."""
    strs = []
    for stmtl in stmts.values():
        for stmt in stmtl:
            for agent in stmt.agent_list():
                if agent is not None:
                    txt = agent.db_refs.get('TEXT')
                    if txt is not None:
                        strs.append(txt)
    return strs


def get_counts(fname):
    """Return entitiy counts from a pickle file of Statements."""
    print('Loading %s' % fname)
    with open(fname, 'rb') as fh:
        stmts = pickle.load(fh)
        strs = get_strings(stmts)
        counts = Counter(strs)
        countsl = sorted(list(counts.items()),
                         key=lambda x: x[1], reverse=True)
        return countsl


def align_lists(l1, l2):
    """Align two lists so that they share the same joint ordering."""
    d = {}
    for k, v in l1:
        d[k] = [v, 0]
    for k, v in l2:
        try:
            d[k][1] = v
        except KeyError:
            d[k] = [0, v]
    dsort = sorted(d.values(), key=lambda x: x[0], reverse=True)
    ll1, ll2 = [[e1 for e1, e2 in dsort],
                [e2 for e1, e2 in dsort]]
    ll1 = numpy.array(ll1) / numpy.sum(ll1)
    ll2 = numpy.array(ll2) / numpy.sum(ll2)
    return ll1, ll2


def plot_hist_subfig(l1, l2, nents):
    plt.figure(figsize=(12,10))
    for i, nent in enumerate(nents):
        plt.subplot(221 + i)
        plt.plot(l2[:nent], color='blue', label='Test')
        plt.plot(l1[:nent], color='red', label='Training')
        plt.ylim([-0.0005, 0.018])
        plt.title('Top %d entities' % nent)
        plt.ylabel('Realtive frequency of occurrence')
        plt.xlabel('Entity (in order of frequency of occurrence in training set)')
        plt.legend()
    plt.savefig('entity_freqs_train_test.png')


if __name__ == '__main__':
    # Cache the counts pickle because getting counts is fairly slow
    if os.path.exists('counts.pkl'):
        with open('counts.pkl', 'rb') as fh:
            tr_counts, te_counts = pickle.load(fh)
    else:
        tr_counts = get_counts('training_pmid_stmts.pkl')
        te_counts = get_counts('test_pmid_stmts.pkl')
        with open('counts.pkl', 'wb') as fh:
            pickle.dump([tr_counts, te_counts], fh)

    # Align the two lists obtained from training and test sets
    trl, tel = align_lists(tr_counts, te_counts)
    # Plot frequencies as subfigures
    plot_hist_subfig(trl, tel, [10, 100, 1000, 10000])
