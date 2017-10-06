import csv
import numpy
import pickle
import matplotlib.pyplot as plt
from collections import defaultdict
from indra.util import plot_formatting as pf
from indra.databases import hgnc_client
from indra.preassembler.hierarchy_manager import entity_hierarchy as eh
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


def get_hgnc_coverage_stats(stmts):
    counts = defaultdict(int)
    for stmt in stmts:
        for agent in stmt.agent_list():
            if agent is not None:
                hgnc_id = agent.db_refs.get('HGNC')
                if hgnc_id:
                    hgnc_name = hgnc_client.get_hgnc_name(hgnc_id)
                    if hgnc_name:
                        counts[hgnc_name] += 1
    return counts


def get_missing_entries(entries, counts):
    return set(entries) - set(counts.keys())


def plot_counts_by_entry(counts):
    pf.set_fig_params()
    plt.figure(figsize=(2.5, 2.5), dpi=300)
    counts_ord = sorted(counts.items(), key=lambda x: x[1], reverse=True)
    names = [cc[0] for cc in counts_ord]
    counts_for_name = numpy.array([float(cc[1]) for cc in counts_ord])
    vals = numpy.cumsum(counts_for_name) / numpy.sum(counts_for_name)
    idx50 = [i for i, v in enumerate(vals) if v >= 0.5][0]

    # Plot absolute counts
    xvals = range(len(names))
    xvals_top, xvals_rest = xvals[:idx50], xvals[idx50:]
    yvals_top, yvals_rest = counts_for_name[:idx50], counts_for_name[idx50:] 
    plt.bar(xvals_top, yvals_top, color=pf.PURPLE, linewidth=0)
    plt.bar(xvals_rest, yvals_rest, color=pf.GREEN, linewidth=0)
    plt.xlim([0, len(names)])
    plt.ylabel('Number of times grounded to in test corpus')
    plt.xlabel('Bioentities entries')
    plt.yscale('log')
    ax = plt.gca()
    pf.format_axis(ax)
    plt.subplots_adjust(left=0.14, bottom=0.11, top=0.93, right=0.95)
    plt.savefig('entity_coverage_test_corpus.pdf')
    plt.show()

    # Plot cumulative percentage
    print('50%% of mentions is covered by the first %d entries' % (idx50))
    plt.figure(figsize=(2.5, 2.5), dpi=300)
    plt.plot(range(len(names)), vals, color=pf.GREEN)
    plt.plot([0, idx50], [0.5, 0.5], color='black', linestyle='dashed')
    plt.plot([idx50, idx50], [0, 0.5], color='black', linestyle='dashed')
    plt.xlim([0, len(names)])
    plt.ylabel('Cumulative rel. frequency of times \n'
               'grounded to in test corpus')
    plt.xlabel('Bioentities entries')
    ax = plt.gca()
    pf.format_axis(ax)
    plt.subplots_adjust(left=0.19, bottom=0.11, top=0.93, right=0.95)
    plt.savefig('entity_coverage_test_corpus_distr.pdf')
    plt.show()

    # Print top table
    ntop = 5
    for name, count in counts_ord[:ntop]:
        print('%s & %d & %.1f\\\\' % (name, count,
                                      100* count / numpy.sum(counts_for_name)))


def get_stacks_groups(tops):
    groups = []
    for top in tops:
        uri = 'http://sorger.med.harvard.edu/indra/entities/%s' % top
        children = eh.get_children(uri)
        group = {'top': [top], 'middle': [], 'bottom': []}
        for child in children:
            name = child.split('/')[-1]
            if 'hgnc.symbol' in child:
                group['bottom'].append(name)
            else:
                group['middle'].append(name)
        groups.append(group)
    return groups


def plot_stacks_groups(stacks_groups, counts, hgnc_counts, labels):
    plt.figure
    tops = []
    middles = []
    bottoms = []
    for group in stacks_groups:
        top_count = sum([counts.get(entry, 0) for entry in group['top']])
        middle_count = sum([counts.get(entry, 0) for entry in group['middle']])
        bottom_count = sum([hgnc_counts.get(entry, 0) for entry in
                            group['bottom']])
        tops.append(top_count)
        middles.append(middle_count)
        bottoms.append(bottom_count)
    xticks = numpy.arange(len(stacks_groups))
    pf.set_fig_params()
    plt.figure(figsize=(3, 2.5), dpi=300)
    bb = plt.bar(xticks, bottoms, color='#3182bd', align='center')
    mb = plt.bar(xticks, middles, bottom=bottoms, color='#9ecae1',
                 align='center')
    tb = plt.bar(xticks, tops, bottom=numpy.array(middles)+numpy.array(bottoms),
                 color='#deebf7', align='center')
    plt.legend((tb[0], mb[0], bb[0]), ('Top (Bioentities)',
                                       'Middle (Bioentities)',
                                       'Bottom (specific gene)'),
                loc='upper right', frameon=False, fontsize=pf.fontsize)
    plt.xticks(xticks, labels)
    plt.ylabel('Number of times grounded to in test corpus')
    plt.subplots_adjust(left=0.15, right=0.97, top=0.95, bottom=0.08)
    ax = plt.gca()
    pf.format_axis(ax)
    plt.savefig('entity_levels_coverage_test_corpus.pdf')
    plt.show()


def plc_groundings(stmts):
    group = get_stacks_groups(['PLC'])[0]
    for element in group['bottom']:
        print(element)
        print('=====================')
        for stmt in stmts:
            for agent in stmt.agent_list():
                if agent is not None:
                    be_id = hgnc_client.get_hgnc_name(agent.db_refs.get('HGNC'))
                    if be_id == element:
                        print(stmt.evidence[0].pmid, stmt.evidence[0].text)


if __name__ == '__main__':
    fname = '../step3_sample_training_test/bioentities_test_stmts_mapped.pkl'
    with open(fname, 'rb') as fh:
        stmts = pickle.load(fh)

    entries = load_entity_list('../../../bioentities/entities.csv')

    counts = get_coverage_stats(stmts)
    hgnc_counts = get_hgnc_coverage_stats(stmts)
    missing_entries = get_missing_entries(entries, counts)
    groups_to_plot = ['AMPK', 'G_protein', 'PPP2', 'PLC', 'Activin']
    stacks_groups = get_stacks_groups(groups_to_plot)
    labels = [g.replace('_', ' ') for g in groups_to_plot]
    plt.ion()
    plot_stacks_groups(stacks_groups, counts, hgnc_counts, labels)
    plot_counts_by_entry(counts)

    plc_groundings(stmts)
