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
    """Return the number of times each FamPlex entry is grounded to."""
    counts = defaultdict(int)
    for stmt in stmts:
        for agent in stmt.agent_list():
            if agent is not None:
                be_id = agent.db_refs.get('FPLX')
                if be_id:
                    counts[be_id] += 1
    return counts


def get_hgnc_coverage_stats(stmts):
    """Return the number of times each HGNC gene is grounded to."""
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
    """Return list of FamPlex entries that weren't grounded to at all."""
    return set(entries) - set(counts.keys())


def plot_counts_by_entry(counts):
    """Plot the absolute counts and cumulative counts of groundings to
    FamPlex entries."""
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
    plt.xlabel('FamPlex entries')
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
    plt.xlabel('FamPlex entries')
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


def get_level_stats(entries, counts, hgnc_counts):
    """Calculate statistics of FamPlex entries that are grounded to at multiple
    levels of the hierarchy."""
    two_level_counts = {}
    multi_level_counts = {}
    for entry in entries:
        uri = 'http://sorger.med.harvard.edu/indra/entities/%s' % entry
        # If the entry has parents, we skip it (it will be added via its top
        # level entry)
        parents = eh.get_parents(uri)
        if parents:
            continue
        # If the entry does not have any children, then obviously any
        # grounding will be to the entry itself, so no level statistics
        # are needed.
        children = eh.get_children(uri)
        if not children:
            continue
        # Make a level map
        group = {'top': [entry], 'middle': [], 'bottom': []}
        for child in children:
            name = child.split('/')[-1]
            if 'hgnc.symbol' in child:
                group['bottom'].append(name)
            else:
                group['middle'].append(name)
        top_count = sum([counts.get(entry, 0) for entry in group['top']])
        bottom_count = sum([hgnc_counts.get(entry, 0) for entry in
                            group['bottom']])
        if not group['middle']:
            two_level_counts[entry] = (top_count, bottom_count)
        else:
            middle_count = sum([counts.get(entry, 0) for entry in
                               group['middle']])
            multi_level_counts[entry] = (top_count, middle_count, bottom_count)

    # Count the actual number of two-grounding / multi-grounding entries
    two_grounded = 0
    multi_grounded = 0
    for k, v in two_level_counts.items():
        if v[0] > 0 and v[1] > 0:
            two_grounded += 1
    for k, v in multi_level_counts.items():
        if v[0] > 0 and v[1] > 0 and v[2] > 0:
            multi_grounded += 1
        else:
            two_grounded += 1

    print(('Top-level FamPlex entries with groundings a 2 levels: %d '
           'and at least 3 levels: %d') % (two_grounded, multi_grounded))

    return two_level_counts, multi_level_counts


def get_stacks_groups(tops):
    """Return the entries at the various levels of the stack for a given list
    of top-level FamPlex entries."""
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
    """Plot stacks of a given group (entries subsumed by a FamPlex entry
    of interest showing statistics of number of groundings at various
    levels of the hierarchy."""
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
    bb = plt.bar(xticks, bottoms, color='#3182bd', align='center', linewidth=1)
    mb = plt.bar(xticks, middles, bottom=bottoms, color='#9ecae1',
                 align='center', linewidth=1)
    tb = plt.bar(xticks, tops, bottom=numpy.array(middles)+numpy.array(bottoms),
                 color='#deebf7', align='center', linewidth=1)
    plt.legend((tb[0], mb[0], bb[0]), ('Top (FamPlex)',
                                       'Middle (FamPlex)',
                                       'Bottom (specific gene)'),
                loc='upper right', frameon=False, fontsize=pf.fontsize)
    plt.xticks(xticks, labels)
    plt.ylabel('Number of times grounded to in test corpus')
    plt.subplots_adjust(left=0.15, right=0.97, top=0.95, bottom=0.08)
    ax = plt.gca()
    pf.format_axis(ax)
    plt.savefig('entity_levels_coverage_test_corpus.pdf')
    plt.show()


def plc_groundings(stmts, counts, hgnc_counts):
    """Print the number of times each PLC member was referred to as a
    percentage of all the groundings to 
    """
    group = get_stacks_groups(['PLC'])[0]
    allc = 0
    for gr, elements in group.items():
        for element in elements:
            if gr == 'bottom':
                count = hgnc_counts[element]
            else:
                count = counts[element]
            allc += count

    for gr, elements in group.items():
        for element in elements:
            if gr == 'bottom':
                count = hgnc_counts[element]
            else:
                count = counts[element]
            print('%s: %.2f' % (element, 100.0*count / allc))


if __name__ == '__main__':
    # Load Statements from test corpus reading output with FamPlex
    fname = '../step3_sample_training_test/famplex_test_stmts_mapped.pkl'
    with open(fname, 'rb') as fh:
        stmts = pickle.load(fh)

    # Load list of FamPlex entries
    entries = load_entity_list('../../famplex/entities.csv')

    # Get FamPlex counts grounded to in Statements
    counts = get_coverage_stats(stmts)
    # Get HGNC counts grounded to iun Statements
    hgnc_counts = get_hgnc_coverage_stats(stmts)

    # Quantify number of entries that are grounded to at multiple levels
    two_level_counts, multi_level_counts = \
        get_level_stats(entries, counts, hgnc_counts)

    # Print specific statistics about PLC and its children being grounded to
    plc_groundings(stmts, counts, hgnc_counts)

    # Get a list of entries that weren't grounded to at all 
    missing_entries = get_missing_entries(entries, counts)

    # Plot statistics of grounding at various levels for some specific groups
    groups_to_plot = ['AMPK', 'G_protein', 'PPP2', 'PLC', 'Activin']
    stacks_groups = get_stacks_groups(groups_to_plot)
    labels = [g.replace('_', ' ') for g in groups_to_plot]
    plt.ion()
    plot_stacks_groups(stacks_groups, counts, hgnc_counts, labels)
    plot_counts_by_entry(counts)
