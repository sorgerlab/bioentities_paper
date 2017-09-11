import os
import sys
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from indra.util import write_unicode_csv, plot_formatting as pf
from collections import defaultdict


def plot_rolling_ratio(curated):
    # Iterate over rows and calculate rolling ratio due to families/complexes
    totals = np.zeros(num_curated)
    families = np.zeros(num_curated)
    categories = defaultdict(list)
    for ix, row in curated.iterrows():
        cat = row['EntityType']
        categories[cat].append(1)
        if ix == 0:
            totals[ix] = 1
            if cat in family_cats:
                families[ix] = 1
        else:
            totals[ix] = totals[ix - 1] + 1
            if cat in family_cats:
                families[ix] = families[ix - 1] + 1
            else:
                families[ix] = families[ix - 1]

    total = 0
    cat_summary = {}
    # Iterate to get the total
    total_counts = np.sum([np.sum(freqs)
                           for freqs in categories.values()])
    for cat, freqs in categories.items():
        cat_freq = np.sum(freqs)
        cat_summary[cat] = (len(freqs), cat_freq,
                            (100 * cat_freq) / float(total_counts))
        print('%s: %d instances, %d occurrences (%.2f%%)' %
              (cat, len(freqs), cat_freq,
               (100 * cat_freq) / float(total_counts)))

    plt.ion()
    index = range(1, num_curated + 1)
    fig = plt.figure(figsize=(3, 3), dpi=150)
    ax = fig.gca()
    ax.plot(index, totals, 'r')
    ax.plot(index, families, 'b')
    pf.format_axis(ax)
    ax.set_ylabel('Number of occurrences')
    ax.set_yscale('log')
    ax.set_xlabel('Entity index')
    plt.subplots_adjust(left=0.21)

    ratio = families / totals
    fig = plt.figure(figsize=(3, 3), dpi=150)
    ax = fig.gca()
    ax.plot(index, ratio)
    pf.format_axis(ax)
    ax.set_xlabel('Entity index')
    ax.set_ylabel('Pct. family/complex occurrences')
    plt.subplots_adjust(left=0.16)


def grounding_stats(data):
    cats = (['P'], ['F', 'C', 'X'], ['S'], ['B'], ['U'], ['M'])
    cat_names = ('Protein/gene', 'Family/complex', 'Small molecule',
                 'Biological process', 'Other/unknown', 'microRNA')
    rows = []
    num_agents = len(data)
    plt.figure(figsize=(2, 2), dpi=300)
    for ix, cat in enumerate(cats):
        cat_rows = data[data.EntityType.apply(lambda et: et in cat)]
        cat_number = len(cat_rows)
        cat_pct = (100 * cat_number / float(num_agents))
        cat_pct_str = '%.1f' % cat_pct
        correct_rows = cat_rows[cat_rows.Grounding == 1]
        correct_number = len(correct_rows)
        correct_pct = (100 * correct_number / float(cat_number)) if \
            cat_number > 0 else 0
        correct_pct_of_total = (100 * correct_number) / float(num_agents)
        correct_pct_str = '%.1f' % correct_pct
        rows.append((cat, cat_number, cat_pct_str, correct_number,
                     correct_pct_str))
        def stderr(k, n):
            return np.sqrt(((k/float(n)) * (1-(k/float(n)))) / float(n))
        stderr_inc = 100 * stderr(cat_number - correct_number, num_agents)
        inc_handle = plt.bar(ix, cat_pct, color='red', align='center',
                             yerr=stderr_inc, linewidth=0.5)
        stderr_corr = 100 * stderr(correct_number, num_agents)
        corr_handle = plt.bar(ix, correct_pct_of_total, color='blue',
                              align='center', yerr=stderr_corr,
                              linewidth=0.5)
    plt.xticks(range(len(cats)), cat_names, rotation=90)
    plt.ylabel('Pct. Curated Entities')
    plt.subplots_adjust(left=0.18, bottom=0.43, top=0.96)
    ax = plt.gca()
    pf.format_axis(ax)
    plt.legend((corr_handle, inc_handle), ('Correct', 'Incorrect'),
               loc='upper right', frameon=False, fontsize=pf.fontsize)
    plt.show()
    print(rows)
    write_unicode_csv('%s_agents_sample_stats.csv' % mode, rows)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: %s [training|test]' % os.path.basename(__file__))
        sys.exit()
    mode = sys.argv[1]
    if mode not in ('training', 'test'):
        print('Usage: %s [training|test]' % os.path.basename(__file__))
        sys.exit()


    pf.set_fig_params()
    plt.ion()
    family_cats = ('F', 'C', 'X')

    fname = 'training_agents_sample_curated.csv' if mode == 'training' else \
            'test_agents_with_be_sample_curated.csv'
    with open(fname, 'rb') as f:
        data = pd.read_csv(f)

    num_curated = 300
    curated = data[0:num_curated]
    #plot_rolling_ratio(curated)
    grounding_stats(curated)
