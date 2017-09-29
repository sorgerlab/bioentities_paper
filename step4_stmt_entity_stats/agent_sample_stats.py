import os
import sys
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from indra.util import write_unicode_csv, plot_formatting as pf
from collections import defaultdict


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
    grounding_stats(curated)
