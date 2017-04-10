if __name__ == '__main__':
    import pandas as pd
    import numpy as np
    from matplotlib import pyplot as plt
    from indra.util import plot_formatting as pf
    from collections import defaultdict

    pf.set_fig_params()

    with open('training_agents.xlsx', 'rb') as f:
        data = pd.read_excel(f, sheetname=0)

    family_cats = ('F', 'C', 'X')
    num_curated = 310
    curated = data[0:num_curated]
    # Iterate over rows and calculate rolling ratio due to families/complexes
    totals = np.zeros(num_curated)
    families = np.zeros(num_curated)
    categories = defaultdict(list)
    for ix, row in curated.iterrows():
        cat = row['CategoryFinal']
        categories[cat].append(row['Freq'])
        if ix == 0:
            totals[ix] = row['Freq']
            if cat in family_cats:
                families[ix] = row['Freq']
        else:
            totals[ix] = row['Freq'] + totals[ix - 1]
            if cat in family_cats:
                families[ix] = families[ix - 1] + row['Freq']
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

