if __name__ == '__main__':
    import pandas as pd
    import numpy as np
    from matplotlib import pyplot as plt
    from indra.util import plot_formatting as pf
    pf.set_fig_params()

    with open('training_agents.xlsx', 'rb') as f:
        data = pd.read_excel(f, sheetname=0)

    family_cats = ('F', 'C', 'X')
    # Get the first 100 rows
    num_curated = 200
    curated = data[0:num_curated]
    # Iterate over rows and calculate rolling ratio due to families/complexes
    totals = np.zeros(num_curated)
    families = np.zeros(num_curated)
    for ix, row in curated.iterrows():
        cat = row['CategoryFinal']
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

    plt.ion()
    index = range(1, num_curated + 1)
    fig = plt.figure(figsize=(3, 3), dpi=150)
    ax = fig.gca()
    ax.plot(index, totals, 'r')
    ax.plot(index, families, 'b')
    pf.format_axis(ax)
    ax.set_ylabel('Number of occurrences')
    ax.set_xlabel('Entity index')


    ratio = families / totals
    fig = plt.figure(figsize=(3, 3), dpi=150)
    ax = fig.gca()
    ax.plot(index, ratio)
    pf.format_axis(ax)
    ax.set_xlabel('Entity index')
    ax.set_ylabel('Pct. family/complex occurrences')
