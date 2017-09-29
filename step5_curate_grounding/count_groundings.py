from indra.util import read_unicode_csv, plot_formatting as pf
from os.path import join
from collections import Counter
from matplotlib import pyplot as plt

pf.set_fig_params()

be_path = join('..', '..', '..', 'bioentities')

# Open entities list
be_entities = [row[0] for row in read_unicode_csv(join(be_path, 'entities.csv'))]
be_counts = dict([(be, 0) for be in be_entities])
# Count groundings for each entity
for row in read_unicode_csv(join(be_path, 'grounding_map.csv')):
    if row[1] == 'BE':
        be_counts[row[2]] += 1
    assert 'BE' not in row[2:]
be_ctr = Counter(be_counts.values())

# Plot
plt.ion()
plt.figure(figsize=(3, 2), dpi=150)
for num, freq in be_ctr.items():
    plt.bar(num, freq, width=0.8, align="center", color='gray')
ax = plt.gca()
max_count = max(be_ctr.keys())
xticks = range(1, max_count + 1)
ax.set_xticks(xticks)
ax.set_xlabel('Number of synonyms')
ax.set_ylabel('Count')
pf.format_axis(ax, tick_padding=2)
plt.subplots_adjust(bottom=0.15)
plt.savefig('num_groundings.pdf')
