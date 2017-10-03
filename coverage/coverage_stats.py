import os
import sys
import numpy
from indra.util import plot_formatting as pf
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
from matplotlib.patches import Circle
from util import *

bepath = '../../../../bioentities'

def find_db_mappings(entity, equivalences):
    dbs = []
    for db, db_id, be_id in equivalences:
        if be_id == entity:
            dbs.append(db)
    return dbs

def get_entries_by_group(db_mappings, groups):
    missing = len([entity for entity, db_eq in db_mappings.items()
                   if not db_eq])
    group_entries = []
    for idx, group in enumerate(groups):
        entries = []
        for entity, db_eq in db_mappings.items():
            if set(db_eq) & set(group):
                entries.append(entity)
        group_entries.append(entries)
    return group_entries, missing

def plot_venn_diagram(group_entries, missing):
    subsets = [set(g) for g in group_entries]
    pf.set_fig_params()
    plt.figure(figsize=(4, 3))
    v3 = venn3(subsets=subsets, set_labels=('Pfam / InterPro / NextProt / GO',
                                            'NCIT / MeSH', 'BEL / Reactome'))
    radius_missing = v3.radii[0] * numpy.sqrt(1.0*missing / len(subsets[0]))
    center_missing = [0.8, 0]
    circle_missing = Circle(center_missing, radius_missing, color='black',
                            alpha=0.6)
    ax = plt.gca()
    ax.add_patch(circle_missing)
    ax.text(center_missing[0]-0.2, center_missing[1]+0.1, 'No mapping')
    ax.text(center_missing[0], center_missing[1], str(missing))
    plt.xlim((-1, 1))
    pf.format_axis(ax)
    plt.savefig('bioentities_mapping.pdf')
    return v3

if __name__ == '__main__':
    # Read all entities in Bioentities
    entities_file = os.path.join(bepath, 'entities.csv')
    entities = load_entity_list(entities_file)
    # Read equivalence table
    equivalences_file = os.path.join(bepath, 'equivalences.csv')
    equivalences = load_equivalences(equivalences_file)
    # Find equivalence DB types for each term
    db_mappings = {entity: find_db_mappings(entity, equivalences)
                   for entity in entities}

    missing = sorted([entity for entity, db_eq in db_mappings.items()
                      if not db_eq])

    print('Equivalences coverage\n=====================')
    print('Missing: %d (%.0f%%)' % (len(missing),
                                    100.0*len(missing)/len(entities)))

    groups = [['IP', 'PF','NXP', 'GO'], ['NCIT', 'MESH'], ['BEL', 'RE']]
    group_entries, missing = get_entries_by_group(db_mappings, groups)
    v3 = plot_venn_diagram(group_entries, missing)
