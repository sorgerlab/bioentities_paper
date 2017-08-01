import os
import sys
import csv
from indra.util import plot_formatting as pf
import matplotlib.pyplot as plt
from matplotlib_venn import venn3

bepath = '../../../bioentities'

def read_csv(fh, delimiter, quotechar):
    if sys.version_info.major < 3:
        csvreader = csv.reader(fh, delimiter=bytes(delimiter),
                               quotechar=bytes(quotechar))
    else:
        csvreader = csv.reader(fh, delimiter=delimiter, quotechar=quotechar)
    rows = [row for row in csvreader]
    return rows

def load_csv(filename):
    with open(filename) as f:
        rows = read_csv(f, ',', '"')
    return rows

def load_entity_list(filename):
    with open(filename) as f:
        rows = read_csv(f, ',', '"')
    entities = [row[0] for row in rows]
    return entities

def load_equivalences(filename):
    equivalences = []
    with open(filename) as f:
        rows = read_csv(f, ',', '"')
        for row in rows:
            equivalences.append((row[0], row[1], row[2]))
    return equivalences

def find_db_mappings(entity, equivalences):
    dbs = []
    for db, db_id, be_id in equivalences:
        if be_id == entity:
            dbs.append(db)
    return dbs

def plot_venn_diagram(db_mappings, groups):
    missing = len([entity for entity, db_eq in db_mappings.items()
                   if not db_eq])
    group_entries = []
    for idx, group in enumerate(groups):
        entries = []
        for entity, db_eq in db_mappings.items():
            if set(db_eq) & set(group):
                entries.append(entity)
        group_entries.append(entries)
    subsets = (len(set(group_entries[0]) - set(group_entries[1])),
               len(set(group_entries[1]) - set(group_entries[0])),
               missing,
               len(set(group_entries[0]) & set(group_entries[1])),
               0,
               0,
               0)
    pf.set_fig_params()
    plt.figure(figsize=(4, 3))
    venn3(subsets=subsets, set_labels=('BEL', 'Other', 'Unmapped'))
    plt.savefig('bioentities_mapping.pdf')
    ax = plt.gca()
    pf.format_axis(ax)

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
    bel_only = sorted([entity for entity, db_eq in db_mappings.items()
                       if db_eq == ['BEL']])

    print('Equivalences coverage\n=====================')
    print('Missing: %d (%.0f%%)' % (len(missing),
                                    100.0*len(missing)/len(entities)))
    print('BEL only: %d (%.0f%%)' % (len(bel_only),
                                     100.0*len(bel_only)/len(entities)))

    groups = [['BEL'], ['IP', 'PF', 'RE', 'NXP', 'NCIT', 'GO']]
    plot_venn_diagram(db_mappings, groups)
