import os
import sys
import csv
from indra.util import plot_formatting as pf
import matplotlib.pyplot as plt
from matplotlib_venn import venn3

bepath = '../../../../bioentities'

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
    plt.savefig('bioentities_mapping.pdf')
    ax = plt.gca()
    pf.format_axis(ax)
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
