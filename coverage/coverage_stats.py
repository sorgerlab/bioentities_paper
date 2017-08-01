import os
import sys
import csv

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

if __name__ == '__main__':
    # Read all entities in Bioentities
    entities_file = os.path.join(bepath, 'entities.csv')
    entities = load_entity_list(entities_file)
    # Read equivalence table
    equivalences_file = os.path.join(bepath, 'equivalences.csv')
    equivalences = load_equivalences(equivalences_file)
    # Find equivalence DB types for each term
    dbs = {entity: find_db_mappings(entity, equivalences)
                   for entity in entities}

    missing = sorted([entity for entity, db_eq in dbs.items() if not db_eq])
    bel_only = sorted([entity for entity, db_eq in dbs.items()
                       if db_eq == ['BEL']])

    print('Equivalences coverage\n=====================')
    print('Missing: %d (%.0f%%)' % (len(missing),
                                100.0*len(missing)/len(entities)))
    print('BEL only: %d (%.0f%%)' % (len(bel_only),
                                 100.0*len(bel_only)/len(entities)))
