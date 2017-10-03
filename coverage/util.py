import csv
import sys


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

