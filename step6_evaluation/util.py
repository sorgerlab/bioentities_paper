"""This file contains various utility functions for other scripts to use,
including loading various tables of FamPlex."""

import csv
import sys
from collections import defaultdict


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


def load_relationships(filename):
    relationships = []
    with open(filename) as f:
        rows = read_csv(f, ',', '"')
    for row in rows:
        relationships.append(((row[0], row[1]), row[2], (row[3], row[4])))
    return relationships


def load_grounding_map(filename):
    gm_rows = load_csv(filename)
    g_map = {}
    for row in gm_rows:
        key = row[0]
        db_refs = {'TEXT': key}
        keys = [entry for entry in row[1::2] if entry != '']
        values = [entry for entry in row[2::2] if entry != '']
        db_refs.update(dict(zip(keys, values)))
        if len(db_refs.keys()) > 1:
            g_map[key] = db_refs
        else:
            g_map[key] = None
    return g_map


def load_gmap_reverse(filename):
    g_map = load_grounding_map(filename)
    g_map_reverse = defaultdict(set)
    for text, groundings in g_map.items():
        fplx_id = groundings.get('FPLX')
        if fplx_id:
            g_map_reverse[fplx_id].add(text)
    return g_map_reverse
