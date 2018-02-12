"""This script processes a set of EKB extractions by TRIPS and builds
up dictionaries to represent the statistics of grounding raw text
entities to database indentifiers."""
import os
import sys
import glob
import pickle
import xml.etree.ElementTree as ET


def get_terms_with_grounding(tree):
    term_grounding = {}
    terms = tree.findall('TERM')
    for term in terms:
        for drum_term in term.findall('drum-terms/drum-term'):
            mname = drum_term.attrib.get('matched-name')
            score = drum_term.attrib.get('match-score')
            name = drum_term.attrib.get('name')
            grounding = drum_term.attrib.get('dbid')
            dbname, dbid = grounding.split(':')
            key = (dbname, dbid, name, score)
            if mname not in term_grounding:
                term_grounding[mname] = {}
            try:
                term_grounding[mname][key] += 1
            except KeyError:
                term_grounding[mname][key] = 1
    return term_grounding


def merge_dicts(d1, d2):
    for txt, key_counts in d2.items():
        if txt not in d1:
            d1[txt] = d2[txt]
        else:
            for key, count in key_counts.items():
                if key not in d1[txt]:
                    d1[txt][key] = d2[txt][key]
                else:
                    d1[txt][key] += d2[txt][key]


if __name__ == '__main__':
    # The folder to search for *.ekb files
    folder = sys.argv[1]
    fnames = glob.glob(os.path.join(folder, '*.ekb'))
    all_term_grounding = {}
    # Extract grounding stats from each file and merge
    for fname in fnames:
        with open(fname, 'rt') as fh:
            tree = ET.parse(fh)
            term_grounding = get_terms_with_grounding(tree)
            merge_dicts(all_term_grounding, term_grounding)
    with open('trips_entities.pkl', 'wb') as fh:
        pickle.dump(all_term_grounding, fh)
