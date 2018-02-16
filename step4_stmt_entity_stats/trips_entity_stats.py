"""This script processes a set of EKB extractions by TRIPS and builds
up dictionaries to represent the statistics of grounding raw text
entities to database indentifiers."""
import os
import sys
import glob
import pickle
import random
import xml.etree.ElementTree as ET


def get_terms_with_grounding(tree):
    term_grounding = []
    terms = tree.findall('TERM')
    for term in terms:
        drum_terms = term.findall('drum-terms/drum-term')
        mnames = set(term.get('matched-name') for term in drum_terms)
        if not mnames:
            name = term.find('name')
            if name is None:
                continue
            term_grounding.append((term.find('name').text, []))
            continue
        # Matched names, if there are more than one, are always just
        # trivially different e.g. PAI-1 vs PAI1 so we can just take
        # the first one here
        mname = list(mnames)[0]
        values = []
        for drum_term in drum_terms:
            grounding = drum_term.attrib.get('dbid')
            dbname, dbid = grounding.split(':')
            score = drum_term.attrib.get('match-score')
            name = drum_term.attrib.get('name')
            value = (dbname, dbid, name, score)
            values.append(value)
        term_grounding.append((mname, values))
    return term_grounding


def sample_spreadsheet(groundings, nsample, out_file):
    sampled = [random.choice(groundings) for _ in range(nsample)]
    with open(out_file, 'w') as fh:
        fh.write('Text,Family,TopCorrect,AnyCorrect,Groundings\n')
        for name, grounding in groundings:
            match_strs = []
            for match in grounding:
                match_str = '%s/%s/%s/%s' % match
                match_strs.append(match_str)
            fh.write('%s,,,,' % name)
            fh.write(','.join(match_strs))
            fh.write('\n')


if __name__ == '__main__':
    # The folder to search for *.ekb files
    folder = sys.argv[1]
    fnames = sorted(list(glob.glob(os.path.join(folder, '*.ekb'))))
    print('Found %s EKB files to process' % len(fnames))
    all_term_grounding = []
    # Extract grounding stats from each file and merge
    for fname in fnames:
        with open(fname, 'rt') as fh:
            tree = ET.parse(fh)
            term_grounding = get_terms_with_grounding(tree)
            all_term_grounding += term_grounding
    random.seed(1)
    sample_spreadsheet(all_term_grounding, 1000, 'trips_entities_with_be.csv')
