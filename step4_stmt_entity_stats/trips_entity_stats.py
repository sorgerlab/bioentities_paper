"""This script processes a set of EKB extractions by TRIPS and builds
up dictionaries to represent the statistics of grounding raw text
entities to database indentifiers."""
import os
import sys
import csv
import glob
import numpy
import pickle
import random
import xml.etree.ElementTree as ET


def get_terms_with_grounding(tree):
    term_grounding = []
    terms = tree.findall('TERM')
    for term in terms:
        drum_terms = term.findall('drum-terms/drum-term')
        mnames = set(term.get('matched-name') for term in drum_terms)
        uttnum = term.attrib['uttnum']
        sentence = get_evidence_sentence(tree, uttnum)
        if not mnames:
            name = term.find('name')
            if name is None:
                continue
            term_grounding.append((term.find('name').text, sentence, []))
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
        term_grounding.append((mname, sentence, values))
    return term_grounding


def get_evidence_sentence(tree, uttnum):
    sentence_tag = tree.find("input/sentences/sentence[@id='%s']" % uttnum)
    if len(sentence_tag.text) > 200:
        return sentence_tag.text[:200]
    return sentence_tag.text


def sample_spreadsheet(groundings, nsample, out_file):
    sampled = [random.choice(groundings) for _ in range(nsample)]
    with open(out_file, 'w') as fh:
        fh.write('Text,Family,Evidence,TopCorrect,AnyCorrect,Groundings\n')
        for name, sentence, grounding in sampled:
            match_strs = []
            for match in grounding:
                match_str = '%s/%s/%s/%s' % match
                match_strs.append(match_str)
            if grounding:
                fh.write('%s,,"%s",,,' % (name, sentence))
            else:
                fh.write('%s,,"%s",0,0,' % (name, sentence))
            fh.write(','.join(match_strs))
            fh.write('\n')


def analyze_curated_spreadsheet(fname):
    ste = lambda k, n: numpy.sqrt(k/n * (1-k/n)/n)
    ncurated = 0
    nfam = 0
    nfam_top_correct = 0
    nfam_any_correct = 0
    top_ncit = 0
    top_be = 0
    with open(fname, 'r') as fh:
        reader = csv.reader(fh)
        for row in reader:
            is_family, top_correct, any_correct = row[1], row[3], row[4]
            if is_family in ('0', '1'):
                ncurated += 1
            if is_family == '1':
                nfam += 1
            if top_correct == '1':
                nfam_top_correct += 1
                if row[5].startswith('NCIT'):
                    top_ncit += 1
                if row[5].startswith('BE'):
                    top_be += 1
            if any_correct == '1':
                nfam_any_correct += 1
    print('Curated: %d' % ncurated)
    print('Top correct: mean=%.1f, ste=%.1f' %
          (100.0 * nfam_top_correct / nfam,
           100.0 * ste(nfam_top_correct, nfam)))
    print('Any correct: mean=%.1f, ste=%.1f' %
          (100.0 * nfam_any_correct / nfam,
           100.0 * ste(nfam_any_correct, nfam)))
    print('Top correct that is NCIT: %.1f' %
          (100.0 * top_ncit / nfam_top_correct))
    print('Top correct that is BE: %.1f' %
          (100.0 * top_be / nfam_top_correct))


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
    #sample_spreadsheet(all_term_grounding, 1000, 'trips_entities_with_fplx.csv')
    sample_spreadsheet(all_term_grounding, 1000, 'trips_entities_no_fplx.csv')
