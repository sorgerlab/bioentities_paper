"""This script generates overall statistics describing the FamPlex resource."""

import numpy
from collections import Counter
import matplotlib.pyplot as plt
from util import *
from indra.literature.pubmed_client import get_ids
from indra.databases import uniprot_client


entities = load_entity_list('../../famplex/entities.csv')
relations = load_relationships('../../famplex/relations.csv')
gmap = load_grounding_map('../../famplex/grounding_map.csv')
gmap_reverse = load_gmap_reverse('../../famplex/grounding_map.csv')

def num_at_levels():
    """Print number of FPLX entries at each level of relations."""
    top_level = 0
    mid_level = 0
    bottom_level = 0
    no_rel = 0
    for entity in entities:
        rel_par = [r for r in relations if
                   r[0][0] == 'FPLX' and r[0][1] == entity]
        rel_child = [r for r in relations if
                     r[2][0] == 'FPLX' and r[2][1] == entity]
        if rel_child and rel_par:
            mid_level += 1
        if rel_child and not rel_par:
            top_level += 1
        if not rel_child and rel_par:
            bottom_level += 1
        if not rel_child and not rel_par:
            no_rel += 1
    print('Total entities: %d, top: %d, mid: %d, bottom: %d, no relations: %d'
          % (len(entities), top_level, mid_level, bottom_level, no_rel))


def get_children(ns, id):
    all_children = set()
    for rel in relations:
        if rel[2][0] == ns and rel[2][1] == id:
            all_children.add(rel[0])
            all_children |= get_children(*rel[0])
    return all_children


def get_level(entry):
    level = 1
    relevant_rels = [r for r in relations if r[2] == ('FPLX', entry)]
    for rel in relevant_rels:
        if rel[0][0] == 'FPLX':
            level = max(1 + get_level(rel[0][1]), level)
        else:
            level = max(2, level)
    return level


def num_child_concepts():
    """Statistics of number of child concepts for each entry."""
    child_nums = []
    for entity in entities:
        children = get_children('FPLX', entity)
        child_nums.append(len(children))

    plt.ion()
    plt.figure()
    plt.hist(child_nums, 50, color='gray')
    plt.xlabel('Number of distinct children in FamPlex')
    plt.ylabel('Number of FamPlex entries')
    plt.savefig('fplx_children_hist.pdf')
    plt.show()

    print('Number of children: %.2f +/- %.2f, median: %d' %
          (numpy.average(child_nums), numpy.std(child_nums),
           numpy.median(child_nums)))


def top_level_depths():
    """Statistics on the depth of concepts below each top-level entry."""
    levels = []
    for entity in entities:
        rel_par = [r for r in relations if
                   r[0][0] == 'FPLX' and r[0][1] == entity]
        if rel_par:
            continue
        levels.append(get_level(entity))
    level_counter = Counter(levels)
    for num_level, count in sorted(level_counter.items(), key=lambda x: x[0]):
        print('Top-level FamPlex entries with %d levels: %d' %
              (num_level, count))


def num_citations():
    cit_nums = {}
    all_pmids = set()
    for entity in entities:
        rel_par = [r for r in relations if
                   r[0][0] == 'FPLX' and r[0][1] == entity]
        if rel_par:
            continue
        lexes = gmap_reverse.get(entity, set())
        lexes.add(entity.replace('_', ' '))
        children = get_children('FPLX', entity)
        for child_ns, child_name in children:
            if child_ns == 'FPLX':
                lexes |= gmap_reverse.get(child_name)
            elif child_ns == 'HGNC':
                lexes.add(child_name)
            elif child_ns == 'UP':
                gene_name = uniprot_client.get_gene_name(child_name)
                if gene_name:
                    lexes.add(gene_name)
        pmids = set()
        for lex in lexes:
            search_term = '%s[Text Word]' % lex
            lex_pmids = get_ids(search_term, retmax=1000000)
            pmids |= set(lex_pmids)
        cit_nums[entity] = len(pmids)
        all_pmids |= pmids
        print('%s: %d' % (entity, len(pmids)))
    cit_num_vals = list(cit_nums.values())
    print('Number of citations: %.2f +/- %.2f, median: %d' %
          (numpy.average(cit_num_vals), numpy.std(cit_num_vals),
           numpy.median(cit_num_vals)))
    print('All PMIDs found: %d' % len(all_pmids))

    plt.ion()
    plt.figure()
    plt.hist(cit_num_vals, 50, color='gray')
    plt.xlabel('Number of citations for FamPlex entry and its children')
    plt.ylabel('Number of FamPlex entries')
    plt.savefig('fplx_pmids_hist.pdf')
    plt.show()

    return cit_nums, all_pmids


if __name__ == '__main__':
    num_at_levels()
    num_child_concepts()
    top_level_depths()
