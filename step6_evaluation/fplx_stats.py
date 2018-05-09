"""This script generates overall statistics describing the FamPlex resource."""

import os
import json
import numpy
from collections import Counter
import matplotlib
matplotlib.use('svg')
import matplotlib.pyplot as plt
from util import *
from indra.literature.pubmed_client import get_ids
from indra.databases import uniprot_client
from indra.util import plot_formatting as pf

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
    """Return the children of a given FamPlex entry."""
    all_children = set()
    for rel in relations:
        if rel[2][0] == ns and rel[2][1] == id:
            all_children.add(rel[0])
            all_children |= get_children(*rel[0])
    return all_children


def get_level(entry):
    """Get the level at which a given FamPlex entry is in the hierarchy."""
    level = 1
    relevant_rels = [r for r in relations if r[2] == ('FPLX', entry)]
    for rel in relevant_rels:
        if rel[0][0] == 'FPLX':
            level = max(1 + get_level(rel[0][1]), level)
        else:
            level = max(2, level)
    return level


def num_child_concepts():
    """Print statistics of number of child concepts for each entry."""
    child_nums = []
    for entity in entities:
        children = get_children('FPLX', entity)
        child_nums.append(len(children))
        if len(children) > 30:
            print('%s: %s' % (entity, len(children)))

    pf.set_fig_params()
    plt.ion()
    plt.figure(figsize=(2.5, 2.5), dpi=300)
    plt.hist(child_nums, 60, color=pf.GREEN)
    plt.xlabel('Number of distinct children in FamPlex')
    plt.ylabel('Number of FamPlex entries')
    pf.format_axis(plt.gca())
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
        ent_level = get_level(entity)
        levels.append(ent_level)
        if ent_level == 4:
            print('%s has 4 levels' % entity)
    level_counter = Counter(levels)
    for num_level, count in sorted(level_counter.items(), key=lambda x: x[0]):
        print('Top-level FamPlex entries with %d levels: %d' %
              (num_level, count))


def num_citations():
    """Calculate the number of citations for each top-level entity and its
    children."""
    def get_pmids(search_terms):
        pmids = set()
        for lex in search_terms:
            search_term = '%s' % lex
            lex_pmids = set(get_ids(search_term, use_text_word=True,
                                    retmax=1000000))
            pmids |= lex_pmids
        return pmids

    cit_nums = {}
    all_pmids = set()
    for entity in entities:
        # Lexicalizations for the entry itself
        lexes = gmap_reverse.get(entity, set())
        lexes.add(entity.replace('_', '-'))
        pmids_entry = get_pmids(lexes)
        # Lexicalizations for the entries children
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
        pmids_inclusive = get_pmids(lexes)
        cit_nums[entity] = (len(pmids_entry), len(pmids_inclusive))
        all_pmids |= pmids_entry
        all_pmids |= pmids_inclusive
        print('%s: %d, %d' % (entity, len(pmids_entry), len(pmids_inclusive)))
    cit_num_entry = [c[0] for c in cit_nums.values()]
    cit_num_children = [c[1] for c in cit_nums.values()]
    for type, vals in zip(['entry', 'children'],
                          [cit_num_entry, cit_num_children]):
        print('Number of citations (%s): %.2f +/- %.2f, median: %d' %
              (type, numpy.average(vals), numpy.std(vals),
               numpy.median(vals)))
        print('All PMIDs found: %d' % len(all_pmids))
    with open('fplx_citations.json', 'w') as fh:
        json.dump(cit_nums, fh, indent=1)
    return cit_nums, all_pmids


def plot_cit_nums():
    """Plot the histogram of citation numbers for each FamPlex entry."""
    if not os.path.exists('fplx_citations.json'):
        num_citations()
    with open('fplx_citations.json', 'r') as fh:
        cit_nums = json.load(fh)
    pf.set_fig_params()
    plt.ion()
    plt.figure(figsize=(3.5, 2.5), dpi=300)
    cit_sort = sorted(cit_nums.values(), key=lambda x: x[1], reverse=True)
    y1 = [c[0] for c in cit_sort]
    y2 = [c[1]-c[0] for c in cit_sort]
    plt.bar(range(len(y1)), y1, color=pf.GREEN, label='Entry PMIDs')
    plt.bar(range(len(y2)), y2, color=pf.ORANGE, bottom=y1,
            label='Children PMIDs')
    plt.ylabel('Number of citations')
    plt.xlabel('FamPlex entries')
    plt.yscale('log')
    plt.xlim(0, len(y1))
    plt.legend()
    pf.format_axis(plt.gca())
    plt.savefig('fplx_pmids_bar.pdf')
    plt.show()


def num_genes_covered():
    """Calculate the overall number of specific genes that are covered by
    FamPlex."""
    genes_covered = set()
    for source, rel, target in relations:
        if source[0] in ('HGNC', 'UP'):
            genes_covered.add(source[1])
    print('Total number of unique genes covered: %d' % len(genes_covered))


if __name__ == '__main__':
    num_at_levels()
    num_child_concepts()
    top_level_depths()
    plot_cit_nums()
