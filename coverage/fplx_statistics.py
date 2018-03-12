"""This script generates overall statistics describing the FamPlex resource."""

import numpy
import matplotlib.pyplot as plt
from util import *

entities = load_entity_list('../../famplex/entities.csv')
relations = load_relationships('../../famplex/relations.csv')

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

def num_child_concepts():
    """Statistics of number of child concepts for each entry."""
    child_nums = []
    for entity in entities:
        children = get_children('FPLX', entity)
        child_nums.append(len(children))

    plt.figure()
    plt.hist(child_nums, 50, color='gray')
    plt.xlabel('Number of distinct children in FamPlex')
    plt.ylabel('Number of FamPlex entries')
    plt.savefig('fplx_children_hist.pdf')
    plt.show()

    print('Number of children: %.2f +/- %.2f, median: %d' %
          (numpy.average(child_nums), numpy.std(child_nums),
           numpy.median(child_nums)))

if __name__ == '__main__':
    num_at_levels()
    num_child_concepts()

