import requests
import logging
from indra.databases import hgnc_client

logger = logging.getLogger('reactome')
#curl -X GET --header 'Accept: application/json'

def query_id(up_id):
    react_search_url = 'http://www.reactome.org/ContentService/search/query'
    params = {'query': up_id, 'cluster': 'true', 'species':'Homo sapiens'}
    headers = {'Accept': 'application/json'}
    res = requests.get(react_search_url, headers=headers, params=params)
    if not res.status_code == 200:
        return None
    json = res.json()
    results = json.get('results')
    if not results:
        print('No results for %s' % up_id)
        return None
    logger.info('len(results) = %d' % len(results))
    stable_ids = []
    for result in results:
        entries = result.get('entries')
        for entry in entries:
            stable_id = entry.get('stId')
            if not stable_id:
                continue
            name = entry.get('name')
            stable_ids.append(stable_id)
            print('%s: %s' % (name, stable_id))

    return stable_ids


def get_parents(stable_id):
    react_data_url = 'http://www.reactome.org/ContentService/data/entity/' + \
                     stable_id + '/componentOf'
    headers = {'Accept': 'application/json'}
    res = requests.get(react_data_url, headers=headers)
    if not res.status_code == 200:
        return []
    json = res.json()
    names = []
    stable_ids = []
    schema_classes = []
    for parent_group in json:
        if not parent_group.get('type') in \
                            ['hasComponent', 'hasMember', 'hasCandidate']:
            continue
        names += parent_group.get('names')
        stable_ids += parent_group.get('stIds')
        schema_classes += parent_group.get('schemaClasses')
    parents_at_this_level = list(zip(names, stable_ids, schema_classes))
    parents_at_next_level_up = []
    for p_name, p_id, sc in parents_at_this_level:
        parents_at_next_level_up += get_parents(p_id)
    return parents_at_this_level + parents_at_next_level_up

def get_all_parents(up_id):
    linked_stable_ids = query_id(up_id)
    parents = []
    for ls_id in linked_stable_ids:
        parents += get_parents(ls_id)
    sets = [tup for tup in parents if tup[2] != 'Complex']
    complexes = [tup for tup in parents if tup[2] == 'Complex']
    return sets, complexes

if __name__ == '__main__':
    import csv
    genes = []
    with open('../../bioentities/relations.csv', 'rt') as fh:
        csvreader = csv.reader(fh, delimiter=',', quotechar='"')
        for row in csvreader:
            if row[0] == 'HGNC':
                genes.append(row[1])
    with open('bioentities_genes.csv', 'wt') as f:
        for gene in genes:
            f.write('%s\n' % gene)

    gene_ids = []
    for hgnc_sym in genes:
        hgnc_id = hgnc_client.get_hgnc_id(hgnc_sym)
        up_id = hgnc_client.get_uniprot_id(hgnc_id)
        gene_ids.append((hgnc_sym, up_id))
    gene_ids = list(set(gene_ids))

    # Iterate over all genes and get all the sets (families)
    parents = []
    for hgnc_sym, up_id in gene_ids:
        parents = get_all_parents(up_id)

    """
    # For a given gene, get parents in Bioentities, along with members
    # For a given set in Reactome, get all members
    # Compare all pairs of sets, Bioentities vs. Reactome
    # Determine which Reactome set, if any, has the closest overlap with
    # the Bioentities set.
    # Could do something similar with NCIT??
    """

