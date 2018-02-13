"""This script samples a subset of the PMIDs and reads
their abstracts with the TRIPS/DRUM reading system."""

import os
import random
from indra.literature import pubmed_client
from indra.sources.trips.drum_reader import DrumReader

def sample_trips_pmids(nabstracts):
    """Return a list of PMIDs to read sampled randomly."""
    # Load dict of PMIDs->statements
    with open('../step1_genes_pmids/combined_pmids.txt', 'r') as fh:
        pmids = [l.strip() for l in fh.readlines()]

    # Sort the PMIDs
    pmids = sorted(pmids)
    # Seed the random number generator for reproducibility
    random.seed(100)
    # Shuffle the PMIDs
    random.shuffle(pmids)
    # Number of abstracts to read with TRIPS
    pmids_to_read = pmids[:nabstracts]
    return pmids_to_read


def save_abstracts(pmids):
    """Download and save the abstracts for a list of PMIDs."""
    for pmid in pmids:
        abstract = pubmed_client.get_abstract(pmid)
        with open('trips_abstracts/%s.txt' % pmid, 'wb') as fh:
            fh.write(abstract.encode('utf-8'))


def read_abstract(pmid):
    """Read a given PMID's abstract with TRIPS/DRUM and save the result."""
    if os.path.exists('trips_abstracts/%s.ekb' % pmid):
        return
    with open('trips_abstracts/%s.txt' % pmid, 'rb') as fh:
        abstract = fh.read().decode('utf-8')
    dr = DrumReader(to_read=[abstract])
    try:
        dr.start()
    except SystemExit:
        print('====EXITED====')
        pass
    with open('trips_abstracts/%s.ekb' % pmid, 'wb') as fh:
        fh.write(dr.extractions[0].encode('utf-8'))

if __name__ == '__main__':
    nabstracts = 100
    pmids = sample_trips_pmids(nabstracts)
    save_abstracts(pmids)
    for pmid in pmids:
        print(pmid)
        read_abstract(pmid)
