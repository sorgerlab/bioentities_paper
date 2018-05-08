"""This script queries S3 to get the content type associated with
each PMID in the corpus and reports statistics.

NOTE: the S3 storage is not publicly accessible.
"""

import time
from multiprocessing import Pool
from collections import Counter
from indra.literature.s3_client import get_full_text

def get_ct(pmid):
    """Return the content type associated with the given PMID"""
    try:
        _, content_type = get_full_text(pmid, True)
    except Exception:
        content_type = None
    return content_type

if __name__ == '__main__':
    # Read the combined PMIDs from the text file
    with open('combined_pmids.txt', 'r') as fh:
        pmids = [l.strip() for l in fh.readlines()]

    # In many parallel threads, get the content type from S3
    # for each of the PMIDs
    p = Pool(64)
    ts = time.time()
    content_types = p.map(get_ct, pmids)
    te = time.time()
    print('%.2fs' % (te-ts))

    # Make a Counter for the content types and print it
    content_type_stats = Counter(content_types)
    print(content_type_stats)
