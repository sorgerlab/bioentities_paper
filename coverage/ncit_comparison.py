"""This script looks up all the mappings from FPLX to NCIT and
checks if the NCIT entry has any child concepts."""

import os
import csv
import requests
import zipfile
from collections import defaultdict
from util import load_equivalences

flat_file = 'https://evs.nci.nih.gov/ftp1/NCI_Thesaurus/Thesaurus_18.02d.FLAT.zip'

if not os.path.exists('Thesaurus.txt'):
    res = requests.get(flat_file)
    with open('Thesaurus_18.02d.FLAT.zip', 'wb') as fh:
        fh.write(res.content)
    with zipfile.ZipFile('Thesaurus_18.02d.FLAT.zip', 'r') as fh:
        fh.extractall('.')

children_dict = defaultdict(set)

with open('Thesaurus.txt', 'r') as fh:
    reader = csv.reader(fh, delimiter='\t')
    for row in reader:
        child = row[0]
        parents = row[2].strip().split('|')
        for parent in parents:
            children_dict[parent].add(child)

equivalences = load_equivalences('../../famplex/equivalences.csv')
ncit_ids = [e[1] for e in equivalences if e[0] == 'NCIT']

has_children = 0
num_ids = len(ncit_ids)
for ncit_id in ncit_ids:
    children = children_dict[ncit_id]
    if children:
        has_children += 1
print('Number of NCIT mappings: %d, has children: %d, %.2f%%' % 
      (num_ids, has_children, 100*has_children/num_ids))
