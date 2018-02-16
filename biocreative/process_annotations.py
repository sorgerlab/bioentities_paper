import os
from collections import namedtuple
from collections import Counter, defaultdict
import csv
from indra.preassembler.grounding_mapper import default_grounding_map
from xml.etree import ElementTree as ET
import random
from indra.sources import reach


random.seed(1)

bioentities_path = os.environ['BIOENTITIES_HOME']
annotations_file = 'BioIDtraining_2/annotations.csv'

def sorted_ctr(items):
    ctr = Counter(items)
    return sorted([(k, v) for k, v in ctr.items()],
                  key=lambda x: x[1], reverse=True)

BioIdRow = namedtuple('BioIdRow',
                      ['thomas_article', 'doi', 'don_article', 'figure',
                       'annot_id', 'paper_id', 'first_left', 'last_right',
                       'length', 'byte_length', 'left_alphanum', 'text',
                       'right_alphanum', 'obj', 'overlap', 'identical_span',
                       'overlap_label_count'])

# Column indices in the annotations.csv file
TEXT_COL = 11
OBJ_COL = 13

mappings = []
bioidrows = []

# Process the annotations CSV file
with open(annotations_file, 'rt') as f:
    csv_reader = csv.reader(f, delimiter=',')
    # Skip the header row
    next(csv_reader)
    for row in csv_reader:
        text = row[TEXT_COL]
        # Split the obj field by the pipe character for multiple entries
        objs = []
        for obj in row[OBJ_COL].split('|'):
            # Split each mapping into namespace and ID on the : character
            obj = tuple([x.strip() for x in obj.strip().split(':', maxsplit=1)])
            # Unfortunately the Cellosaurus annotations don't have an explicit
            # namespace, so we add one for completeness
            if len(obj) < 2:
                assert(obj[0].startswith('CVCL'))
                obj = ('Cellosaurus', obj[0])
            objs.append(obj)
        mapping = (text, tuple(objs))
        mappings.append(mapping)
        row[OBJ_COL] = mapping
        bioidrow = BioIdRow(*row)
        bioidrows.append(bioidrow)

namespace_list = [mapping[0] for text, mapping_list in mappings
                             for mapping in mapping_list]
namespaces_ctr = sorted_ctr(namespace_list)

# Find exact matches between text entries and grounding map
matches = defaultdict(list)
match_count = 0
for text, mapping_list in mappings:
    if text in default_grounding_map:
        match_count += 1
        matches[text].append(frozenset(mapping_list))

for text, mapping_list in matches.items():
    matches[text] = Counter(mapping_list)

# Sort text instances mapped to namespace 'protein' and 'gene'
gp_texts = [text for text, mapping_list in mappings
                 for mapping in mapping_list
                 if mapping[0] in ('protein', 'gene', 'NCBI gene', 'Uniprot')]
gp_texts_ctr = sorted_ctr(gp_texts)

gp_rows = [row for row in bioidrows
               for mapping in row.obj[1]
               if mapping[0] in ('protein', 'gene', 'NCBI gene', 'Uniprot')]
random.shuffle(gp_rows)


html = """
<html>
<head>
<meta charset="UTF-8" />
<title>Biocreative VI annotations for curation</title></head>
<body>
<table border=1>
<tr>
<th>Text</th><th>In BE?</th><th>Mappings</th><th>Caption text</th>
</tr>
"""
def format_id(ns, id):
    label = '%s:%s' % (ns, id)
    label = label.replace(' ', '_')
    if ns == 'Uniprot':
        url = 'http://identifiers.org/uniprot/%s' % id
    elif ns == 'NCBI gene':
        url = 'http://identifiers.org/ncbigene/%s' % id
    else:
        url = None
    return (label, url)

for ann in gp_rows[0:100]:
    in_bioent = '<mark>False</mark>'
    if ann.text in default_grounding_map:
        in_bioent = 'True'
    xml_path = os.path.join('BioIDtraining_2/caption_bioc',
                            ann.don_article + '.xml')
    with open(xml_path, 'rt') as f:
        tree = ET.fromstring(f.read())
    docs = tree.findall('./document')
    for doc in docs:
        pmc_id = doc.find('./infon[@key="pmc_id"]').text
        figure_id = doc.find('./infon[@key="sourcedata_figure_dir"]').text
        if pmc_id == ann.don_article and figure_id == ann.figure:
            #print("--------------------------------")
            #print(pmc_id)
            #print(figure_id)
            passage = doc.find('./passage/text').text
            fl = int(ann.first_left)
            lr = int(ann.last_right)
            hl_passage = ('%s<mark>%s</mark>%s' %
                          (passage[0:fl], passage[fl:lr], passage[lr:]))
            #print("Mapping: %s" % str(ann.obj))
            #print()
            #print(hl_passage)
            #print("Processing with REACH")
            #rp = reach.process_text(passage)
    identifier_links = '<ul>\n'
    for ns, id in ann.obj[1]:
        label, url = format_id(ns, id)
        if url is None:
            entry = '<li>%s</li>\n' % label
        else:
            entry = '<li><a href="%s">%s</a></li>\n' % (url, label)
        identifier_links += entry
    identifier_links += '</ul>\n'
    tr = ('<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td>\n' %
          (ann.text, in_bioent, identifier_links, hl_passage))
    html += tr
html += "</body></html>"

with open('ann_sample.html', 'wt', encoding='utf8') as f:
    f.write(html)

#print("Matches:", match_count)
#print(gp_texts)
#print(namespaces)
#print(mappings[0:20])
