import os
from collections import namedtuple
from collections import Counter, defaultdict
import csv
from indra.preassembler.grounding_mapper import load_grounding_map
from xml.etree import ElementTree as ET
import random
from indra.sources import reach
from indra.databases import hgnc_client, uniprot_client
from indra.util import write_unicode_csv

def norm_text(text):
    """Normalize text for case-insensitive matches."""
    normed = text.replace('‐', '-')
    normed = normed.upper()
    return normed


def sorted_ctr(items):
    """Return a list of unique items in descending order by frequency."""
    ctr = Counter(items)
    return sorted([(k, v) for k, v in ctr.items()],
                  key=lambda x: x[1], reverse=True)


def format_id(ns, id):
    """Format a namespace/ID pair for display and curation."""
    label = '%s:%s' % (ns, id)
    label = label.replace(' ', '_')
    if ns == 'Uniprot':
        url = 'http://identifiers.org/uniprot/%s' % id
    elif ns == 'NCBI gene':
        url = 'http://identifiers.org/ncbigene/%s' % id
    else:
        url = None
    return (label, url)


BioIdRow = namedtuple('BioIdRow',
                      ['thomas_article', 'doi', 'don_article', 'figure',
                       'annot_id', 'paper_id', 'first_left', 'last_right',
                       'length', 'byte_length', 'left_alphanum', 'text',
                       'right_alphanum', 'obj', 'overlap', 'identical_span',
                       'overlap_label_count'])


def process_annotation_file(filename):
    """Process the annotations file into a list of BioIDRow objects."""
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
                obj = tuple([x.strip()
                             for x in obj.strip().split(':', maxsplit=1)])
                # Unfortunately the Cellosaurus annotations don't have an
                # explicit namespace, so we add one for completeness
                if len(obj) < 2:
                    assert(obj[0].startswith('CVCL'))
                    obj = ('Cellosaurus', obj[0])
                objs.append(obj)
            mapping = (text, tuple(objs))
            mappings.append(mapping)
            row[OBJ_COL] = mapping
            bioidrow = BioIdRow(*row)
            bioidrows.append(bioidrow)
    return bioidrows


def get_annotation_text(ann):
    # Get the caption text for this annotation
    xml_path = os.path.join('BioIDtraining_2/caption_bioc',
                            ann.don_article + '.xml')
    with open(xml_path, 'rt') as f:
        tree = ET.fromstring(f.read())
    docs = tree.findall('./document')
    for doc in docs:
        pmc_id = doc.find('./infon[@key="pmc_id"]').text
        figure_id = doc.find('./infon[@key="sourcedata_figure_dir"]').text
        if pmc_id == ann.don_article and figure_id == ann.figure:
            passage = doc.find('./passage/text').text
    return passage


def write_curation_tsv(annotations, be_strings, output_file):
    rows = []
    for ix, ann in enumerate(annotations):
        ix += 1
        in_bioent = True if norm_text(ann.text) in be_strings_norm \
                         else False
        passage = get_annotation_text(ann)
        fl = int(ann.first_left)
        lr = int(ann.last_right)
        hl_passage = ('%s[[[%s]]]%s' %
                      (passage[0:fl], passage[fl:lr], passage[lr:]))
        identifier_links = []
        for ns, id in ann.obj[1]:
            label, url = format_id(ns, id)
            if url is None:
                identifier_links.append(label)
            else:
                identifier_links.append(url)
        identifier_links = ' | '.join(identifier_links)
        doi = 'https://dx.doi.org/%s' % ann.doi
        row = (ix, ann.text, str(in_bioent), identifier_links, doi, hl_passage)
        rows.append(row)
    write_unicode_csv(output_file, rows, delimiter='\t')


def write_curation_html(annotations, be_strings, output_file):
    # Header for the HTML file
    html = """
        <html>
        <head>
        <meta charset="UTF-8" />
        <title>Biocreative VI annotations for curation</title></head>
        <body>
        <table border=1>
        <tr>
        <th>ID</th><th>Text</th><th>In BE?</th><th>Mappings</th>
        <th>Caption text</th>
        </tr>
    """
    for ix, ann in enumerate(annotations[0:100]):
        ix += 1
        in_bioent = 'False'
        # Check if this text is matched in the grounding map
        if norm_text(ann.text) in be_strings_norm:
            in_bioent = 'True'
        # Get the annotation text
        passage = get_annotation_text(ann)
        fl = int(ann.first_left)
        lr = int(ann.last_right)
        hl_passage = ('%s<mark>%s</mark>%s' %
                      (passage[0:fl], passage[fl:lr], passage[lr:]))
        # Add hyperlinks to grounding
        identifier_links = '<ul>\n'
        for ns, id in ann.obj[1]:
            label, url = format_id(ns, id)
            if url is None:
                entry = '<li>%s</li>\n' % label
            else:
                entry = '<li><a href="%s">%s</a></li>\n' % (url, label)
            identifier_links += entry
        identifier_links += '</ul>\n'
        # Add the row to the table
        tr = ('<tr><td>%d</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td>\n' %
              (ix, ann.text, in_bioent, identifier_links, hl_passage))
        html += tr
    html += "</body></html>"

    with open(output_file, 'wt', encoding='utf8') as f:
        f.write(html)


def filter_to_namespaces(annotations, namespaces):
    """Filter annotations to those grounded to one of the given namespaces."""
    filt = []
    for ann in annotations:
        gene_protein = False
        for mapping in ann.obj[1]:
            if mapping[0] in namespaces:
                gene_protein = True
        if gene_protein:
            filt.append(ann)
    return filt


def process_gene_prefixes(filename):
    prefixes = set()
    with open(filename, 'rt') as f:
        csvreader = csv.reader(f)
        for row in csvreader:
            if row[1] == 'experimental context':
                text = row[0]
                text = text.replace('{Gene name}', '')
                text = text.replace('-', '')
                text = text.strip()
                prefixes.add(text)
    return list(prefixes)


def filter_human(annotations):
    filt = []
    for ann in annotations:
        non_human = False
        for ns, id in ann.obj[1]:
            if ns == 'NCBI gene':
                hgnc_id = hgnc_client.get_hgnc_from_entrez(id)
                if not hgnc_id:
                    non_human = True
            elif ns == 'Uniprot':
                if not uniprot_client.is_human(id):
                    non_human = True
        if not non_human:
            filt.append(ann)
    return filt


if __name__ == '__main__':
    # Load and normalize the grounding map
    bioentities_path = os.environ['BIOENTITIES_HOME']
    be_grounding_map = load_grounding_map(os.path.join(bioentities_path,
                                         'grounding_map.csv'))
    be_strings_norm = [norm_text(text) for text in be_grounding_map]
    be_gene_prefix_file = os.path.join(bioentities_path, 'gene_prefixes.csv')
    norm_prefixes = process_gene_prefixes(be_gene_prefix_file)

    # Get the annotation data
    annotations_file = 'BioIDtraining_2/annotations.csv'
    annotations = process_annotation_file(annotations_file)

    # Filter the annotations to those identified as genes/proteins
    gp_grounded = filter_to_namespaces(annotations, ('NCBI gene', 'Uniprot'))
    # Filter to human
    gp_human = filter_human(gp_grounded)
    print("Grounded:", len(gp_grounded))
    print("Grounded human:", len(gp_human))
    gp_ungrounded = filter_to_namespaces(annotations, ('protein', 'gene'))
    gp_ungrounded_non_exp = []
    gp_ungrounded_exp = []
    for ann in gp_ungrounded:
        if ann.text in norm_prefixes:
            gp_ungrounded_exp.append(ann)
        else:
            gp_ungrounded_non_exp.append(ann)

    # Optional: filter to only those that are either ungrounded or are
    # grounded to multiple IDs
    #gp_fam = [ann for ann in gp_grounded if len(ann.obj[1]) > 1]
    anns_for_curation = gp_human + gp_ungrounded_non_exp

    # Shuffle the annotations
    random.seed(1)
    random.shuffle(anns_for_curation)

    #write_curation_html(anns_for_curation, be_strings_norm,
    #                    'ann_sample_test.html')
    write_curation_tsv(anns_for_curation, be_strings_norm,
                      'ann_sample_test.tsv')

