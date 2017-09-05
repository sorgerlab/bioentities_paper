import os
import sys
import csv
from indra.databases import get_identifiers_url

def get_link(db_name, db_ref):
    url = get_identifiers_url(db_name, db_ref)
    link = '<a href="%s" target="_blank">%s:%s</a>' % (url, db_name, db_ref)
    return link

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: %s [training|test]' % os.path.basename(__file__))
        sys.exit()
    mode = sys.argv[1]
    if mode not in ('training', 'test'):
        print('Usage: %s [training|test]' % os.path.basename(__file__))
        sys.exit()

    sample_file = 'training_agents_sample_curated.csv' if mode == 'training' \
        else 'test_agents_with_be_sample.csv'
    out_file = 'training_agents_sample_links.html' if mode == 'training' \
            else 'test_agents_with_be_sample.html'

    with open(sample_file, 'r') as fh1, open(out_file, 'w') as fh:
        reader = csv.reader(fh1)
        next(reader)
        fh.write('<body><table border=1>')
        for row in reader:
            fh.write('<tr><td>%s</td>' % row[2])
            link = 'https://www.ncbi.nlm.nih.gov/pubmed/%s' % row[9]
            fh.write('<td><a href="%s" target="_blank">%s</a></td>' %
                     (link, row[9]))
            for i in range(3):
                db_name = row[2*i+3]
                if db_name == 'IPR':
                    db_name = 'IP'
                db_id = row[2*i+4]
                if not db_name:
                    fh.write('<td>&nbsp;</td>')
                    continue
                link = get_link(db_name, db_id)
                fh.write('<td>%s</td>' % link)
            fh.write('</tr>')
        fh.write('</table></body>')
