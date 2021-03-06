from collections import defaultdict, namedtuple, Counter
from indra.util import read_unicode_csv

def sorted_ctr(items):
    """Return a list of unique items in descending order by frequency."""
    ctr = Counter(items)
    return sorted([(k, v) for k, v in ctr.items()],
                  key=lambda x: x[1], reverse=True)


def load_curated_annotations(filename):
    def process_row(row):
        # Convert all nones in the boolean fields to 0s, and convert strings
        # '1' and '0' to ints so boolean operations will work
        def norm_int_field(entry):
            if entry == '' or entry == '0':
                return 0
            elif entry == '1':
                return 1
            else:
                raise ValueError('Entry %s is not an integer' % entry)

        if row[1] == '':
            is_family = None
        else:
            is_family = int(row[1])
        proc_row = ((int(row[0]), is_family, row[2], norm_int_field(row[3]),
                     norm_int_field(row[4]), norm_int_field(row[5])) +
                    tuple(row[6:]))
        return proc_row

    cols = ['id', 'is_family', 'entity_text', 'in_be_map', 'in_be',
            'false_positive', 'grounding', 'notes', 'doi', 'passage']
    CuratedRow = namedtuple('CuratedRow', cols)
    rows = []
    results = defaultdict(lambda: 0)
    for row in read_unicode_csv(filename, delimiter='\t', skiprows=1):
        proc_row = process_row(row)
        cr = CuratedRow(*proc_row)
        rows.append(cr)
    return rows


def load_curated_unmatched(filename):
    cols = ['text', 'count', 'is_protein', 'famplex_id']
    UnmatchedRow = namedtuple('UnmatchedRow', cols)
    rows = []
    results = defaultdict(lambda: 0)
    for row in read_unicode_csv(filename, delimiter='\t', skiprows=1):
        ur = UnmatchedRow(row[0], int(row[1]), int(row[2]), row[3])
        rows.append(ur)
    return rows


def score_results(cur_rows):
    results = {'total_curated': [0, []],
               'is_family': [0, []],
               'in_be_map': [0, []],
               'in_be': [0, []],
               'missing': [0, []],
               }
    def add_row(keyname, row):
        results[keyname][0] += 1
        results[keyname][1].append(row)

    for cr in cur_rows:
        if cr.is_family is None:
            continue
        add_row('total_curated', cr)
        if cr.is_family == 1:
            add_row('is_family', cr)
            if cr.in_be_map:
                add_row('in_be_map', cr)
            else:
                if cr.in_be:
                    add_row('in_be', cr)
                else:
                    add_row('missing', cr)
                    #print(cr)
    return results


def entity_freqs(results):
    freqs = {}
    for key, (count, rows) in results.items():
        freqs[key] = sorted_ctr([cr.entity_text for cr in rows])
    return freqs


def pct_str(num, den):
    pct = round(100 * num / den, 1)
    return '%s / %s (%s)' % (num, den, pct)

if __name__ == '__main__':
    cur_rows = load_curated_annotations('annotations_curated.tsv')
    results = score_results(cur_rows)
    print([(k, v[0]) for k, v in results.items()])
    print("-- Manual eval --")
    print("In grounding map: %s" %
          pct_str(results['in_be_map'][0], results['is_family'][0]))
    print("Has ID: %s" %
          pct_str(results['in_be_map'][0] + results['in_be'][0],
                  results['is_family'][0]))

    um_rows = load_curated_unmatched(
                    'automated_eval_unmatched_texts_curated.tsv')
    total_unmatched_entities = len(um_rows)
    total_unmatched_freq = sum([r.count for r in um_rows])
    unmatched_prot_rows = [r for r in um_rows if r.is_protein]
    unmatched_prot_entities = len(unmatched_prot_rows)
    unmatched_prot_freq = sum([r.count for r in unmatched_prot_rows])
    unmatched_has_id_rows = [r for r in unmatched_prot_rows if r.famplex_id]
    unmatched_has_id_entities = len(unmatched_has_id_rows)
    unmatched_has_id_freq = sum([r.count for r in unmatched_has_id_rows])

    num_matched = 1908
    total_unmatched_prot = num_matched + unmatched_prot_freq
    print("\n-- Automated eval --")
    print("In grounding map: %s" % pct_str(num_matched, total_unmatched_prot))
    print("Has ID: %s" % pct_str(num_matched + unmatched_has_id_freq,
                                 total_unmatched_prot))
