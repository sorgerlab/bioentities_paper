"""This script calculates ranked lists of entities directly from the REACH
output to measure whether there is difference in constructing the same
statistics from event-extraction-based INDRA Statements."""

from indra.literature import s3_client

def get_entities(reach_json):
    """Return raw entities and event-extracted entities from REACH output.

    Parameters
    ----------
    reach_json : str
        The output of REACH for a given PMID.

    Returns
    -------
    raw_entities, event_entities : tuple
        raw_entities: the full list of entities found in text irrespective
            of their appearence in events
        event_entities: the list of entities that are arguments of extracted
            events
    """
    raw_entities = reach_json['entities']['frames']
    raw_entity_map = {ent['frame-id']: ent for ent in raw_entities}
    events = reach_json['events']['frames']
    event_entities = []
    for event in events:
        args = event.get('arguments', [])
        for arg in args:
            if arg['argument-type'] == 'entity':
                entity = raw_entity_map[arg['arg']]
                event_entities.append(entity)
    print('%d/%d: %d/%d entities from %s' %
          (i, len(pmids), len(raw_entities), len(event_entities), pmid))
    return raw_entities, event_entities


def update_entity_dict(entity_dict, entities):
    """Update a summary dict of entities and groundings given new entities."""
    for entity in entities:
        txt = entity.get('text')
        if not txt:
            continue
        if txt not in entity_dict:
            entity_dict[txt] = {}
        xrefs = entity['xrefs']
        if not xrefs:
            key = (None, None)
            try:
                entity_dict[txt][key] += 1
            except KeyError:
                entity_dict[txt][key] = 1
        else:
            for xref in xrefs:
                if xref['namespace'] == 'uaz':
                    key = ('uaz', 'UAZ')
                else:
                    key = (xref['namespace'], xref['id'])
                try:
                    entity_dict[txt][key] += 1
                except KeyError:
                    entity_dict[txt][key] = 1


if __name__ == '__main__':
    raw_entities_all = {}  # All raw entities and groundings
    event_entities_all = {}  # All event argument entities and groundings

    # Read in list of PMIDs
    with open('../step1_genes_pmids/combined_pmids.txt', 'r') as fh:
        pmids = [l.strip() for l in fh.readlines()]

    # Download and prcess output for each PMID and update dicts of entities
    for i, pmid in enumerate(pmids):
        # Get REACH JSON string from S3 for given PMID
        reach_json = s3_client.get_reader_output('reach_no_famplex', pmid)
        if not reach_json:
            continue
        # Extract the entities from the JSON and update dicts
        raw_entities, event_entities = get_entities(reach_json)
        update_entity_dict(raw_entities_all, raw_entities)
        update_entity_dict(event_entities_all, event_entities)
    # Rank list of entities for each case according to number of times
    # extracted
    raw_ranked = sorted(raw_entities_all.items(),
                        key=lambda x: sum([v for k, v in x[1].items()]),
                        reverse=True)
    event_ranked = sorted(event_entities_all.items(),
                          key=lambda x: sum([v for k, v in x[1].items()]),
                          reverse=True)
