__doc__ = """
This script is actually independent of settings.
"""

import os

out_folder = 'clinvar_split'
os.mkdir(out_folder)

filename = 'variant_summary.txt.gz'

import gzip, csv, re
from collections import defaultdict

variants = defaultdict(list)

with gzip.open(filename) as fh:
    split = lambda line: line.decode().strip().replace('#', '').split('\t')
    headers = split(next(fh))
    parse = lambda line: dict(zip(headers, split(line)))
    for line in fh:
        data = parse(line)
        # print(data)
        # extract mutation (to residue may be jibberish)
        pattern = r'\(p\.(?P<from_residue>\w{3})(?P<residue_index>\d+)(?P<to_residue>.*)\)'
        mutation_match = re.search(pattern, data['Name'])
        if not mutation_match:
            continue
        mutation = mutation_match.groupdict()
        if mutation['to_residue'] == '=':
            continue
        mutation['residue_index'] = int(mutation['residue_index'])
        mapping = {'AlleleID':             ('identifier', int),
                   'Type':                 ('consequence', str),
                   'Name':                 ('name', str),
                   'ClinicalSignificance': ('impact', str),
                   'RS (dbSNP)':           ('RS', int),
                   'PhenotypeList':        ('phenotype', lambda s: s.split('|')),
                   'NumberSubmitters':     ('N_sub', int)}
        variants[data['GeneSymbol']] = dict(symbol=data['GeneSymbol'],
                                            **mutation,
                                            **{masterkey: fun(data[k]) for k, (masterkey, fun) in mapping.items()}
                                            )
# https://www.ncbi.nlm.nih.gov/snp/rs80357346)
import json

with open('clinvar.json', 'w') as fh:
    json.dump(variants, fh)