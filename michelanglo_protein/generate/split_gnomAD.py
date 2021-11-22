# This script is actually independent of settings.

import re
import vcf
from collections import defaultdict
import json
import os
from typing import *

folder: str = '/well/gel/HICF2/ref/alleleFrequencies/gnomAD/v3.1.1/'

vcf_filenames: List[str] = [os.path.join(folder, file) for file in os.listdir(folder) if all([
    os.path.splitext(file)[1] == '.gz',
    '.PASS.vcf.gz' in file
])]


class NonSNPError(Exception):
    pass


class gnomADSplitter:
    totals = [  # https://gnomad.broadinstitute.org/help/what-populations-are-represented-in-the-gnomad-data
        {'name': 'African/African American', 'counts': 4554, 'abbreviation': 'afr'},
        {'name': 'Amish', 'counts': 30, 'abbreviation': 'ami'},
        {'name': 'Latino/Admixed American', 'counts': 2345, 'abbreviation': 'amr'},
        {'name': 'Ashkenazi Jewish', 'counts': 68, 'abbreviation': 'asj'},
        {'name': 'East Asian', 'counts': 1215, 'abbreviation': 'eas'},
        {'name': 'European (Finnish)', 'counts': 2750, 'abbreviation': 'fin'},
        {'name': 'Middle Eastern', 'counts': 123, 'abbreviation': 'mid'},
        {'name': 'European (non-Finnish)', 'counts': 3427, 'abbreviation': 'nfe'},
        {'name': 'South Asian', 'counts': 1558, 'abbreviation': 'sas'},
        {'name': 'Other', 'counts': 395, 'abbreviation': 'oth'},
        {'name': 'XX', 'counts': 6717, 'abbreviation': 'XX'},
        {'name': 'XY', 'counts': 9748, 'abbreviation': 'XY'},
        {'name': 'Total', 'counts': 16465, 'abbreviation': '-'}
    ]

    def __init__(self,
                 vcf_filename: str,
                 out_folder: str = 'gnomAD_split',
                 save_frequency: int = 10_000,
                 non_existent: bool = False):
        vcf_reader = vcf.Reader(filename=vcf_filename,
                                compressed=True
                                )
        self.vep_headers = re.search(r'Format\:\s+(.*)', vcf_reader.infos['vep'].desc).group(1).split('|')
        self.out_folder = out_folder
        self.check_out_folder(non_existent)
        i = 0
        self.parsed_buffer = defaultdict(list)  # buffer
        for record in vcf_reader:
            try:
                summary = self.summarize_record(record)
            except NonSNPError:
                continue
            except Exception as error:
                print(f'{error.__class__.__name__}: {error} for {record} in {vcf_filename}')
            self.parsed_buffer[summary['symbol']].append(summary)
            i += 1
            if i == save_frequency:
                print('Saving ', *self.parsed_buffer.keys())
                self.save()
                i = 0

    def check_out_folder(self, non_existent=True):
        if not os.path.exists(self.out_folder):
            os.mkdir(out_folder)
        elif non_existent:
            raise FileExistsError
        else:
            pass

    def save(self):
        for gene in list(self.parsed_buffer.keys()):
            filepath = os.path.join(self.out_folder, gene + '.json')
            if os.path.exists(filepath):
                with open(filepath, 'r') as fh:
                    old = json.load(fh)
            else:
                old = []
            entries = sorted(old + self.parsed_buffer[gene], key=lambda entry: entry['residue_index'])
            with open(filepath, 'w') as fh:
                json.dump(entries, fh)
            del self.parsed_buffer[gene]

    def get_canonical(self, vep_entries: List[str]) -> dict:
        vep_data = [dict(zip(self.vep_headers, entry.split('|'))) for entry in vep_entries if entry]

        canonicals = [datum for datum in filter(lambda e: len(e) > 1, vep_data) if 'CANONICAL' in datum and \
                      datum['CANONICAL'] == 'YES' and \
                      datum['Consequence'] in ('missense_variant', 'stop_gained', 'stop_lost', 'frameshift_variant')
                      ]
        if not canonicals:
            raise NonSNPError
        else:
            return canonicals[0]

    def summarize_vep(self, record: vcf.model._Record) -> dict:
        canonical_vep = self.get_canonical(record.INFO['vep'])
        return dict(symbol=canonical_vep['SYMBOL'],
                    from_residue=canonical_vep['Amino_acids'].split('/')[0],
                    residue_index=int(canonical_vep['Protein_position'].split('-')[0]),
                    to_residue=canonical_vep['Amino_acids'].split('/')[1] if '/' in canonical_vep[
                        'Amino_acids'] else 'X',
                    impact=canonical_vep['IMPACT'],
                    consequence=canonical_vep['Consequence']
                    )

    def summarize_info(self, record: vcf.model._Record) -> dict:
        key_mapping = dict(homozygous='nhomalt_controls_and_biobanks',
                           count='AC_controls_and_biobanks',
                           frequency='AF_controls_and_biobanks')
        summary = {new: int(record.INFO[old][0]) if old in record.INFO else 0 for new, old in key_mapping.items()}
        # ethnicity subs
        for total in self.totals:  # totals: List[Dict] . Used for abbr only...
            abbr = total['abbreviation']
            if total['name'] in ('XX', 'XY', 'Total', 'Amish'):
                continue
            elif 'AF_controls_and_biobanks_' + abbr not in record.INFO:
                continue
            summary[f'frequency {abbr}'] = record.INFO['AF_controls_and_biobanks_' + abbr][0]
        return summary

    def summarize_record(self, record: vcf.model._Record) -> dict:
        # if record.FILTER != 'PASS': raise NonSNPError
        return {'identifier': record.ID,
                **self.summarize_vep(record),
                **self.summarize_info(record)
                }

if __name__ == '__main__':
    from multiprocessing import Pool

    with Pool(15) as p:
        p.map(gnomADSplitter, vcf_filenames)