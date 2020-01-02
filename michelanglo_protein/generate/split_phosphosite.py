import os, csv, gzip, json
from ..settings_handler import global_settings
from collections import defaultdict

licence_note = """
Unfortunately, you have to agree to the CC-by licence at https://www.phosphosite.org/staticDownloads at phosphosite.
Then you have to manually download all the files with _site_dataset.gz.
Afterwards run the settings method retrieve_references.
"""

__doc__ = """
The old code was slow:

    def get_PTM(self):
        assert self.uniprot, 'Uniprot Acc. required. Kind of.'
        modified_residues = []
        for f in os.listdir(self.settings.reference_folder):
            if '_site_dataset' in f and '.gz' not in f:  # it's a Phosphosite plus entry.
                with open(os.path.join(self.settings.reference_folder, f)) as fh:
                    next(fh)  # date
                    next(fh)  # licence
                    next(fh)  # blankline
                    for row in csv.DictReader(fh, delimiter='\t'):
                        if row['ACC_ID'] == self.uniprot: ## this will not pick up mice!
                            modified_residues.append(row["MOD_RSD"])
        self.features['PSP_modified_residues'] = modified_residues ## list of str (e.g. 'K30-m2')
        
        """+licence_note

class Phosphosite:
    def __init__(self):
        self.settings = global_settings
        self.sources = []
        self.data = defaultdict(list)
        for f in os.listdir(self.settings.reference_folder):
            if '_site_dataset' in f:  # it's a Phosphosite plus entry.
                self.sources.append(os.path.join(self.settings.reference_folder,f))
        assert len(self.sources), "No source files found. "+licence_note

    def split(self):
        for f in self.sources:
            if self.settings.verbose:
                print(f'Parsing: {f}')
            with gzip.open(f, 'rt') as fh:
                for row in fh:
                    # there is junk at the top.
                    if row.strip() == '':
                        break
                for row in csv.DictReader(fh, delimiter='\t'):
                    if row['ORGANISM'] == 'human' and '-' in row['MOD_RSD']:
                        ## PSP contains residue info that may not be changes in human, while are in mice.
                        res, ptm = row['MOD_RSD'].split('-')
                        self.data[row['ACC_ID']].append({'symbol': row['GENE'],
                                                    'residue_index': int(res[1:]),
                                                    'from_residue': res[0],
                                                    'ptm': ptm,
                                                    'count': sum([int(row[k]) if row[k] != '' else 0 for k in ('LT_LIT', 'MS_LIT', 'MS_CST')])})
                    #if row['ACC_ID'] == self.uniprot:  ## this will not pick up mice!
                    #    modified_residues.append(row["MOD_RSD"])
        return self

    def write(self,folder='phosphosite'):
        full=os.path.join(global_settings.temp_folder,folder)
        if not os.path.exists(full):
            os.mkdir(full)
        for gene in self.data:
            with open(os.path.join(full,gene+'.json'),'w') as w:
                json.dump(self.data[gene],w)


