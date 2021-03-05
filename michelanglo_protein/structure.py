import pickle, os, re, json
from datetime import datetime
from .settings_handler import global_settings #the instance not the class.
import gzip, requests
from michelanglo_transpiler import PyMolTranspiler
from collections import defaultdict
import pymol2

from warnings import warn
from .metadata_from_PDBe import PDBMeta
from typing import Dict

class Structure:
    #lolz. a C++ coder would hate this name. Sturcture as in "protein structure"
    #that is not funny. Why I did I think it was?
    #Why am I talking to my past self?!
    """
    No longer a namedtuple.
    Stores the structural data for easy use by FeatureViewer and co. Can be converted to StructureAnalyser
    type = rcsb | swissmodel | homologue
    """
    settings = global_settings

    #__slots__ = ['id', 'description', 'x', 'y', 'url','type','chain','offset', 'coordinates', 'extra']
    def __init__(self, id, description, x:int, y:int, code, type='rcsb',chain='*',offset:int=0, coordinates=None, extra=None, url=''):
        """
        Stores the structural data for easy use by FeatureViewer and co. Can be converted to StructureAnalyser
        type = rcsb | swissmodel | homologue | www | local
        """
        self.id = id #: RCSB code
        self.description = description #: description
        self.x = int(x)  #: resi in the whole uniprot protein
        self.y = int(y)  #: end resi in the whole uniprot protein
        self.offset = int(offset) #: offset is the number *subtracted* from the PDB index to make it match the position in Uniprot.
        self.offsets = {} if chain == '*' else {chain: int(offset)} ## this is going to be the only one.
        self.pdb_start = None  # no longer used. TO be deleted.
        self.pdb_end = None   # ditto.
        self.resolution = 0 #: crystal resolution. 0 or lower will trigger special cases
        self.code = code
        self.chain_definitions = [] #filled by SIFT. This is a list with a Dict per chain.
        self.type = type.lower() #: str: rcsb | swissmodel | homologue | www | local | custom
        self.chain = chain #: type str: chain letter or * (all)
        if extra is None:
            self.extra = {}
        else:
            self.extra = extra
        self.coordinates = coordinates #: PDBblock
        self.url = url  # for type = www or local or swissmodel
        # https://files.rcsb.org/download/{self.code}.pdb does not work (often) while the url is something odd.

    def is_satisfactory(self, resi:int):
        with pymol2.PyMOL() as pymol:
            pymol.cmd.read_pdbstr(self.coordinates, 'given_protein')
            residex = defaultdict(list)
            note = 'custom protein'
            for atom in pymol.cmd.get_model('name CA').atom:
                residex[atom.chain].append(atom.resi)
            if len(residex) == 1 and 'A' not in residex:
                note += ' - chain moved to A'
                pymol.cmd.alter('given_protein', 'chain="A"')
                pymol.cmd.sort()
                move = list(residex.values())[0]
                residex = {'A': move}
                self.coordinates = pymol.cmd.get_pdbstr()
            if not self.chain_definitions:
                self.chain_definitions = [{'chain': chain,
                                      'uniprot': "XXX",
                                      'x': min(residex[chain]),
                                      'y': max(residex[chain]),
                                      'offset': 0,
                                      'range': f'0-9999',
                                      'name': self.code,
                                      'description': note} for chain in residex]
            assert pymol.cmd.select('given_protein'), 'Given protein had no valid data to load'
            assert pymol.cmd.select('chain A'), 'Given protein has no chain A'
            assert pymol.cmd.select(f'chain A and resi {resi}'), f'Given protein has no residue {resi} in chain A'
            for name in ('N', 'CA', 'C'):
                assert pymol.cmd.select(
                    f'chain A and resi {resi} and name {name}'), f'Given protein has no {name} atom in residue {resi} in chain A'

    def to_dict(self) -> Dict:
        return {'x': self.x, 'y': self.y, 'id': self.id, 'description': self.description}

    def __str__(self):
        return str(self.to_dict())

    def get_coordinates(self) -> str:
        """
        Gets the coordinates (PDB block) based on ``self.url`` and ``self.type``
        :return: coordinates
        :rtype: str
        """
        if self.type == 'rcsb':
            r = requests.get(f'https://files.rcsb.org/download/{self.code}.pdb')
        elif self.type == 'swissmodel':
            assert self.url, 'No URL provided for SWISSMODEL retrieval'
            r = requests.get(self.url, allow_redirects=True)
        elif self.type == 'www':
            assert self.url, 'No URL provided for www retrieval'
            r = requests.get(self.url)
        elif self.type == 'local':
            assert self.url, 'No filepath provided for local retrieval'
            self.coordinates = open(self.url).read()
            return self.coordinates
        elif self.type == 'custom': # provided.
            assert self.coordinates, 'No coordinates provided for custom retrieval'
            return self.coordinates
        else:
            raise ValueError(f'Model type {self.type}  for {self.id} could not be recognised.')
        if r.status_code == 200:
            self.coordinates = r.text
        else:
            warn(f'Model {self.code} ({self.url}) failed.')
        return self.coordinates

    def get_offset_coordinates(self):
        """
        Gets the coordinates and offsets them.
        :return:
        """
        if not self.chain_definitions:
            self.lookup_sifts()
        self.coordinates = PyMolTranspiler().renumber(self.get_coordinates(), self.chain_definitions, make_A=self.chain).raw_pdb
        if self.chain != 'A':
            ## fix this horror.
            for i, c in enumerate(self.chain_definitions):
                if self.chain_definitions[i]['chain'] == 'A':
                    self.chain_definitions[i]['chain'] = 'XXX'
                    break
            for i, c in enumerate(self.chain_definitions):
                if self.chain_definitions[i]['chain'] == self.chain:
                    self.chain_definitions[i]['chain'] = 'A'
                    break
            for i, c in enumerate(self.chain_definitions):
                if self.chain_definitions[i]['chain'] == 'XXX':
                    self.chain_definitions[i]['chain'] = self.chain
                    break
        return self.coordinates

    def includes(self, position, offset=0):
        """
        Generally there should not be an offset as x and y are from Uniprot data so they are already fixed!
        :param position:
        :param offset:
        :return:
        """
        if self.x + offset > position:
            return False
        elif self.y + offset < position:
            return False
        else:
            return True


    def lookup_sifts(self):
        """
        SIFTS data. for PDBe query see elsewhere.
        There are four start/stop pairs that need to be compared to get a good idea of a protein.
        For a lengthy discussion see https://blog.matteoferla.com/2019/09/pdb-numbering-rollercoaster.html
        Also for a good list of corner case models see https://proteopedia.org/wiki/index.php/Unusual_sequence_numbering
        :return: self
        """
        if self.type != 'rcsb':
            return self #it is probably clean.

        if not self.chain_definitions:
            details = self._get_sifts()
            for detail in details:
                # clean rows
                for k in ('PDB_BEG','PDB_END', 'RES_END', 'RES_BEG', 'SP_BEG','SP_END'):
                    if k == 'None' or k is None:
                        detail[k] = None
                    elif isinstance(detail[k], int):
                        pass # this means so test is being done.
                    else:
                        r = re.search('(-?\d+)', detail[k]) #str().isdigit() does not like negatives.
                        if r is None:
                            detail[k] = None
                        else:
                            detail[k] = int(r.group(1)) #yes. py int is signed
                # get offset
                if detail['PDB_BEG'] is not None:  #nice.
                    offset = detail['SP_BEG'] - detail['PDB_BEG']
                elif detail['PDB_END'] is not None:
                    offset = detail['SP_BEG'] - ( detail['PDB_END'] - (detail['SP_END'] - detail['SP_BEG']))
                elif detail['SP_BEG']:
                    offset = 0
                else:
                    offset = 0
            self.chain_definitions = [{'chain': d['CHAIN'],
                                       'uniprot': d['SP_PRIMARY'],
                                       'x': d["SP_BEG"],
                                       'y': d["SP_END"],
                                       'offset': offset,
                                       'range': f'{d["SP_BEG"]}-{d["SP_END"]}',
                                       'name': None,
                                       'description': None} for d in details]
        try:
            if self.chain != '*':
                detail = next(filter(lambda x: self.chain == x['chain'], self.chain_definitions))
                self.offset = detail['offset']
        except StopIteration:
            warn(f'{self.code} {self.chain} not in {self.chain_definitions}')
            return self
        self.offsets = {d['chain']: d['offset'] for d in self.chain_definitions}
        return self

    def _get_sifts(self, all_chains=True): #formerly called .lookup_pdb_chain_uniprot
        details = []
        headers = 'PDB     CHAIN   SP_PRIMARY      RES_BEG RES_END PDB_BEG PDB_END SP_BEG  SP_END'.split()
        with self.settings.open('pdb_chain_uniprot') as fh:
            for row in fh:
                if self.code.lower() == row[0:4]:
                    entry = dict(zip(headers, row.split()))
                    if self.chain == entry['CHAIN'] or all_chains:
                        details.append(entry)
        return details

    def get_offset_from_PDB(self, chain_detail: Dict, sequence:str) -> int:
        """
        This gets the offset for a chain in the code given a sequence.
        It takes 30 ms to run. However, the sequence is problematic.

        :param chain: see SIFTs {'PDB': '5l8o', 'CHAIN': 'C', 'SP_PRIMARY': 'P51161', 'RES_BEG': 1, 'RES_END': 128, 'PDB_BEG': None, 'PDB_END': None, 'SP_BEG': 1, 'SP_END': 128}
        :type chain: Dict
        :param sequence: sequence of unirpot
        :type sequence: str
        :return: offset
        :rtype: int
        """
        aa = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
              'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
              'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
              'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
        begin = int(chain_detail['SP_BEG'])
        end = int(chain_detail['SP_BEG'])
        assert isinstance(chain_detail, dict), 'Chain detail is a Dict of the specific chain. Not whole protein.'
        debugprint = lambda x: None
        with pymol2.PyMOL() as pymol:
            # Load file
            pymol.cmd.set('fetch_path', os.path.join(self.settings.temp_folder, 'PDB'))
            pymol.cmd.fetch(self.code)
            # Try different windows
            for begin_offset in range(0,len(sequence) - begin, 10):
                # Try full size window
                window = 50
                target = sequence[begin - 1 + begin_offset : begin + window + begin_offset]
                if len(target) == 0:
                    # Under what conditions does this happen???
                    debugprint(f'sequence is {len(sequence)}, while range is {begin+begin_offset}-{end+begin_offset}')
                    continue
                sele_target = f"chain {chain_detail['CHAIN']} and pepseq {target} and name CA"
                # Shrink window to account for failed selection due to weird atoms or short peptides
                while pymol.cmd.select(sele_target) == 0:
                    window -= 10
                    if window < 10:
                        debugprint(f'{self.code} ({chain_detail}) does not contain this sequence ({target})')
                        break #double continue
                    target = sequence[begin +  begin_offset - 1: begin +  begin_offset + window]
                    sele_target = f"chain {chain_detail['CHAIN']} and pepseq {target} and name CA"
                # double continue
                if window < 10:
                    continue
                # Iterate
                atoms = pymol.cmd.get_model(sele_target)
                prev = 'X'
                prev_i = -999
                for atom in atoms.atom:
                    if int(atom.resi) == prev_i:
                        continue
                    if atom.resn in aa:
                        # Simplest case:
                        # if aa[atom.resn] == target[0]:
                        #     return chain_detail["SP_BEG"] - atom.resi
                        # In case there are missing parts.
                        # for i in range(0, 20):
                        #     if aa[atom.resn] == target[i]:
                        #         return (i + chain_detail["SP_BEG"]) - atom.resi
                        # In case there are missing parts and repeated residues.
                        #print(atom.resn, atom.resi, target)
                        for i in range(1, window): #in case there are missing parts.
                            if aa[atom.resn] == target[i] and target[i - 1] == prev:
                                #print(f'YATTA {(i + chain_detail["SP_BEG"]) - int(atom.resi)}')
                                return (i + int(chain_detail["SP_BEG"])) - int(atom.resi)
                        prev = aa[atom.resn]
                    else:
                        prev = 'X'
                    prev_i = int(atom.resi)
                else: #more than 50 aa without coordinates at the N terminus? Engineered residues.
                    debugprint(f'{self.code} More than {window} AA without a match!! {target} {prev}')
                    continue
            warn(f'UTTER FAILURE for {self.code}')
            return 0

    def lookup_resolution(self):
        if self.type != 'rcsb':
            return self
        with self.settings.open('resolution') as fh:
            resolution = json.load(fh)
            for entry in resolution:
                if entry['IDCODE'] == self.code:
                    if entry['RESOLUTION'].strip():
                        self.resolution = float(entry['RESOLUTION'])
                    break
            else:
                warn(f'No resolution info for {self.code}')
        return self

    def lookup_ligand(self):
        warn('TEMP! Returns the data... not self')
        return PDBMeta(self.code+'_'+self.chain).data
