import pickle, os, re, json
from datetime import datetime
from .settings_handler import global_settings #the instance not the class.
import gzip, requests
from michelanglo_transpiler import PyMolTranspiler, pymol2

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
        self.offsets = {} if chain == '*' else {chain: int(offset)} ### this is going to be the only one.
        self.pdb_start = None  # no longer used. TO be deleted.
        self.pdb_end = None   # ditto.
        self.resolution = 0 #: crystal resolution. 0 or lower will trigger special cases
        self.code = code
        self.chain_definitions = [] #filled by SIFT. This is a list with a Dict per chain.
        self.type = type.lower() #: str: rcsb | swissmodel | homologue | www | local
        self.chain = chain #: type str: chain letter or * (all)
        if extra is None:
            self.extra = {}
        else:
            self.extra = extra
        self.coordinates = coordinates #: PDBblock
        self.url = url  ## for type = www or local or swissmodel
        # https://files.rcsb.org/download/{self.code}.pdb does not work (often) while the url is something odd.

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
            r = requests.get(self.url)
        elif self.type == 'www':
            r = requests.get(self.url)
        elif self.type == 'local':
            self.coordinates = open(self.url).read()
            return self.coordinates
        else:
            warn(f'Model type {self.type}  for {self.id} could not be recognised.')
            return None
        if r.status_code == 200:
            self.coordinates = r.text
        else:
            warn(f'Model {self.code} failed.')
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
            ### fix this horror.
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
                ## clean rows
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
                ## get offset
                if detail['PDB_BEG'] is not None:  ##nice.
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
        with pymol2.PyMOL() as pymol:
            pymol.cmd.set('fetch_path', os.path.join(self.settings.temp_folder, 'PDB'))
            pymol.cmd.fetch(self.code)
            target = sequence[chain_detail['SP_BEG'] - 1: chain_detail['SP_BEG'] + 10]
            if len(target) == 0:
                return 0
            atoms = pymol.cmd.get_model(f"chain {chain_detail['CHAIN']} and pepseq {target}")
            prev = 'X'
            for atom in atoms.atom:
                if atom.resn in aa:
                    for i in range(1, 10):
                        if aa[atom.resn] == target[i] and target[i - 1] == prev:
                            return (i + chain_detail["SP_BEG"]) - atom.resi

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
