import pickle, os, re, json
from datetime import datetime
from .settings_handler import global_settings  # the instance not the class.
import gzip, requests
from michelanglo_transpiler import PyMolTranspiler  # called by get_offset_coordinates
from collections import defaultdict
import pymol2
from collections import Counter

from warnings import warn
from .metadata_from_PDBe import PDBMeta
from typing import *

import logging

log = logging.getLogger()


class Structure:
    # Structure as in "protein structure" not C++ structure
    """
    No longer a namedtuple.
    """
    settings = global_settings
    important_attributes = ['x', 'y', 'id', 'description', 'resolution', 'extra', 'alignment']
    temporary_folder = 'temp'

    # __slots__ = ['id', 'description', 'x', 'y', 'url','type','chain','offset', 'coordinates', 'extra', 'offset_corrected']
    def __init__(self, id, description, x: int, y: int, code, type='rcsb', chain='*', offset: int = 0, coordinates=None,
                 extra=None, url=''):
        """
        Stores the structural data for easy use by FeatureViewer and co. Can be converted to StructureAnalyser
        type = rcsb | swissmodel | homologue | www | local | custom

        The type ``rcsb`` isnt called pdb as that would be ambiguous w/ the format
        """
        self.id = id  #: RCSB code
        self.description = description  #: description
        self.x = int(x)  #: resi in the whole uniprot protein
        self.y = int(y)  #: end resi in the whole uniprot protein
        self.offset_corrected = False  # prevents a double trip
        # TODO these do not seem to used as are overwritten by chain definitions.
        if offset is None:
            self.offset = None
        else:
            self.offset = int(offset)  #: offset is the number *subtracted* from the PDB index
            # to make it match the position in Uniprot.
        self.offsets = {} if chain == '*' else {chain: self.offset}  ## this is going to be the only one.
        self.pdb_start = None  # no longer used. TO be deleted.
        self.pdb_end = None  # ditto.
        self.resolution = 0  #: crystal resolution. 0 or lower will trigger special cases
        self.code = code
        self.chain_definitions = []  # filled by SIFT. This is a list with a Dict per chain.
        self.type = type.lower()  #: str: rcsb | swissmodel | homologue | www | local | custom | alphafold2
        self.chain = chain  #: type str: chain letter or * (all)
        self.alignment = {}
        if extra is None:
            self.extra = {}
        else:
            self.extra = extra
        self.coordinates = coordinates  #: PDBblock
        self.url = url  # for type = www or local or swissmodel
        # https://files.rcsb.org/download/{self.code}.pdb does not work (often) while the url is something odd.

    def is_satisfactory(self, resi: int):
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
                assert pymol.cmd.select(f'chain A and resi {resi} and name {name}'), \
                    f'Given protein has no {name} atom in residue {resi} in chain A'

    def to_dict(self, full=False) -> Dict:
        if full:
            extras = {key: getattr(self, key) for key in self.important_attributes if hasattr(self, key)}
        else:
            extras = {}
        return {'x': self.x, 'y': self.y, 'id': self.id,
                'type': self.type, 'description': self.description, **extras}

    def __str__(self):
        return str(self.to_dict())

    def get_coordinates(self) -> str:
        """
        Gets the coordinates (PDB block) based on ``self.url`` and ``self.type``
        :return: coordinates
        :rtype: str
        """
        # custom/w coordinates
        if self.coordinates:
            return self.coordinates
        elif self.type == 'custom':  # provided.
            # self.coordinates was empty
            raise ValueError('No coordinates provided for custom retrieval')
        # url is filepath
        elif self.url and self.type == 'local':
            self.coordinates = open(self.url).read()
            return self.coordinates
        elif self.type == 'local':
            # self.url was empty
            raise ValueError('No filepath provided for local retrieval')
        # url present
        elif self.url:  # regardless of type/
            r = requests.get(self.url, allow_redirects=True)
        elif self.type in ('www', 'alphafold2'):
            assert self.url, 'No URL provided for www retrieval'
            r = requests.get(self.url)
        # other
        elif self.type == 'rcsb':
            r = requests.get(f'https://files.rcsb.org/download/{self.code}.pdb')
        elif self.type == 'swissmodel':
            assert self.url, 'No URL provided for SWISSMODEL retrieval'
            r = requests.get(self.url, allow_redirects=True)
        else:
            raise ValueError(f'Model type {self.type}  for {self.id} could not be recognised.')
        # --- read reply
        if r.status_code == 200:
            self.coordinates = r.text
        else:
            raise ConnectionError(f'Model {self.code} ({self.url}) failed.')
        return self.coordinates

    def get_offset_coordinates(self, sequence: Optional[str] = None):
        """
        Gets the coordinates and offsets them.
        :return:
        """
        # if self.offset_corrected:
        #     return self.coordinates
        log.debug(self.chain_definitions)
        if not self.chain_definitions:
            self.lookup_sifts()
        self.coordinates = PyMolTranspiler().renumber(self.get_coordinates(),
                                                      self.chain_definitions,
                                                      sequence=sequence,
                                                      make_A=self.chain,
                                                      remove_solvent=True).raw_pdb
        self.fix_renumbered_annotation()
        return self.coordinates

    def fix_renumbered_annotation(self):
        # this should have logging.
        if self.chain != 'A':
            ## fix this horror.
            log.warning('Chain A is not the target chain in the definitions!')
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
        # just in case there is a double trip! (Non server usage)
        self.chain = 'A'
        self.offset = 0
        self.offset_corrected = True
        for i, c in enumerate(self.chain_definitions):
            if c['chain'] == 'A':
                self.chain_definitions[i]['offset'] = 0

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
            return self  # it is probably clean.
        log.debug(f'Looking up sifts if this is empty: {self.chain_definitions}')
        if not self.chain_definitions:
            details = self._get_sifts()
            offset = 0
            for detail in details:
                # clean rows
                for k in ('PDB_BEG', 'PDB_END', 'RES_END', 'RES_BEG', 'SP_BEG', 'SP_END'):
                    if k == 'None' or k is None:
                        detail[k] = None
                    elif isinstance(detail[k], int):
                        pass  # this means so test is being done.
                    else:
                        r = re.search('(-?\d+)', detail[k])  # str().isdigit() does not like negatives.
                        if r is None:
                            detail[k] = None
                        else:
                            detail[k] = int(r.group(1))  # yes. py int is signed
                # get offset
                if detail['PDB_BEG'] is not None:  # nice.
                    offset = detail['SP_BEG'] - detail['PDB_BEG']
                elif detail['PDB_END'] is not None:
                    offset = detail['SP_BEG'] - (detail['PDB_END'] - (detail['SP_END'] - detail['SP_BEG']))
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

    def _get_sifts(self, all_chains=True):  # formerly called .lookup_pdb_chain_uniprot
        details = []
        headers = 'PDB     CHAIN   SP_PRIMARY      RES_BEG RES_END PDB_BEG PDB_END SP_BEG  SP_END'.split()
        with self.settings.open('pdb_chain_uniprot') as fh:
            for row in fh:
                if self.code.lower() == row[0:4]:
                    entry = dict(zip(headers, row.split()))
                    if self.chain == entry['CHAIN'] or all_chains:
                        details.append(entry)
        return details

    def get_offset_from_PDB(self, chain_detail: Dict, sequence: str) -> int:
        """
        This is used by sandbox. if transition to swissmodel data works this will be removed.
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
            for begin_offset in range(0, len(sequence) - begin, 10):
                # Try full size window
                window = 50
                target = sequence[begin - 1 + begin_offset: begin + window + begin_offset]
                if len(target) == 0:
                    # Under what conditions does this happen???
                    debugprint(
                        f'sequence is {len(sequence)}, while range is {begin + begin_offset}-{end + begin_offset}')
                    continue
                sele_target = f"chain {chain_detail['CHAIN']} and pepseq {target} and name CA"
                # Shrink window to account for failed selection due to weird atoms or short peptides
                while pymol.cmd.select(sele_target) == 0:
                    window -= 10
                    if window < 10:
                        debugprint(f'{self.code} ({chain_detail}) does not contain this sequence ({target})')
                        break  # double continue
                    target = sequence[begin + begin_offset - 1: begin + begin_offset + window]
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
                        # print(atom.resn, atom.resi, target)
                        for i in range(1, window):  # in case there are missing parts.
                            if aa[atom.resn] == target[i] and target[i - 1] == prev:
                                # print(f'YATTA {(i + chain_detail["SP_BEG"]) - int(atom.resi)}')
                                return (i + int(chain_detail["SP_BEG"])) - int(atom.resi)
                        prev = aa[atom.resn]
                    else:
                        prev = 'X'
                    prev_i = int(atom.resi)
                else:  # more than 50 aa without coordinates at the N terminus? Engineered residues.
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
        # code not used anywhere.
        return PDBMeta(self.code + '_' + self.chain).data

    @classmethod
    def from_swissmodel_query(cls, structural_data: dict, uniprot: str):
        # called by retrieve_structures_from_swissmodel of ProteinCore
        # structural_data is an entry in the list data['result']['structures'] from a JSON query to swissmodel
        keepers = ['coverage', 'created_date', 'from', 'gmqe', 'identity', 'in_complex_with', 'ligand_chains', 'method',
                   'oligo-state', 'qmean', 'similarity', 'template', 'to']
        # offset when downloaded from Expasy not PDB
        if structural_data['provider'] == 'PDB':
            offset = structural_data['chains'][0]['segments'][0]['uniprot']['from'] - \
                     structural_data['chains'][0]['segments'][0]['pdb']['from']
        else:
            offset = 0
        chain_id = 0  # is this ever not the zeroth?
        structure = cls(
            # these two do ziltch:
            id=structural_data['md5'],
            description=structural_data['coordinates'].split('/')[-1],
            # this is shown by venus
            code=structural_data['template'] if structural_data[
                                                    'provider'] == 'PDB' else f'based upon {structural_data["template"]}',
            # range in uniprot
            x=structural_data['from'],
            y=structural_data['to'],
            # offset not provided by Swissmodel. get_offset_from_PDB is the only option.
            offset=offset,
            type='swissmodel' if structural_data['provider'] == 'SWISSMODEL' else 'rcsb',
            url=structural_data['coordinates'],
            # structural_data['template'][-1] works only for swissmodel for chain.
            chain=structural_data['chains'][chain_id]['id'],
            extra={k: structural_data[k] for k in keepers if k in structural_data}
        )

        # do not fill the structure.chain_definitions via SIFT
        def get_chain_def(info) -> List[dict]:
            chain = info['id']
            chained = [dict(chain=chain,
                            x=structure.x,
                            y=structure.y,
                            offset=structure.offset,
                            range=f'{structure.x}-{structure.y}',
                            uniprot=uniprot)]
            if 'in_complex_with' in structural_data:
                for partner_chain, partner_def in structural_data['in_complex_with'].items():
                    chained.append(dict(chain=partner_chain,
                                offset=0,
                                description=partner_def[0]['description'],
                                uniprot=partner_def[0]['uniprot_ac'] if 'uniprot_ac' in partner_def[0] else 'P00404',
                                ))
            return chained

        structure.chain_definitions = [chained for info in structural_data['chains'] for chained in get_chain_def(info)]

        # ---- sequence

        def get_sequence(part):
            # the gap needs to be for 'uniprot' only!
            seq = '-' * (structural_data['chains'][chain_id]['segments'][0]['uniprot']['from'] - 1)
            seq += structural_data['chains'][chain_id]['segments'][0][part]['aligned_sequence']
            return seq

        if structure.type == 'swissmodel':
            # herein the template is
            structure.alignment = {'template': get_sequence('smtl'), 'uniprot': get_sequence('uniprot')}
        elif structure.type == 'rcsb':
            structure.alignment = {'template': get_sequence('pdb'), 'uniprot': get_sequence('uniprot')}
        return structure

    # def Xget_coordinates_w_template_extras(self, sequence: Optional[str] = None, monomer=True):
    #     """
    #     Sequence is the Uniprot sequence. For safety/debug
    #     Biological assembly is a tricky one. Therefore it is often safer to just use a monomer.
    #     in_struct_asyms and in_chains are different if bio assembly is smaller than async assembly
    #     """
    #     assert self.type == 'swissmodel'
    #     template_code, chain = re.search('(\w{4})\.\w+\.(\w)', self.code).groups()
    #     pdbblock = self.get_coordinates()
    #     meta = PDBMeta(template_code, chain)
    #     other_chains = meta.get_other_chains(chain, first_only=monomer) - set(self.extra['in_complex_with'].keys())
    #     if 'ligand_chains' in self.extra:
    #         present_ligands = {e['hetid'] for e in self.extra['ligand_chains'] if 'hetid' in e}
    #     else:
    #         present_ligands = set()
    #     other_ligands = meta.get_interesting_ligand_names() - present_ligands
    #     if len(other_chains) + len(other_ligands) > 0:
    #         with pymol2.PyMOL() as pymol:
    #             pymol.cmd.read_pdbstr(self.coordinates, 'threaded')
    #             pymol.cmd.fetch(template_code, 'template', file=None)
    #             pymol.cmd.remove(f'template and chain {chain} and polymer')
    #             interesting_lig = ' or '.join([f'resi {name3}' for name3 in other_ligands])
    #             wanted_chains = ' or '.join([f'chain {chain}' for chain in other_chains])
    #             chain_sele = f'(template and ({wanted_chains}))'
    #             lig_sele = f'(template and ({interesting_lig}))'
    #             if wanted_chains and interesting_lig:
    #                 pymol.cmd.create('combo', f'threaded or {chain_sele} or {lig_sele}')
    #             elif wanted_chains:
    #                 pymol.cmd.create('combo', f'threaded or {chain_sele}')
    #             elif interesting_lig:
    #                 pymol.cmd.create('combo', f'threaded or {lig_sele}')
    #             else:
    #                 pymol.cmd.create('combo', 'threaded')  # this should be a break statement with logging
    #                 # log.critical('Impossible')
    #             pymol.cmd.remove('(byres polymer around 5) and not polymer')  # remove asymmetric non bio ligands.
    #             pdbblock = pymol.cmd.get_pdbstr('combo')
    #             # update chain definitions
    #             for entity in meta.get_other_proteins(chain):
    #                 for chain in entity['in_chains']:
    #                     if chain not in other_chains:
    #                         continue
    #                     if pymol.cmd.select(f'combo and chain {chain}') == 0:
    #                         continue
    #                     if 'mappings' in entity['source'][0]:
    #                         x = entity['source'][0]['mappings'][0]['start']['residue_number']
    #                         y = entity['source'][0]['mappings'][0]['end']['residue_number']
    #                     else:
    #                         x = 1
    #                         y = 1_000_000
    #                     self.chain_definitions.append(dict(chain=chain,
    #                                                        x=x, y=y,
    #                                                        offset=0,  # who cares.
    #                                                        range=f'{x}-{y}',
    #                                                        uniprot='P00404',  # unknown.
    #                                                        name=entity['molecule_name'][0],
    #                                                        transplanted=True
    #                                                        ))
    #                     if monomer:
    #                         break  # only first entity['in_chains'] was added
    #             for entity in meta.get_other_polymers(chain):
    #                 if meta.is_peptide(entity):
    #                     continue
    #                 for chain in entity['in_chains']:
    #                     if pymol.cmd.select(f'combo and chain {chain}') == 0:
    #                         continue
    #                     self.chain_definitions.append(dict(chain=chain,
    #                                                        x=1, y=entity['length'],
    #                                                        offset=0,  # who cares.
    #                                                        range=f"1-{entity['length']}",
    #                                                        uniprot='P00404',  # not valid
    #                                                        name=entity['molecule_name'][0],
    #                                                        transplanted=True
    #                                                        ))
    #     # pro-forma.
    #     self.fix_renumbered_annotation()
    #     return self.coordinates

    def get_coordinates_w_template_extras(self, sequence: Optional[str] = None):
        """
        Sequence is the Uniprot sequence. For safety/debug
        Biological assembly is a tricky one. Therefore it is often safer to just use a monomer.

        in_struct_asyms and in_chains are different if bio assembly is smaller than async assembly

        The Swissmodel template library contains segi and chain renumbered structures
        preferably as biounits. This is great except for the fact that the data of the chains is not available.
        Hence the hybrid PDBe approach.
        This means that for 6pax, whereas the (segi, chains) are

        ([('B', 'A'), ('C', 'B'), ('A', 'C'), ('B', 'D'), ('C', 'E'), ('A', 'F')]

        In the SMTL 6pax.1 it is:

        [('-', ''), ('A', ''), ('B', ''), ('C', '')])

        """
        assert self.type == 'swissmodel'
        template_code, chain = re.search('(\w{4})\.\w+\.(\w)', self.code).groups()
        log.debug(f'template_code={template_code} chain={chain} current={self.chain}')
        pdbblock = self.get_coordinates()
        meta = PDBMeta(template_code, chain)
        with pymol2.PyMOL() as pymol:
            pymol.cmd.set('fetch_path', self.temporary_folder)
            pymol.cmd.read_pdbstr(pdbblock, 'threaded')
            try:
                pymol.cmd.fetch(template_code, 'template', file=None, type='pdb1')
            except Exception as err:
                # Nowadays pymol2.pymol.CmdException inherits Exception and not BaseException directly
                log.warning(f'Error caused with pdb1. {err.__class__.__name__}')
                pymol.cmd.fetch(template_code, 'template', file=None)
            log.debug('Merging with template...')
            pymol.cmd.remove('solvent')
            overlap_iter = pymol.cmd.get_model('template and (threaded around 0.1)').atom
            data = Counter([(atom.chain, atom.segi) for atom in overlap_iter])
            removed_chains = []
            for (chain, segi), count in data.most_common():
                if count < 10:
                    pass  # a couple of overlapping atoms are fine.
                else:
                    if chain and segi:
                        pymol.cmd.remove(f'template and chain {chain} and segi {segi} and polymer')
                    elif chain:
                        pymol.cmd.remove(f'template and chain {chain} and polymer')
                    else:
                        pymol.cmd.remove(f'template and chain "" and and polymer')
                    removed_chains.append(chain)
                    meta.remove_chain(chain)
            log.debug(f'The following chains were removed {removed_chains}')
            pymol.cmd.remove(f'template and (threaded around 0.1) and not polymer')
            boring = '+'.join(meta.boring_ligands)
            if boring:
                pymol.cmd.remove(f'template and resn {boring}')
            # correct template chain A.
            # thanks to the fact that SM do not have segi.
            if self.chain != 'A' and pymol.cmd.select(f'threaded and chain A'):
                pymol.cmd.alter(f'threaded and chain A', 'chain="X"')
                pymol.cmd.sort()
                for defi in self.chain_definitions:
                    if defi['chain'] == 'A':
                        defi['chain'] = 'X'
                        break
                    else:
                        log.debug(f'Chain absent in definition')
            if self.chain != 'A':
                pymol.cmd.alter(f'threaded and chain {self.chain}', 'chain="A"')
                pymol.cmd.sort()
                for defi in self.chain_definitions:
                    if defi['chain'] == self.chain:
                        defi['chain'] = 'A'
                        break
                    else:
                        log.debug(f'... Chain absent in definition')
                self.chain = 'A'
            transpiler = PyMolTranspiler()  # using parts.
            transpiler.pymol = pymol
            chaingen = transpiler.get_new_letter()
            present_chains = {atom.chain for atom in pymol.cmd.get_model('threaded').atom}
            to_be_added_chains = {atom.chain for atom in pymol.cmd.get_model('template').atom}
            for chain in to_be_added_chains:
                if chain in present_chains:
                    old_chain = chain
                    new_chain = next(chaingen)
                    pymol.cmd.alter(f"template and chain {old_chain}", f'chain="{new_chain}"')
                    pymol.cmd.sort()
                    meta.move_chain(old_chain, new_chain)
                elif chain == '':
                    pymol.cmd.alter(f"template and chain ''", f'chain="{new_chain}"')
                    pymol.cmd.sort()
            pymol.cmd.create('combo', 'template or threaded', 1)  # state 1.
            self.coordinates = pymol.cmd.get_pdbstr('combo')
            if False:  # extreme debug time
                pymol.cmd.save('temp.pse')
                raise SystemExit
            # update chain definitions
            for entity in meta.get_proteins():
                for chain in entity['in_chains']:
                    if chain not in to_be_added_chains:
                        continue
                    elif pymol.cmd.select(f'combo and chain {chain}') == 0:
                        continue
                    elif 'mappings' in entity['source'][0]:
                        x = entity['source'][0]['mappings'][0]['start']['residue_number']
                        y = entity['source'][0]['mappings'][0]['end']['residue_number']
                    else:
                        x = 1
                        y = 1_000_000
                    log.debug(f'Adding chain {chain}')
                    self.chain_definitions.append(dict(chain=chain,
                                                       x=x, y=y,
                                                       offset=0,  # who cares.
                                                       range=f'{x}-{y}',
                                                       uniprot='P00404',  # unknown.
                                                       name=entity['molecule_name'][0],
                                                       transplanted=True
                                                       ))
            for entity in meta.get_polymers():
                if meta.is_peptide(entity):
                    continue
                for chain in entity['in_chains']:
                    if pymol.cmd.select(f'combo and chain {chain}') == 0:
                        continue
                    self.chain_definitions.append(dict(chain=chain,
                                                       x=1, y=entity['length'],
                                                       offset=0,  # who cares.
                                                       range=f"1-{entity['length']}",
                                                       uniprot='P00404',  # not valid
                                                       name=entity['molecule_name'][0].replace('*P', '')
                                                                                      .replace('P*', '')
                                                                                      .replace('*', ''),
                                                       transplanted=True
                                                       ))
        # end of pymol context manager
        # pro-forma:
        self.fix_renumbered_annotation()
        log.debug('... Merged')
        return self.coordinates
