from typing import *

import re
import logging
import requests
from collections import Counter


# not std lib imports in methods --> bad typehinting

class Consurfer:
    """
    Get and parse consurf data.
    It takes 1.7 seconds to fetch a grades file off the web.
    I do not know 100% if this is within the T&Cs of Consurf, but I believe it is.
    It definitely isnt if commercial use though.

    Either from web

    >>> cp = ConsurfPoser.from_web('1UBQ', 'A')

    or filename

    >>> cp = ConsurfPoser.from_filename('grades.txt')

    Data is in self.data (dict of ``MET1:A`` to dict of values (see ``ConsurfPoser.keys`` class attribute))
    and can be converted to pandas:

    >>> cp = ConsurfPoser.from_web('1UBQ', 'A')
    >>> cp.to_pandas()

    If a residue appears in SEQPOS but no ATOM records are present, they will be like ``cp.data['___1:A']``.
    The key is the ``3LATOM`` field, this is ATOM numbering, while POS is the SEQPOS numbering.
    ``apply_offset_from_swissmodel`` uses the latter and makes both Uniprot numbering.

    Also can add a consurf conservation to a PyRosetta pose in place.

    >>> cp.add_bfactor_to_pose(pose)

    Or a pymol object

    also, if chain number differs, e.g. V in consurf grades file and A in pose:

    >>> cp.remap_chains({'B': 'A'})

    Likewise with offset

    >>> cp.offset({'A': 10})

    If the Uniprot id is known, the offset can be taken from Swissmodel

    >>> self.apply_offset_from_swissmodel(uniprot, code, chain)

    Potentially support multi-chain operations, but not tested.

    >>> cp = ConsurfPoser.merge([cp1, cp2, cp3])

    """
    log = logging.getLogger()

    # ----- init methods ---------

    def __init__(self):
        self.request = requests.session()
        self.grades_block = ''
        self.data = {}
        self.present_chain = 'A'  # fallback
        self.code = None

    @classmethod
    def from_web(cls, code: str, chain: str):
        """
        Consurf DB. Two requests.
        """
        self = cls()
        self.fetch(code, chain)
        return self

    @classmethod
    def from_filename(cls, grades_filename: str):
        """
        dot grades file.
        """
        self = cls()
        self.read(grades_filename)
        return self

    @classmethod
    def merge(cls, consufers: List['self']):
        """
        Merging into one.
        """
        self = cls()
        # shoddy chaining.
        self.data = {key: values for cs in consufers for key, values in cs.items()}
        return self

    keys = ['POS', 'SEQ', '3LATOM', 'SCORE', 'COLOR', 'CONFIDENCE INTERVAL', 'CONFIDENCE COLORS', 'MSA DATA',
            'RESIDUE VARIETY']

    # ----- inner methods: parse ---------

    def parse(self) -> dict:
        """
        parses a grades block and fills self.conscores
        """
        # self.data = {} # outer key is MET1:A
        assert len(self.grades_block), 'self.grades_block is empty'
        for row in self.grades_block.split('\n'):
            row = row.strip()
            # if pos is not a number it is not a row.
            if not len(row) or not row[0].isdigit():
                continue
            parts = dict(zip(self.keys, [v.strip() for v in row.split('\t') if len(v)]))
            name = parts['3LATOM'].strip()
            if name == '-':
                # missing density
                # note: uses fallback present_chain
                name = '___' + parts['POS'] + ':' + self.present_chain
            else:
                self.present_chain = name[-1]
            self.data[name] = parts
        return self.data

    # ----- inner methods: getters ---------

    def get_conscore(self, key) -> float:
        """
        Key is the ATOM resn resi : Chain
        """
        if key in self.data:
            return float(self.data[key]['SCORE'])
        else:
            return float('nan')

    def get_residue_chain(self, res: str) -> str:
        return res[-1]

    def get_residue_name(self, res: str) -> str:
        return res[:3]

    def get_residue_index(self, res: str) -> int:
        return int(res[3:-2])

    def get_color(self, res: str) -> int:
        # COLOR is an integer... A wee more resolution with confidence range!
        conf_colors = self.data[res]['CONFIDENCE COLORS']
        color = sum(map(int, conf_colors.strip().split(','))) / 2
        # 9 is conserved 1 is not.
        return color

    def get_key(self, index: int, chain: Optional[str] = None) -> str:
        for entry in self.data:
            if int(index) != self.get_residue_index(entry):
                continue
            elif chain is None or chain == self.get_residue_chain(entry):
                return entry
            else:
                continue
        else:
            raise ValueError(f'Absent {index}:{chain}')

    def get_variety(self, res: str) -> List[str]:
        return self.data[res]['RESIDUE VARIETY'].split(',')

    # ----- inner methods: utils ---------

    def to_pandas(self) -> 'pd.DataFrame':
        import pandas as pd
        return pd.DataFrame(self.data).transpose()

    def add_bfactor_to_pose(self, pose: 'pyrosetta.Pose', strict: bool = True):
        """
        Adds the conscores dictionary to bfactor in place.

        If ``strict`` is True, an error is raised if a residue mismatches.
        """
        pdb_info = pose.pdb_info()
        assert not pdb_info.obsolete(), 'The pdb_info is obsolete'
        pdb2pose = pdb_info.pdb2pose
        for con_res in self.data:
            chain = self.get_residue_chain(con_res)
            resi = self.get_residue_index(con_res)
            resn = self.get_residue_name(con_res)
            pose_res = pdb2pose(res=resi, chain=chain)
            if pose_res == 0:
                continue
            residue = pose.residue(pose_res)
            if strict:
                assert residue.name3() == resn, f'{residue.name3()} ≠ {resn}'
            for i in range(residue.natoms()):
                pdb_info.bfactor(pose_res, i, self.get_conscore(con_res))

    def add_bfactor_to_pymol(self, pymol, add_score: bool = False, missing: float = 1.):
        """
        pymol is a pymol2.PyMOL object.

        >>> with pymol2.PyMOL() as pymol:
        >>>     pymol.cmd.fetch('3CMD', 'model')
        >>>     cp.offset_pymol()
        >>>     pymol.cmd.save('test.pse')

        ``add_score`` True adds the score (zero centered), else the 10-1 color
        """
        pymol.cmd.alter(f'*', f'b={missing}')

        for key, values in self.data.items():
            resi = self.get_residue_index(key)
            chain = self.get_residue_chain(key)
            # COLOR is an integer... A wee more resolution with confidence range!
            color = sum(map(int, values['CONFIDENCE COLORS'].strip().split(','))) / 2
            # 9 is conserved 1 is not.
            if add_score:
                color = self.get_conscore(key)
            else:
                # 9 is conserved 1 is not.
                color = 10 - self.get_color(key)
            pymol.cmd.alter(f'resi {resi} and chain {chain}', f'b={color}')
        pymol.cmd.sort()

    def add_bfactor_via_pymol(self, coordinates: str, add_score: bool = False):
        import pymol2
        with pymol2.PyMOL() as pymol:
            pymol.cmd.read_pdbstr(coordinates, 'model')
            self.add_bfactor_to_pymol(pymol, add_score=add_score)
            return pymol.cmd.get_pdbstr()

    # ----- dependent methods: web

    def fetch(self, code: str, chain: str) -> dict:
        """Uses https://consurfdb.tau.ac.il/ so make sure you fall within its usage."""
        first = self._fetch_initial(code, chain)
        self.code = first
        self.grades_block = self._fetch_final(first)
        self.parse()

    def _fetch_initial(self, code: str, chain: str) -> str:
        reply = self.request.get('https://consurfdb.tau.ac.il/scripts/chain_selection.php',
                                 params=dict(pdb_ID=code.upper()))
        self.assert_reply(reply, msg=f'matching {code}')
        mapping = dict(re.findall('option value="(\w) (\w{5})"', reply.text))
        if 'No chains found for' in reply.text:
            self.log.debug(f'Reply: {reply.text.strip()}')
            raise KeyError(f'{code} has no chains according to Consurf')
        elif len(mapping) == 1:
            return list(mapping.values())[0]
        elif chain not in mapping:
            self.log.debug(f'Reply: {reply.text.strip()}')
            raise KeyError(f'Chain {chain} is absent in {code} according to Consurf (SMTL changed the chain number)')
        else:
            return mapping[chain]

    def _fetch_final(self, final: str):
        url = f'https://consurfdb.tau.ac.il/DB/{final}/consurf_summary.txt'
        reply = self.request.get(url)
        self.assert_reply(reply, msg=f'retrieval of suggestion {final}')
        return reply.text

    def assert_reply(self, reply, msg):
        if 'The requested URL was rejected' in reply.text:
            requests.exceptions.ConnectionError(
                f'{msg} gave an error code 200 with Consurf.',
                request=reply.request,
                response=reply)
        if reply.status_code == 200:
            pass
        elif reply.status_code == 404:
            raise requests.exceptions.InvalidURL(f'{msg} could not be found in Consurf',
                                                 request=reply.request,
                                                 response=reply)
        else:
            raise requests.exceptions.ConnectionError(
                f'{msg} gave a code {reply.status_code} with Consurf ({reply.text})',
                request=reply.request,
                response=reply)

    # ----- dependent methods: file

    def read(self, grades_filename: str) -> dict:
        with open(grades_filename) as r:
            self.grades_block = r.read()
        self.parse()

    # ---- offset correction

    def remap_chains(self, chain_map: Dict[str, str]) -> None:
        """
        Say conscores is chain A based, but the model is not.
        """
        new_data = {}
        for res, values in list(self.data.items()):
            old_chain = self.get_residue_chain(res)
            if old_chain in chain_map:
                new_name = res[:-2] + ':' + chain_map[old_chain]
                new_data[new_name] = values
            else:
                new_data[res] = values
        if self.present_chain in chain_map:
            self.present_chain = chain_map[self.present_chain]
        self.data = new_data

    def get_consurf_chain(self):
        # corrects self.present_chain
        chains = [self.get_residue_chain(res) for res in self.data]
        self.present_chain = Counter(chains).most_common(1)[0][0]
        return self.present_chain

    def offset_atom(self, offset_map: Dict[str, int]) -> None:
        """
        offset the index (the ATOM record numbering)
        """
        new_data = {}
        for res, values in list(self.data.items()):
            chain = self.get_residue_chain(res)
            resi = self.get_residue_index(res)
            resn = self.get_residue_name(res)
            if chain in offset_map:
                offset = offset_map[chain]
                new_name = f'{resn}{resi + offset}:{chain}'
                new_data[new_name] = values
            else:
                new_data[res] = values
        self.data = new_data

    def offset_seqpos(self, offset_map: Dict[str, Union[int, List[int]]]) -> None:
        """
        Offset the index (the SEQPOS record numbering and use that.
        """
        new_data = {}
        for ri, (res, values) in enumerate(list(self.data.items())):
            chain = self.get_residue_chain(res)
            resi = self.get_residue_index(res)
            resn = self.get_residue_name(res)
            if chain in offset_map:
                if isinstance(offset_map[chain], int):
                    offset = offset_map[chain]
                elif isinstance(offset_map[chain], list):
                    offset = offset_map[chain][ri]
                else:
                    raise TypeError
                if isinstance(values['POS'], str):
                    values['POS'] = int(values['POS'].strip())
                if offset is None:
                    continue
                new_resi = values['POS'] + offset
                if new_resi < 1:
                    continue
                values['POS'] = new_resi
                new_name = f'{resn}{new_resi}:{chain}'
                new_data[new_name] = values
            else:
                new_data[res] = values
        self.data = new_data

    def get_offset_from_swissmodel(self, uniprot, code, chain) -> int:
        """
        Swissmodel has the offset from the SEQPOS, SIFTS from the ATOM,
        but if there is no ATOM, it is none, so this is safe.
        """
        sm_data = self.request.get(f'https://swissmodel.expasy.org/repository/uniprot/{uniprot}.json').json()
        segs_data = [structure['chains'] for structure in sm_data['result']['structures'] if
                     structure['template'].upper() == code]
        chain_data = [seg for seg in segs_data[0] if seg['id'] == chain][0]['segments'][0]
        pdb_from = chain_data['pdb']['from']
        uniprot_from = chain_data['uniprot']['from']
        return uniprot_from - pdb_from

    def apply_offset_from_swissmodel(self, uniprot, code, chain) -> int:
        offset = self.get_offset_from_swissmodel(uniprot, code, chain)
        self.offset_seqpos({chain: offset})
        return offset

    @property
    def sequence(self):
        return ''.join(row['SEQ'].strip() for row in self.data.values())

    def align(self, ref_sequence: str) -> Tuple[str, str]:
        from Bio import pairwise2
        # alignment = pairwise2.align.globalxx(ref_sequence, data_sequence)[0]
        # 3 sequential mismatches make a gap favourable:
        alignment = pairwise2.align.globalms(ref_sequence, self.sequence,
                                             2,  # match
                                             -.5,  # mismatch
                                             -1,  # open
                                             0  # extend
                                             )[0]
        ref_al, con_al, _, _, _ = alignment
        self.log.debug('ref al ' + ref_al)
        self.log.debug('con al ' + con_al)
        return ref_al, con_al

    def get_offset_by_alignment(self, ref_sequence: str) -> int:
        """
        Offset by alignment. SEQPOS alignment
        """
        ref_al, con_al = self.align(ref_sequence)
        offset = 0
        for r, c in zip(ref_al, con_al):
            if r == '-':
                offset += 1  # number added to con numbering.
            elif c == '-':
                offset -= 1  # number added to con numbering.
            else:
                break
        return -offset

    def get_offset_vector_by_alignment(self, ref_sequence: str) -> List[int]:
        """
        Offset by alignment. SEQPOS alignment
        """
        ref_al, con_al = self.align(ref_sequence)
        c2r = []
        ri = 0
        for r, c in zip(ref_al, con_al):
            if c == '-' and r == '-':
                raise ValueError('Impossible at ``get_offset_map_by_alignment``')
            elif c != '-' and r != '-':
                c2r.append(ri)
                ri += 1
            elif c != '-' and r == '-':
                c2r.append(None)  # no match.
            elif c == '-' and r != '-':
                ri += 1  # no match for R.
        return [ri - ci if ri is not None else None for ci, ri in enumerate(c2r)]

    def apply_offset_by_alignment(self, ref_sequence: str, chain: Optional[str] = None) -> int:
        if chain is None:
            chain = self.present_chain
        offset_vector = self.get_offset_vector_by_alignment(ref_sequence)
        self.log.debug(f'offset {offset_vector}')
        self.offset_seqpos({chain: offset_vector})
        return offset_vector
