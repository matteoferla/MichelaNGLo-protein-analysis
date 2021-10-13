from ..structure import Structure
from ..mutation import Mutation
from michelanglo_transpiler import PyMolTranspiler
from .consurf import Consurfer
import pymol2
from itertools import product
import math, re, time, logging
import numpy as np
from typing import *


class StructureAnalyser:
    """
    3D details of a mutation based on the structure which is stored in ``.structure`` (which is a ``Bio.PDB.PDBParser`` obj.)
    """
    # Wikipedia Tien et al. 2013 (emp.)
    maxASA = {'A': 121.0, 'R': 265.0, 'N': 187.0, 'D': 187.0, 'C': 148.0, 'E': 214.0, 'Q': 214.0, 'G': 97.0, 'H': 216.0,
              'I': 195.0, 'L': 191.0, 'K': 230.0, 'M': 203.0, 'F': 228.0, 'P': 154.0, 'S': 143.0, 'T': 163.0,
              'W': 264.0, 'Y': 255.0, 'V': 165.0}
    # normal heavy atom count done with the script genenerate_AA_MCS.py This includes OXT atom.
    normal_HA = {'A': 6, 'R': 12, 'N': 9, 'D': 9, 'C': 7, 'E': 10, 'Q': 10, 'G': 5, 'H': 11, 'I': 9, 'L': 9, 'K': 10,
                 'M': 9, 'F': 12, 'P': 8, 'S': 7, 'T': 8, 'W': 15, 'Y': 13, 'V': 8}
    # I think PyMOL has H,S,L only?
    ss_types = {'H': 'Helix', 'S': 'Sheet', 'L': 'Loop', 'G': '3_10 helix', 'I': 'Pi helix', 'T': 'Turn', 'C': 'Coil',
                'E': 'Sheet', 'B': 'Beta bridge', '-': 'Unassigned'}
    log = logging.getLogger()
    error_on_missing_conservation = False

    def __init__(self, structure: Structure, mutation: Mutation, sequence: Optional[str] = None, no_conservation=False):
        """

        :param structure: a instance of Structure, a former namedtuple and is in core.py
        :param mutation: a instance of mutation.
        """
        self.mutation = mutation
        self.position = mutation.residue_index
        self.structure = structure
        self.model = None
        self.chain = 'A'  # structure.chain get offset will change the chain to A.
        self.code = structure.code
        self.old_chain = str(structure.chain)
        self.sequence = sequence
        self.has_conservation = False
        # -------------- fix for non serverside usage...
        # structure type str: rcsb | swissmodel | homologue | www | local | custom | alphafold2
        if self.structure.type is None or self.structure.type == '':
            self.log.critical(f'Empty structure type: {self.structure}')
            if len(self.code) == 4:
                self.structure.type = 'rcsb'
            else:
                self.structure.type = 'custom'
        # ------------- get coordinates
        if self.structure.coordinates:
            self.coordinates = self.structure.coordinates
        elif self.structure.type == 'rcsb':
            self.coordinates = structure.get_offset_coordinates()
        elif self.structure.type == 'swissmodel':
            # NB. sequence optional argument for ``structure.get_coordinates_w_template_extras``
            # has no a serverside usage
            self.coordinates = structure.get_coordinates_w_template_extras()
        elif self.structure.type == 'alphafold2':  # kept separate for clarity...
            self.coordinates = structure.get_coordinates()
            self.history = {'code': 'AlphaFold2'}
        else:
            self.coordinates = structure.get_coordinates()
        # correct issues.............................
        # these two are very much for the ajax.
        self.chain_definitions = structure.chain_definitions
        # self.chain_definitions seems redundant but str(structure) does not give these.
        self.history = {'code': self.code, 'changes': 'offset and made chain A'}
        # do operations on structure
        self.target_selection = f'(resi {self.position} and chain {self.chain})'
        self.pymol = None
        self._obj_name = 'myprotein'
        with pymol2.PyMOL() as self.pymol:
            self.pymol.cmd.read_pdbstr(self.coordinates, self._obj_name)
            self.N_atoms = self.pymol.cmd.select(self.target_selection)
            # --- verify okay
            if self.N_atoms == 0:
                raise ValueError('Residue is in missing density')
            elif self.N_atoms == 1:
                raise ValueError('Structure is likely an Calpha trace')
            # --- analyses
            # the normal number of heavy atoms (`.normal_HA`) is counting OXT.
            self.has_all_heavy_atoms = self.N_atoms >= self.normal_HA[self.mutation.from_residue] - 1
            self.pymol.cmd.h_add()
            self.neighbours = self.get_neighbours()
            self.SASA = self.get_SASA()
            if self.mutation.from_residue != 'G':
                self.SASA_sidechain = self.get_SASA(f'{self.target_selection} and not name N+H+C+CA+HA+O+OXT')
            else:
                self.SASA_sidechain = self.get_SASA(f'{self.target_selection} and name CA')
            self.RSA = self.SASA / self.maxASA[self.mutation.from_residue]
            self.simbai_ddG = self.get_simbai_ddG()
            self.SS = self.get_SS()
            self.buried = self.RSA <= 0.2
            self.bfactor = self.pymol.cmd.get_model(f'{self.target_selection} and name CA').atom[0].b
            self.ligand_list = self.get_ligand_list()
            t = self.get_distance_to_closest_ligand()
            self.closest_ligand = t['closest']
            self.distance_to_closest_ligand = t['distance']
            if self.structure.type in ('rcsb', 'swissmodel') and not no_conservation:
                self.add_conservation()  # to neighbours
        self.pymol = None

    def get_SS(self, sele=None):
        assert self.pymol is not None, 'Can only be called within a PyMOL session'
        if not sele:
            sele = self.target_selection
        my_dict = {'residues': []}
        CA = f'{sele} and name CA'
        if not self.pymol.cmd.select(CA):
            return 'Unknown'
        else:
            self.pymol.cmd.iterate(CA, "residues.append(ss)", space=my_dict)
            ss = my_dict['residues'][0]
            if ss in self.ss_types:
                return self.ss_types[ss]
            else:
                'Unknown'

    def get_SASA(self, sele=None):
        assert self.pymol is not None, 'Can only be called within a PyMOL session'
        self.pymol.cmd.set('dot_solvent', 1)
        self.pymol.cmd.set('dot_density', 3)
        if not sele:
            sele = self.target_selection
        return self.pymol.cmd.get_area(sele)

    def get_simbai_ddG(self,
                       rsa: Optional[float] = None,
                       from_residue: Optional[str] = None,
                       to_residue: Optional[str] = None,
                       positive_destabilized: Optional[bool] = True) -> float:
        """
        Metric from Caldararu, Blundell, Kepp, J. Chem. Inf. Model. 2021, 61, 4
        https://pubs.acs.org/doi/10.1021/acs.jcim.1c00201

        Note that in the original paper the RSA was calcuated via Naccess.
        (a 1992 software that requires a licence mailed in)
        Here PyMOL is used. In the tests ``pyrosetta.rosetta.core.scoring.sasa.rel_per_res_sc_sasa``
        was also used: faster, but equivalent to PyMOL RSA.

        Do note that statistically Simba-I scores well in benchmarks, but that is mainly driven by
        the fact that the values are always mild.
        Assuming a buried residue is RSA=0. and a surface one is RSA=1.(One can have 1.1 for GXG in a loop),
        one can get the extremes:

        - stabilising surface... LEU to ASP: 2.2 kcal/mol
        - destabilising surface... ASP to LEU: -0.99 kcal/mol
        - stabilising buried... ASP to TRP: 1.3 kcal/mol
        - destabilising buried... TRP to ASP: -4.5 kcal/mol

        Namely: surface[-2.2 to +1], buried [+1.3 to +4.5]

        Also, silent mutations will be off by ±1 kcal/mol.

        :param rsa:
        :param from_residue: single letter code... from Bio.SeqUtils import seq1
        :param to_residue: single letter code
        :param positive_destabilized: ProTherm and Blundell papers have it going up the funnel... This flips it.
        :return: float
        """
        # defaults
        if rsa is None:
            rsa = self.RSA
        if from_residue is None:
            from_residue = self.mutation.from_residue
        if to_residue is None:
            to_residue = self.mutation.to_residue
        # AA constants
        volumes = {'A': 88.6, 'R': 173.4, 'N': 114.1, 'D': 111.1, 'C': 108.5, 'Q': 143.8, 'E': 138.4, 'G': 60.1,
                   'H': 153.2, 'I': 166.7, 'L': 166.7, 'K': 168.6, 'M': 162.9, 'F': 189.9, 'P': 112.7, 'S': 89.0,
                   'T': 116.1, 'W': 227.8, 'Y': 193.6, 'V': 140.0}
        hydrophobicities = {'A': 8.1, 'R': 10.5, 'N': 11.6, 'D': 13.0, 'C': 5.5, 'Q': 10.5, 'E': 12.3, 'G': 9.0,
                            'H': 10.4, 'I': 5.2, 'L': 4.9, 'K': 11.3, 'M': 5.7, 'F': 5.2, 'P': 8.0, 'S': 9.2, 'T': 8.6,
                            'W': 5.4, 'Y': 6.2, 'V': 5.9}
        # ∆∆G = a + b * RSA + c* V + d * H + e * RSA * V + f RSA * H
        offset = -1.6
        rsa_factor = +2.2
        dV_factor = +0.59
        dH_factor = -0.29
        rsa_dV_factor = -0.54
        rsa_dH_factor = +0.49
        # calculate
        dV = (volumes[to_residue] - volumes[from_residue]) / 100
        dH = (hydrophobicities[to_residue] - hydrophobicities[from_residue])
        simba = offset + \
                rsa_factor * rsa + dV_factor * dV + dH_factor * dH + \
                rsa_dV_factor * rsa * dV + rsa_dH_factor * rsa * dH
        if positive_destabilized:
            return -simba
        else:
            return simba

    def get_neighbours(self, threshhold: float = 12, sele: str = None):
        """
        Throughout the modules neighbors/neighbours is accidentally written in British English.
        Do not correct.

        :param threshhold: &Aring;ngstrom distance.
        :return:
        """
        assert self.pymol is not None, 'Can only be called within a PyMOL session'
        if not sele:
            # none of the paretheses are required, but help in legibility...
            sele = f'name CA and byres (({self.target_selection} and name CA) expand {threshhold})'
            # sele = f'(byres ({self.target_selection} around {threshhold})) and name CA'
        my_dict = {'residues': []}
        self.pymol.cmd.iterate(sele, "residues.append({'resi': resi, 'resn': resn, 'chain': chain})", space=my_dict)
        neighbours = my_dict['residues']
        # pymol.cmd.distance gets the last atom distance
        # pymol.cmd.get_distance is atom to atom. 23 s
        # atom.coord gives coord. 0.2 s
        # selector = lambda atom: f'resi {atom.resi} and chain {atom.chain} and name {atom.name}'
        target_coords = [atom.coord for atom in self.pymol.cmd.get_model(self.target_selection).atom]
        for res in neighbours:
            print(res)
            near_coords = [atom.coord for atom in self.pymol.cmd.get_model(self.neigh2selection(res)).atom]
            p_iter = product(map(np.array, target_coords), map(np.array, near_coords))
            res['distance'] = np.min([np.linalg.norm([aa - bb]) for aa, bb in p_iter])
        print(neighbours)
        return neighbours

    def neigh2selection(self, neighbour, name: Optional[str] = None) -> str:
        if name is None:
            return f'resi {neighbour["resi"]} and chain {neighbour["chain"]}'
        else:
            return f'resi {neighbour["resi"]} and chain {neighbour["chain"]} and name {name}'

    def add_conservation(self):
        # get code and chain
        if self.structure.code is None or len(self.structure.code) == 0:
            raise ValueError('No code')
        if self.structure.type == 'swissmodel':  # based upon 5v3j.2.C C
            code, chain = re.search(' (\w{4})\.\d+\.(\w)', self.structure.code).groups()
            # I cannot guarantee self.structure.chain is the same.
        elif len(self.structure.code) == 4:
            code = self.structure.code
            chain = self.old_chain
        else:
            raise ValueError(f'What is {self.structure.code}')
        if not self.error_on_missing_conservation:
            cought = Exception
        else:
            cought = ()
        try:
            con = Consurfer().from_web(code, chain)
            con.apply_offset_by_alignment(self.sequence)
            self.has_conservation = True
            self.pymol.cmd.delete('*')
            self.pymol.cmd.read_pdbstr(self.coordinates, 'mod_')
            con_chain = con.get_consurf_chain()
            if con_chain != 'A':
                con.remap_chains({con_chain: 'A'})
            con.add_bfactor_to_pymol(self.pymol)
            # No need for: self.pymol.cmd.remove('element H')
            self.coordinates = self.pymol.cmd.get_pdbstr()
            self.log.debug(con.data)
        except cought as error:  # () for crappy debug
            self.log.debug(f'{error.__class__.__name__}: {error} at consurf at line ({error.__traceback__.tb_lineno})')
            self.has_conservation = False
            return
            # raise BaseException(str(error))  # no catch
        # {'GLY51:A': {'POS': '1', 'SEQ': '   G', '3LATOM': '   GLY51:A', 'SCORE': ' 1.313', 'COLOR': '  2',
        # 'CONFIDENCE INTERVAL': ' 0.264, 1.652', 'CONFIDENCE COLORS': '    4,1', 'MSA DATA': '  11/300',
        # 'RESIDUE VARIETY': 'A,G,R,V,K,I,E'},
        for neigh_data in self.neighbours:
            try:
                key = con.get_key(neigh_data['resi'])
                neigh_data['variety'] = con.get_variety(key)
                neigh_data['color'] = con.get_color(key)
                neigh_data['conscore'] = con.get_conscore(key)
            except ValueError as error:
                self.log.debug(f'{error.__class__.__name__}: {error} for {neigh_data["resi"]}')
                pass  # absent.

    def get_distance_to_closest_chain(self):
        # this is CA only.
        assert self.pymol is not None, 'Can only be called within a PyMOL session'
        my_dict = {'chains': set(), 'target': {}, 'others': {}}
        self.pymol.cmd.iterate(f'name CA', "chains.add(chain)", space=my_dict)
        if len(my_dict['chains']) == 0:
            return {'closest': None, 'distance': None}
        else:
            self.pymol.cmd.iterate_state(1, f'{self.target_selection} and name CA', "target['CA'] =(x,y,z)",
                                         space=my_dict)
            closest_d = 99999
            closest_c = ''
            for c in my_dict['chains']:
                self.pymol.cmd.iterate_state(1, f'chain {c} and name CA', "others[resi]=(x,y,z)", space=my_dict)
                for a in my_dict['others']:
                    d = self.euclidean(my_dict['target']['CA'], my_dict['others'][a])
                    if d < closest_d:
                        closest_d = d
                        closest_c = f'{a}:{c}'
            return {'closest': closest_c, 'distance': closest_d}

    def get_distance_to_closest_ligand(self):
        """

        :return: {'target': target atom, 'closest': ligand atom, 'distance': Ang distance}
        """
        assert self.pymol is not None, 'Can only be called within a PyMOL session'
        my_dict = {'target': {}, 'ligands': {}}
        # assert self.pymol.cmd.select(self.target_selection) != 0, f'selection {self.target_selection} is invalid'
        self.pymol.cmd.iterate_state(1, f'{self.target_selection}', "target[resi+'.'+name+':'+chain] =(x,y,z)",
                                     space=my_dict)
        if not self.ligand_list:
            return {'target': None, 'closest': None, 'distance': None}
        # self.ligand_list = {'GDP'}
        # my_dict['ligands'] = {'[GDP]4.PA:_': (2.3580000400543213, -11.11400032043457, 49.12099838256836), ...}
        for lig in self.ligand_list:
            # NGL selection
            filler = "ligands['['+resn+']'+resi+'.'+name+':'+chain] = (x,y,z)"
            self.pymol.cmd.iterate_state(1, f'resn {lig}', filler, space=my_dict)
        closest_t = ''
        closest_l = ''
        closest_d = 99999
        for l in my_dict['ligands']:
            for target_name, target_coords in my_dict['target'].items():
                d = self.euclidean(my_dict['ligands'][l], target_coords)
                if d < closest_d:
                    closest_t = target_name
                    closest_l = l
                    closest_d = d
        return {'target': closest_t, 'closest': closest_l, 'distance': closest_d}

    @staticmethod
    def euclidean(a, b):
        return math.sqrt(sum([(a[i] - b[i]) ** 2 for i in range(3)]))

    def get_ligand_list(self):
        assert self.pymol is not None, 'Can only be called within a PyMOL session'
        my_dict = {'residues': set()}
        self.pymol.cmd.iterate(f'all', "residues.add(resn)", space=my_dict)
        exclusion = PyMolTranspiler.boring_ligand + PyMolTranspiler.water_ligand + PyMolTranspiler.aa_ligand
        return {lig.upper() for lig in my_dict['residues'] if lig.upper() not in exclusion}
