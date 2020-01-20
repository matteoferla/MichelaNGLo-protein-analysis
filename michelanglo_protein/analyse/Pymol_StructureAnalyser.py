from ..structure import Structure
from ..mutation import Mutation
from michelanglo_transpiler import pymol2, PyMolTranspiler
import math

class StructureAnalyser:
    """
    3D details of a mutation based on the structure which is stored in ``.structure`` (which is a ``Bio.PDB.PDBParser`` obj.)
    """
    #Wikipedia Tien et al. 2013 (emp.)
    maxASA = {'A': 121.0, 'R': 265.0, 'N': 187.0, 'D': 187.0, 'C': 148.0, 'E': 214.0, 'Q': 214.0, 'G': 97.0, 'H': 216.0, 'I': 195.0, 'L': 191.0, 'K': 230.0, 'M': 203.0, 'F': 228.0, 'P': 154.0, 'S': 143.0, 'T': 163.0, 'W': 264.0, 'Y': 255.0, 'V': 165.0}
    # normal heavy atom count done with the script genenerate_AA_MCS.py
    normal_HA = {'A': 6, 'R': 12, 'N': 9, 'D': 9, 'C': 7, 'E': 10, 'Q': 10, 'G': 5, 'H': 11, 'I': 9, 'L': 9, 'K': 10, 'M': 9, 'F': 12, 'P': 8, 'S': 7, 'T': 8, 'W': 15, 'Y': 13, 'V': 8}
    # I think PyMOL has H,S,L only?
    ss_types = {'H': 'Helix','S': 'Sheet', 'L': 'Loop', 'G': '3_10 helix', 'I': 'Pi helix', 'T': 'Turn', 'C': 'Coil', 'E': 'Sheet', 'B': 'Beta bridge', '-': 'Unassigned'}

    def __init__(self, structure: Structure, mutation: Mutation):
        """

        :param structure: a instance of Structure, a former namedtuple and is in core.py
        :param mutation: a instance of mutation.
        """
        self.mutation = mutation
        self.position = mutation.residue_index
        self.structure = structure
        self.model = None
        self.chain = 'A' #structure.chain get offset will change the chain to A.
        self.coordinates = structure.get_offset_coordinates()
        self.code = structure.code
        self.target_selection = f'(resi {self.position} and chain {self.chain})'
        self.pymol = None
        self._obj_name = 'myprotein'
        with pymol2.PyMOL() as self.pymol:
            self.pymol.cmd.read_pdbstr(self.coordinates, self._obj_name)
            self.N_atoms = self.pymol.cmd.select(self.target_selection)
            self.has_all_heavy_atoms = self.N_atoms >= self.normal_HA[self.mutation.from_residue]
            self.pymol.cmd.h_add()
            self.neighbours = self.get_neighbours()
            self.SASA = self.get_SASA()
            if self.mutation.from_residue != 'G':
                self.SASA_sidechain = self.get_SASA(f'{self.target_selection} and not name N+H+C+CA+HA+O+OXT')
            else:
                self.SASA_sidechain = self.get_SASA(f'{self.target_selection} and name CA')
            self.RSA = self.SASA / self.maxASA[self.mutation.from_residue]
            self.SS = self.get_SS()
            self.buried = self.RSA >= 0.2
            self.ligand_list = self.get_ligand_list()
            t = self.get_distance_to_closest_ligand()
            self.closest_ligand = t['closest']
            self.distance_to_closest_ligand = t['distance']
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

    def get_neighbours(self, threshhold: float = 3, sele: str = None):
        """

        :param threshhold: &Aring;ngstrom distance.
        :return:
        """
        assert self.pymol is not None, 'Can only be called within a PyMOL session'
        if not sele:
            sele = f'(not {self.target_selection}) and (byres ({self.target_selection} around {threshhold})) and name CA'
        my_dict = {'residues': []}
        self.pymol.cmd.iterate(sele, "residues.append((resi,resn))", space=my_dict)
        return my_dict['residues']


    def get_distance_to_closest_chain(self):
        ## this is CA only.
        assert self.pymol is not None, 'Can only be called within a PyMOL session'
        my_dict = {'chains': set(), 'target': {}, 'others': {}}
        self.pymol.cmd.iterate(f'name CA', "chains.add(chain)", space=my_dict)
        if len(my_dict['chains']) == 0:
            return {'closest': None, 'distance': None}
        else:
            self.pymol.cmd.iterate_state(1, f'{self.target_selection} and name CA', "target['CA'] =(x,y,z)", space=my_dict)
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
        self.pymol.cmd.iterate_state(1, f'{self.target_selection}', "target[resi+'.'+name+':'+chain] =(x,y,z)", space=my_dict)
        if not self.ligand_list:
            return {'target': None, 'closest': None, 'distance': None}
        for lig in self.ligand_list:
            self.pymol.cmd.iterate_state(1, f'resn {lig}', "ligands['['+resn+']'+resi+'.'+name+':'+chain] = (x,y,z)", space=my_dict)
        closest_t = ''
        closest_l = ''
        closest_d = 99999
        for l in my_dict['ligands']:
            for t in my_dict['target']:
                d = self.euclidean(my_dict['ligands'][l], my_dict['target'][t])
                if d < closest_d:
                    closest_t = t
                    closest_l = l
                    closest_d = d
        return {'target': closest_t, 'closest': closest_l, 'distance': closest_d}


    @staticmethod
    def euclidean(a, b):
        return math.sqrt(sum([(a[i]-b[i])**2 for i in range(3)]))

    def get_ligand_list(self):
        assert self.pymol is not None, 'Can only be called within a PyMOL session'
        my_dict = {'residues': set()}
        self.pymol.cmd.iterate(f'all', "residues.add(resn)", space=my_dict)
        exclusion = PyMolTranspiler.boring_ligand + PyMolTranspiler.water_ligand + PyMolTranspiler.aa_ligand
        return {lig.upper() for lig in my_dict['residues'] if lig.upper() not in exclusion}

