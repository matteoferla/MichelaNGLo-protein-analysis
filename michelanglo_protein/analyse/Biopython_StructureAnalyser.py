"""
This module uses Biopython for the calculation of HSExpose

"""

class StructureAnalyser:
    """
    3D details of a mutation based on the structure which is stored in ``.structure`` (which is a ``Bio.PDB.PDBParser`` obj.)
    """
    def __init__(self, structure, mutation):
        """

        :param structure: a instance of Structure, a former namedtuple and is in core.py
        :param mutation: a instance of mutation.
        """
        self.mutation = mutation
        self.position = mutation.residue_index
        self.structure = structure
        self.model = PDBParser(quiet=False).get_structure('model', io.StringIO(structure.get_offset_coordinates()))
        print('hello')
        self.chain = structure.chain
        self.code = structure.code
        self.target_residue = self.model[0][self.chain][self.position]
        self._target_hse = None

    @property
    def target_hse(self):
        if not self._target_hse:
            hse = HSExposureCB(self.model)
            self._target_hse = hse[(self.model[0][self.chain].get_id(), self.target_residue.id)]
        return self._target_hse


    def is_full_surface(self):
        """
        Uses half solvent exposure within PDB module.
        """
        return self.target_hse[0] < 1

    def is_core(self):
        """
        Uses half solvent exposure within PDB module.
        """
        return self.target_hse[0] > 15

    def is_partial_surface(self):
        """
        Uses half solvent exposure within PDB module.
        """
        return 0 < self.target_hse[0] < 15

    def get_superficiality(self):
        if self.target_hse[0] == 0:
            return 'surface'
        elif self.target_hse[0] < 15:
            return 'partially buried'
        else:
            return 'buried'

    def get_structure_neighbours(self, threshhold:float=3):
        """

        :param threshhold: &Aring;ngstrom distance.
        :return:
        """
        # the longest amino acid is 10&Aring; long.
        neighbours = set()
        overthreshhold = 10+threshhold # no atom, even within a arginine, will be within threshhold of a given atom if any atom is greater than this.
        double_overthreshhold = 20 + threshhold #no atom, even within a arginine, will be within threshhold of a any atom in a given residue if any atom is greater than this.
        doublebreak = False
        for residue in self.model[0][self.chain]:
            for atom in residue:
                if doublebreak==True:
                    doublebreak = False
                    break
                for ref_atoms in self.target_residue:
                    if atom - ref_atoms > double_overthreshhold:
                        doublebreak=True # no point in checking the distances to the whole residue
                        break
                    elif atom - ref_atoms > overthreshhold:
                        break # no point in checking the distances to this atom
                    elif atom - ref_atoms < 4:
                        if residue.id[0] == 'W':
                            break #water
                        neighbours.add(residue.id[1])
        return neighbours