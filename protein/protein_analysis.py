__description__ = """
The class ProteinAnalyser builds upon the ProteinLite core and expands it with 
"""

from .core import ProteinCore
from .mutation import Mutation
import re
from Bio.PDB import PDBParser
from Bio.PDB.HSExposure import HSExposureCB
import io, os

class ProteinAnalyser(ProteinCore):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        ### other ###
        self.structure = None
        ### mutation ###
        self._mutation = None
        ## structure
        self.model = None #Structural instance

    ############## elm
    _elmdata = []

    @property
    def elmdata(self):
        if not len(self._elmdata):
            with open(os.path.join(self.settings.reference_folder, 'elm_classes.tsv')) as fh:
                header = ("Accession", "ELMIdentifier", "FunctionalSiteName", "Description", "Regex", "Probability",
                          "#Instances", "#Instances_in_PDB")
                for line in fh:
                    if line[0] == '#':
                        continue
                    if "Accession" in line:
                        continue
                    self._elmdata.append(dict(zip(header, line.replace('"', '').split('\t'))))
            self.__class__._elmdata = self._elmdata  ## change the class attribute too!
            return self._elmdata

    def _set_mutation(self, mutation):
        if isinstance(mutation, str):
            self._mutation = Mutation(mutation)
        else:
            self._mutation = mutation

    mutation = property(lambda self:  self._mutation, _set_mutation)

    # decorator no longer used.
    def _sanitise_position(fun):
        """
        Decorator that makes sure that position is a number. It is a bit unnecassary for a one job task...
        :return: int,
        """
        def sanitiser(self, position):
            if isinstance(position, str):
                position = int(position)
            elif isinstance(position, int):
                pass
            elif not isinstance(position, int):
                position = position.residue_index
            elif position == None:
                position = self.mutation.residue_index
            else:
                position = position
            return fun(self, position)

        return sanitiser


    def _neighbours(self, midresidue, position, marker='*', span=10):
        """
        Gets the 10 AA stretch for mutant or not.
        :param midresidue: what is the letter to put in middle. Used for wt and mutant.
        :param position: number
        :param marker: '*' to surround the midresidue.
        :param span: length of span. default 10.
        :return: 10 aa span.
        """
        halfspan = int(span / 2)
        if position < 5:
            neighbours = '{pre}{m}{i}{m}{post}'.format(pre=self.sequence[:position - 1],
                                                       i=midresidue,
                                                       post=self.sequence[position:position + halfspan],
                                                       m=marker)
        elif position > 5 and len(self.sequence) > position + halfspan:
            neighbours = '{pre}{m}{i}{m}{post}'.format(pre=self.sequence[position - 1 - halfspan:position - 1],
                                                       i=midresidue,
                                                       post=self.sequence[position:position + halfspan],
                                                       m=marker)
        elif len(self.sequence) < position + 5:
            neighbours = '{pre}{m}{i}{m}{post}'.format(pre=self.sequence[position - 1 - halfspan:position - 1],
                                                       i=midresidue,
                                                       post=self.sequence[position:],
                                                       m=marker)
        else:
            neighbours = 'ERROR.'
        return neighbours

    ################### mutant related
    def predict_effect(self):
        assert self.mutation, 'No mutation specified.'
        if self.mutation:
            if not self.check_mutation():
                raise ValueError(self.mutation_discrepancy())
        self.check_elm()
        self.get_features_at_position()

    def check_mutation(self):
        if len(self.sequence) > self.mutation.residue_index and self.sequence[self.mutation.residue_index - 1] == self.mutation.from_residue:
            return True
        else:
            return False  # call mutation_discrepancy to see why.

    def mutation_discrepancy(self):
        # returns a string explaining the `check_mutation` discrepancy error
        neighbours = ''
        if len(self.sequence) < self.mutation.residue_index:
            return 'Uniprot {g} is {l} amino acids long, while user claimed a mutation at {i}.'.format(
                g=self.uniprot,
                i=self.mutation.residue_index,
                l=len(self.sequence)
            )
        else:
            neighbours = self._neighbours(midresidue=self.sequence[self.mutation.residue_index - 1],
                                          position=self.mutation.residue_index,
                                          marker='*')
            return 'Residue {i} is {n} in Uniprot {g}, while user claimed it was {f}. (neighbouring residues: {s}'.format(
                i=self.mutation.residue_index,
                n=self.sequence[self.mutation.residue_index - 1],
                g=self.uniprot,
                f=self.mutation.from_residue,
                s=neighbours
            )

    ################################# ELM

    def _rex_elm(self, neighbours, regex):
        rex = re.search(regex, neighbours)
        if rex:
            return (rex.start(), rex.end())
        else:
            return False

    def check_elm(self):
        assert self.sequence, 'No sequence defined.'
        position = self.mutation.residue_index
        neighbours = self._neighbours(midresidue=self.sequence[position - 1], position=position, span=10, marker='')
        mut_neighbours = self._neighbours(midresidue=self.mutation.to_residue, position=position, span=10, marker='')
        results = []
        elm = self.elmdata
        for r in elm:
            w = self._rex_elm(neighbours, r['Regex'])
            m = self._rex_elm(mut_neighbours, r['Regex'])
            if w != False or m != False:
                match = {'name': r['FunctionalSiteName'],
                         'description': r['Description'],
                         'regex': r['Regex'],
                         'probability': float(r['Probability'])}
                if w != False and m != False:
                    match['x'] = w[0] + position - 5
                    match['y'] = w[1] + position - 5
                    match['status'] = 'kept'
                elif w != False and m == False:
                    match['x'] = w[0] + position - 5
                    match['y'] = w[1] + position - 5
                    match['status'] = 'lost'
                else:
                    match['x'] = m[0] + position - 5
                    match['y'] = m[1] + position - 5
                    match['status'] = 'gained'
                results.append(match)
        self.mutation.elm = sorted(results, key=lambda m: m['probability'] + int(m['status'] == 'kept'))
        return self

    ##################### Position queries.

    def get_features_at_position(self, position=None):
        """
        :param position: mutation, str or position
        :return: list of gNOMAD mutations, which are dictionary e.g. {'id': 'gNOMAD_19_19_rs562294556', 'description': 'R19Q (rs562294556)', 'x': 19, 'y': 19, 'impact': 'MODERATE'}
        """
        position = position if position is not None else self.mutation.residue_index
        return self.get_features_near_position(position, wobble=0)

    def get_features_near_position(self, position=None, wobble=10):
        position = position if position is not None else self.mutation.residue_index
        valid = [{**f, 'type': g} for g in self.features for f in self.features[g] if f['x'] - wobble < position < f['y'] + wobble]
        svalid = sorted(valid, key=lambda v: int(v['y']) - int(v['x']))
        return svalid

    def get_gNOMAD_near_position(self, position=None, wobble=5):
        """
        :param position: mutation, str or position
        :param wobble: int, number of residues before and after.
        :return: list of gNOMAD mutations, which are dictionary e.g. {'id': 'gNOMAD_19_19_rs562294556', 'description': 'R19Q (rs562294556)', 'x': 19, 'y': 19, 'impact': 'MODERATE'}
        """
        position = position if position is not None else self.mutation.residue_index
        valid = [g for g in self.gNOMAD if g.x - wobble < position < g.y + wobble]
        svalid = sorted(valid, key=lambda v: v.y - v.x)
        return svalid

    #### THE FUTURE

    def analyse_structure(self, position=None):
        position = position if position is not None else self.mutation.residue_index
        structures = self._get_structures_with_position(position)
        if not structures:
            return self
        self.structure_code = structures[0]['id']
        self.model = StructureAnalyser(position, structure=self.get_structure(self, structures[0]['id']), chain=structures[0]['chain'], code=structures[0]['id'])
        self._get_structure_neighbours(self.structure, position)


    def _get_structures_with_position(self, position):
        """
        Fetches structures that exists at a given position.
        :param position: mutation, str or position
        :return: list of self.pdbs+self.swissmodel+self.pdb_matches...
        """
        return [pdb for pdb in self.pdbs + self.swissmodel + self.pdb_matches if int(pdb['x']) < position < int(pdb['y'])]


    def get_structure(self, code):  #a function to fetch pdb?
        pass


    # conservation score
    # disorder

class StructureAnalyser:

    def __init__(self, position, structure, chain, code):
        self.position = position
        self.structure = PDBParser().get_structure('model', io.StringIO(structure))
        self.chain = chain
        self.code = code
        self.target_residue = self.structure[0][self.chain][self.position]
        self._target_hse = None

    @property
    def target_hse(self):
        if not self._target_hse:
            hse = HSExposureCB(self.structure)
            self._target_hse = hse[(self.structure[0][self.chain].get_id(), self.target_residue.id)]
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
        overthreshhold = 10+threshhold #no atom in a arginine will be within threshhold of a given atom if any atom is greater than this.
        double_overthreshhold = 20 + threshhold #no atom in a arginine will be within threshhold of a any atom in a given residue if any atom is greater than this.
        doublebreak = False
        for residue in self.structure[0][self.chain]:
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


