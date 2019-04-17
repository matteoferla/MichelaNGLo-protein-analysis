__description__ = """
The class ProteinAnalyser builds upon the ProteinLite core and expands it with 
"""

from .core import ProteinCore
import re

class ProteinAnalyser(ProteinCore):

    # decorator
    def _sanitise_position(fun):
        """
        Decorator that makes sure that position is a number.
        :return: int,
        """

        def sanitiser(self, position):
            if isinstance(position, str):
                position = int(position)
            elif not isinstance(position, int):
                position = position.residue_index
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
    def predict_effect(self, mutation=None):
        if mutation:
            self.mutation = mutation
        assert self.mutation, 'No mutation specified.'
        self.check_elm()

    def check_mutation(self, mutation):
        if len(self.sequence) > mutation.residue_index and self.sequence[mutation.residue_index - 1] == mutation.from_residue:
            return True
        else:
            return False  # call mutation_discrepancy to see why.

    def mutation_discrepancy(self, mutation):
        # returns a string explaining the `check_mutation` discrepancy error
        neighbours = ''
        if len(self.sequence) < mutation.residue_index:
            return 'Uniprot {g} is {l} amino acids long, while user claimed a mutation at {i}.'.format(
                g=self.uniprot,
                i=mutation.residue_index,
                l=len(self.sequence)
            )
        else:
            neighbours = self._neighbours(midresidue=self.sequence[mutation.residue_index - 1],
                                          position=mutation.residue_index,
                                          marker='*')
            return 'Residue {i} is {n} in Uniprot {g}, while user claimed it was {f}. (neighbouring residues: {s}'.format(
                i=mutation.residue_index,
                n=self.sequence[mutation.residue_index - 1],
                g=self.uniprot,
                f=mutation.from_residue,
                s=neighbours
            )

    ################################# ELM

    def _rex_elm(self, neighbours, regex):
        rex = re.search(regex, neighbours)
        if rex:
            return (rex.start(), rex.end())
        else:
            return False

    def check_elm(self, mutation=None):
        if mutation:
            self.mutation = mutation
        assert self.sequence, 'No sequence defined.'
        position = self.mutation.residue_index
        neighbours = self._neighbours(midresidue=self.sequence[position - 1], position=position, span=10, marker='')
        mut_neighbours = self._neighbours(midresidue=self.mutation.to_residue, position=position, span=10, marker='')
        results = []
        for r in self.elmdata:
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

    @_sanitise_position
    def get_features_at_position(self, position):
        """
        :param position: mutation, str or position
        :return: list of gNOMAD mutations, which are dictionary e.g. {'id': 'gNOMAD_19_19_rs562294556', 'description': 'R19Q (rs562294556)', 'x': 19, 'y': 19, 'impact': 'MODERATE'}
        """
        return self.get_features_near_position(position, wobble=0)

    @_sanitise_position
    def get_features_near_position(self, position, wobble=10):
        valid = [{**f, 'type': g} for g in self.features for f in self.features[g] if f['x'] - wobble < position < f['y'] + wobble]
        svalid = sorted(valid, key=lambda v: int(v['y']) - int(v['x']))
        return svalid

    @_sanitise_position
    def get_gNOMAD_near_position(self, position, wobble=5):
        """
        :param position: mutation, str or position
        :param wobble: int, number of residues before and after.
        :return: list of gNOMAD mutations, which are dictionary e.g. {'id': 'gNOMAD_19_19_rs562294556', 'description': 'R19Q (rs562294556)', 'x': 19, 'y': 19, 'impact': 'MODERATE'}
        """

        valid = [g for g in self.gNOMAD if g['x'] - wobble < position < g['y'] + wobble]
        svalid = sorted(valid, key=lambda v: int(v['y']) - int(v['x']))
        return svalid

    #### THE FUTURE

    @_sanitise_position
    def get_structures_with_position(self, position):
        """
        Fetches structures that exists at a given position.
        :param position: mutation, str or position
        :return: list of self.pdbs+self.swissmodel+self.pdb_matches...
        """
        return [pdb for pdb in self.pdbs + self.swissmodel + self.pdb_matches if int(pdb['x']) < position < int(pdb['y'])]

    def is_core(self):
        pass

    @_sanitise_position
    def get_structure_neighbours(self, position):
        """
        Finds residues nearby
        :param position:
        :return:
        """
        pdbs = self.get_structures_with_position(position)
        if not pdbs:
            return None
        pdb = pdbs[0]
        #a function to fetch pdb?


    # conservation score
    # disorder
