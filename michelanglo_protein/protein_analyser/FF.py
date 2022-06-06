# pyrosetta can throw segfaults:
from typing import Union, Dict, Optional, Callable, List

from .base import ProteinAnalyserBase
# pyrosetta can throw segfaults:
from .subprocess import run_subprocess
from ..analyse import Mutator
from ..mutation import Mutation
from ..gnomad_variant import Variant


class ProteinAnalyserFF(ProteinAnalyserBase):

    def analyse_FF(self, spit_process=True,
                   scaling_factor: Optional[float] = None,
                   **mutator_options) -> Union[Dict, None]:
        """
        Calls the pyrosetta, which tends to raise segfaults, hence the whole subpro business.

        :param spit_process: run as a separate process to avoid segfaults?
        :params scaling_factor: multiplied to fix overestimated ddG
        :params mutator_options: neighbour_only_score, outer_constrained for debug
        :return:
        """
        if scaling_factor:
            Mutator.scaling_factor = scaling_factor

        if self.pdbblock is None:
            # to do remember what kind of logging happens down here...
            return {'error': 'ValueError', 'msg': 'no PDB block'}
        ### prepare.
        init_settings = {**self._init_settings, **mutator_options}

        def analysis(to_resn, init_settings):
            mut = Mutator(**init_settings)
            return {**mut.analyse_mutation(to_resn),
                    'scaling_factor': mut.scaling_factor,
                    'neighbor_description': mut.describe_neighbors()}

        if not spit_process:
            msg = analysis(init_settings=init_settings, to_resn=self.mutation.to_residue)
        else:
            msg = run_subprocess(analysis, to_resn=self.mutation.to_residue, init_settings=init_settings)
        self.energetics = msg
        return msg

    def analyse_gnomad_FF(self, spit_process=True,
                          scaling_factor: Optional[float] = None,
                          **mutator_options) -> Union[Dict, None]:
        """
        Calls the pyrosetta, which tends to raise segfaults, hence the whole subpro business.

        :param spit_process: run as a separate process to avoid segfaults?
        :params scaling_factor: multiplied to fix overestimated ddG
        :params mutator_options: neighbour_only_score, outer_constrained for debug
        :return:
        """
        if scaling_factor:
            Mutator.scaling_factor = scaling_factor
        if self.pdbblock is None:
            return {'error': 'ValueError', 'msg': 'no PDB block'}
        ### perpare.
        init_settings = {**self._init_settings, **mutator_options}

        def analysis(gnomads, init_settings):
            mut = Mutator(**init_settings)
            return mut.score_gnomads(gnomads)

        close_variants: List[Variant] = []
        for neigh in self.structural.neighbours:
            close_variants.extend([var for var in self.clinvar if var.residue_index == neigh['resi']])
            close_variants.extend([var for var in self.gnomAD if var.residue_index == neigh['resi']])
        if not spit_process:
            msg = analysis(init_settings=init_settings, gnomads=close_variants)
        else:
            msg = run_subprocess(analysis, init_settings=init_settings, gnomads=close_variants)
        self.energetics_gnomAD = msg
        return msg

    def analyse_other_FF(self,
                         mutation: Union[Mutation, str],
                         algorithm,
                         spit_process=True,
                         scaling_factor: Optional[float] = None,
                         **mutator_options) -> Union[Dict, None]:
        if scaling_factor:
            Mutator.scaling_factor = scaling_factor
        # sort out mutation
        if isinstance(mutation, str):
            mutation = Mutation(mutation)
        elif isinstance(mutation, Mutation):
            pass
        else:
            raise TypeError(f'whats {mutation}?')
        if self.pdbblock is None:
            return {'error': 'ValueError', 'msg': 'no PDB block'}
        ### perpare.
        init_settings = {**self._init_settings, **mutator_options}
        init_settings['target_resi'] = mutation.residue_index

        def relax(resi, from_resn, to_resn, init_settings):
            mut = Mutator(**init_settings)  # altered target_residue from that of the mutation!
            results = mut.analyse_mutation(to_resn)
            return {'coordinates': results['mutant'], 'ddg': results['ddG']}

        def repack(resi, from_resn, to_resn, init_settings):
            mut = Mutator(**init_settings)
            return mut.repack_other(resi, from_resn, to_resn)

        analysis: Callable
        if algorithm == 'relax':
            analysis = relax
        elif algorithm == 'repack':
            analysis = repack
        else:
            ValueError(f'What is this {algorithm}')
        if not spit_process:
            msg = analysis(init_settings=init_settings,
                           to_resn=mutation.to_residue,
                           from_resn=mutation.from_residue,
                           resi=mutation.residue_index)
        else:

            msg = run_subprocess(analysis,
                                 init_settings=init_settings,
                                 to_resn=mutation.to_residue,
                                 from_resn=mutation.from_residue,
                                 resi=mutation.residue_index
                                 )
        return msg

    def phosphorylate_FF(self, spit_process=True) -> Union[str, None]:
        """
                Calls the pyrosetta, which tends to raise segfaults, hence the whole subpro business.

                :param spit_process: run as a separate process to avoid segfaults?
                :return: a PDB block not a score dict!
                """
        if self.pdbblock is None:
            # print('no self.pdbblock')
            return None  # it returns a string not a dictionary.
        elif 'PSP_modified_residues' not in self.features:
            # print('no features')
            return None
        elif not self.features['PSP_modified_residues']:
            # print('no features2')
            return None
        ### perpare.
        init_settings = self._init_settings

        def analysis(ptms, init_settings):
            mut = Mutator(**init_settings)
            return mut.make_phospho(ptms)

        if not spit_process:
            msg = analysis(init_settings=init_settings, ptms=self.features['PSP_modified_residues'])
        else:

            msg = run_subprocess(analysis, init_settings=init_settings, ptms=self.features['PSP_modified_residues'])
        # self.phosphorylated_pdbblcok = msg
        return msg
