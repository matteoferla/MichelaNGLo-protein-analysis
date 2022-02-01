from typing import (List, Dict, Optional)

import pyrosetta

from .base import MutatorBase, log, Target

# mro: MutatorBase -> MutatorInit -> MutatorNeighbors -> MutatorCon -> MutatorRelax -> Mutator
class MutatorInit(MutatorBase):
    def __init__(self,
                 pdbblock: str,
                 target_resi: int,
                 target_chain: str = 'A',
                 cycles: int = 2,
                 radius: int = 12,
                 params_filenames: List[str] = (),
                 scorefxn_name: str = 'ref2015',
                 use_pymol_for_neighbours: bool = False,
                 neighbour_only_score: bool = False,
                 outer_constrained: bool = False,
                 remove_ligands: bool = False,
                 single_chain: bool = False,
                 prevent_acceptance_of_incrementor: bool = True):
        """

        Load.

        :param pdbblock: PDB block
        :type pdbblock: str
        :param target_resi: mutate residue PDB index
        :type target_resi: int
        :param target_chain: chain id
        :type target_chain: str
        :param cycles: (opt) cycles of relax.
        :param radius: (opt) angstrom to expand around
        :param params_filenames: list of filenames of params files (rosetta topology files)
        :param scorefxn_name: scorefunction to use, some modifying options are enabled if not 'ref2015'
        :type scorefxn_name: str
        :param use_pymol_for_neighbours: Pyrosetta neighbourhood is calculated by centroid to centroid (C-&beta; generally)
            thus being okay for mutations. However there is the annoyance that ligands may have weirdly defined centroids.
            for example in 5CE4 F128 is over 12Å away from OLE, which it hydrogen bonds with...
        :type use_pymol_for_neighbours: bool
        :param neighbour_only_score: score the neighbourhood only (correctly w/ h-bond correction)
        :type neighbour_only_score: bool
        :param outer_constrained: prevent the loss of interactions with residues beyond the neighbour with constraints.
        :type outer_constrained: bool
        :param remove_ligands: remove ligands?
        :type remove_ligands: bool
        :param single_chain: Remove all bar the first chain...
            assumes there is some biological assembly error. But also good for testing...
        :type single_chain: bool
        :param prevent_acceptance_of_incrementor: there is a peculiarity that in some cases, a score worse
            than the original is accepted. This prevents it... but often fails to so no relaxation is done.
        :type prevent_acceptance_of_incrementor:bool
        """
        _keys = [cycles, radius, params_filenames, scorefxn_name, neighbour_only_score, outer_constrained,
                 remove_ligands, single_chain, prevent_acceptance_of_incrementor]
        log.debug('Initialising')
        log.debug(f'Mutator {_keys}')
        self.scorefxn = self.get_scorefunction(scorefxn_name)
        ap_st = pyrosetta.rosetta.core.scoring.ScoreType.atom_pair_constraint
        self.scorefxn.set_weight(ap_st, 10)
        # correct for split as per https://www.rosettacommons.org/node/11245
        self.scorefxn.weights()  # Create the EnergyMap
        emopts = pyrosetta.rosetta.core.scoring.methods.EnergyMethodOptions(self.scorefxn.energy_method_options())
        emopts.hbond_options().decompose_bb_hb_into_pair_energies(True)
        self.scorefxn.set_energy_method_options(emopts)
        self.scores = {}  # gets filled by .mark()
        self.neighbour_only_score = bool(neighbour_only_score)
        self.cycles = int(cycles)
        self.radius = int(radius)  # radius cannot be float.
        # Load
        self.target = Target(target_resi, target_chain)
        self.pdbblock = pdbblock
        self.params_filenames = list(params_filenames) + self.default_params
        log.debug(self.params_filenames)
        self.pose = self.load_pose(single_chain, remove_ligands)  # self.pose is intended as the damageable version.
        self.native = None
        log.debug('Pose loaded...')
        # Find neighbourhood (pyrosetta.rosetta.utility.vector1_bool)
        if use_pymol_for_neighbours:
            neighbours = self.calculate_neighbours_in_pymol(self.radius)
            self.neighbour_vector = self.targets2vector(neighbours)
        else:
            self.neighbour_vector = self.calculate_neighbours_in_pyrosetta(self.radius)
        log.debug('About to run raw...')
        self.mark('raw')  # mark scores the self.pose
        log.debug(f'Raw scored: {self.scores}')
        # Ready relax
        self.relax = None  # filled by self.ready_relax
        self.movemap = None  # filled by self.ready_relax
        self.ready_relax(self.cycles)
        self.n_preventions = 0
        self.prevent_acceptance_of_incrementor = prevent_acceptance_of_incrementor
        if outer_constrained:
            self.explicit_cons()

    def explicit_cons(self):
        raise NotImplementedError('Filled by constraints.MutatorCon')

    def ready_relax(self, *args, **kwargs):
        raise NotImplementedError('Filled by relax.MutatorRelax')

    def calculate_neighbours_in_pymol(self, *args, **kwargs):
        raise NotImplementedError('Filled by neighbors.MutatorNeighbors and actually very much debug only')

    def calculate_neighbours_in_pyrosetta(self, *args, **kwargs):
        raise NotImplementedError('Filled by neighbors.MutatorNeighbors')

    def targets2vector(self, *args, **kwargs):
        raise NotImplementedError('Filled by neighbors.MutatorNeighbors')

    def load_pose(self, single_chain=False, remove_ligands=False) -> pyrosetta.Pose:
        """
        Loading from str is a bit messy. this simply does that and returns a Pose

        :return: self.pose
        """
        pose = pyrosetta.Pose()
        if self.params_filenames:
            params_paths = pyrosetta.rosetta.utility.vector1_string()  # noqa -- it's there.
            params_paths.extend(self.params_filenames)
            rts = pyrosetta.generate_nonstandard_residue_set(pose, params_paths)
            # the filename AFAIK does nothing.
            pyrosetta.rosetta.core.import_pose.pose_from_pdbstring(pose, self.pdbblock, rts, 'temp.pdb')
        else:
            pyrosetta.rosetta.core.import_pose.pose_from_pdbstring(pose, self.pdbblock)
        if single_chain:
            pose = pose.split_by_chain(1)
        if remove_ligands:
            pyrosetta.rosetta.core.pose.remove_nonprotein_residues(pose)
        return pose

    def mark(self, label: str) -> Dict:
        """
        Save the score to ``.scores`` (inplace)

        :param label: scores is Dict. label is key.
        :return: {ddG: float, scores: Dict[str, float], native:str, mutant:str, rmsd:int}
        """

        ap_st = pyrosetta.rosetta.core.scoring.ScoreType.atom_pair_constraint
        self.scorefxn.set_weight(ap_st, 0)
        if self.neighbour_only_score:
            self.scorefxn(self.pose)
            self.scores[label]: float = (self.scorefxn.get_sub_score(self.pose, self.neighbour_vector) *
                                         self.scaling_factor)
        else:
            self.scores[label]: float = self.scorefxn(self.pose) * self.scaling_factor
        self.scorefxn.set_weight(ap_st, 10)
        return self.scores


    def output_pdbblock(self, pose: Optional[pyrosetta.Pose] = None) -> str:
        """
        This is weird. I did not find the equivalent to ``pose_from_pdbstring``.
        But using buffer works.

        :return: PDBBlock
        """
        if pose is None:
            pose = self.pose
        buffer = pyrosetta.rosetta.std.stringbuf()  # noqa -- it's there.
        pose.dump_pdb(pyrosetta.rosetta.std.ostream(buffer))  # noqa -- it's there.
        return buffer.str()
