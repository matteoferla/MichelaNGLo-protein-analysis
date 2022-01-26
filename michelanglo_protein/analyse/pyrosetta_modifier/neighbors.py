import pyrosetta
import re
from .initialize import MutatorInit, log, Target
from typing import (List)


# mro: MutatorBase -> MutatorInit -> MutatorNeighbors -> MutatorCon -> MutatorRelax -> Mutator
class MutatorNeighbors(MutatorInit):

    def calculate_neighbours_in_pyrosetta(self, radius: int = 12) -> pyrosetta.rosetta.utility.vector1_bool:  # noqa
        """
        Gets the residues within the radius of target. THis method uses pyrosetta.
        It is for filling self.neighbour_vector

        :return: self.neighbour_vector
        :rtype: pyrosetta.rosetta.utility.vector1_bool
        """
        neigh_sele = self._get_neigh_sele(radius)
        return neigh_sele.apply(self.pose)

    def _get_neigh_sele(self, radius):
        r = self.target_pdb2pose(self.target)
        resi_sele = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(r)
        return pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector(resi_sele, radius, True)

    def calculate_neighbours_in_pymol(self, radius: int = 4) -> List[Target]:
        """
        DEBUG ONLY. TODO switch to:

        cc_sele = pyrosetta.rosetta.core.select.residue_selector.CloseContactResidueSelector()
        cc_sele.central_residue_group_selector(resi_sele)
        cc_sele.threshold(3)

        Gets the residues within the radius of target. THis method uses PyMOL!
        It is for filling self.neighbour_vector, but via self.targets2vector()

        :return: the targets
        :rtype: List[Target]
        """
        import pymol2
        with pymol2.PyMOL() as pymol:
            pymol.cmd.read_pdbstr(self.pdbblock, 'blockprotein')
            sele = f"name CA and (byres chain {self.target.chain} and resi {self.target.resi} around {radius})"
            atoms = pymol.cmd.get_model(sele)
            neighbours = []
            for atom in atoms.atom:
                res = int(re.match(r'\d+', atom.resi).group())
                neighbours.append(Target(resi=res, chain=atom.chain))
        return neighbours

    def targets2vector(self, targets: List[Target]) -> pyrosetta.rosetta.utility.vector1_bool:  # noqa -- it's there.
        neighbours = pyrosetta.rosetta.utility.vector1_bool(self.pose.total_residue())  # noqa -- it's there.
        for target in targets:
            r = self.target_pdb2pose(target)
            neighbours[r] = True
        return neighbours

    def target_pdb2pose(self, target: Target) -> int:
        return self._pdb2pose(chain=target.chain, res=target.resi)

    @property
    def _pdb2pose(self):
        return self.pose.pdb_info().pdb2pose

    def get_pdb_neighbours(self):
        neighs = pyrosetta.rosetta.core.select.residue_selector.ResidueVector(self.neighbour_vector)
        pose2pdb = self.pose.pdb_info().pose2pdb
        return [pose2pdb(r) for r in list(neighs)]

    def is_ligand_in_sele(self):
        lig_sele = pyrosetta.rosetta.core.select.residue_selector.ResiduePropertySelector(
            pyrosetta.rosetta.core.chemical.ResidueProperty.LIGAND)
        lig_resi = pyrosetta.rosetta.core.select.residue_selector.ResidueVector(lig_sele.apply(self.pose))
        neigh_resi = pyrosetta.rosetta.core.select.residue_selector.ResidueVector(self.neighbour_vector)
        if set(neigh_resi) & set(lig_resi):
            return True
        else:
            return False
