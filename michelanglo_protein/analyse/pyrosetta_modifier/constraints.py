import pyrosetta
from .neighbors import MutatorNeighbors

# mro: MutatorBase -> MutatorInit -> MutatorNeighbors -> MutatorCon -> MutatorRelax -> MutatorDescribe -> Mutator
class MutatorCon(MutatorNeighbors):

    def explicit_cons(self, proximity_threshold=3) -> pyrosetta.rosetta.core.scoring.constraints.ConstraintSet:
        """
        Adds contraints between atoms closer than 3 Å between the neighbourhood residues and the other residues.
        Therefore preventing residues (such as bad ligands) from exploding during repacking.
        The constraint generator works with CA only, or without only with peptides,
        hydrogen bonding generator would not capture clashing atom, such as incorrectly protonated ligands.
        This method takes up less than a second or two, it's the relax that takes up more time.
        """
        # sele neigh
        neigh_sele = self._get_neigh_sele(self.radius)
        # get outers
        cc_sele = pyrosetta.rosetta.core.select.residue_selector.CloseContactResidueSelector()
        cc_sele.central_residue_group_selector(neigh_sele)
        cc_sele.threshold(proximity_threshold)
        not_sele = pyrosetta.rosetta.core.select.residue_selector.NotResidueSelector(neigh_sele)
        outer_sele = pyrosetta.rosetta.core.select.residue_selector.AndResidueSelector(cc_sele, not_sele)
        outers = pyrosetta.rosetta.core.select.residue_selector.ResidueVector(outer_sele.apply(self.pose))
        outer_cc_sele = pyrosetta.rosetta.core.select.residue_selector.CloseContactResidueSelector()
        outer_cc_sele.threshold(proximity_threshold)
        # outer_map
        outer_map = {}
        for out in list(outers):
            # get border
            out_sele = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(out)
            outer_cc_sele.central_residue_group_selector(out_sele)
            border_sele = pyrosetta.rosetta.core.select.residue_selector.AndResidueSelector(outer_cc_sele, neigh_sele)
            # border are the border residues of a given outer residue.
            borders = pyrosetta.rosetta.core.select.residue_selector.ResidueVector(border_sele.apply(self.pose))
            outer_map[out] = borders
        # distance
        cons = pyrosetta.rosetta.core.scoring.constraints.ConstraintSet()
        HarmonicFunc = pyrosetta.rosetta.core.scoring.func.HarmonicFunc
        AtomPairConstraint = pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint
        for out, borders in outer_map.items():
            out_res = self.pose.residue(out)
            out_xyzs = {o: out_res.xyz(o) for o in range(1, out_res.natoms() + 1)}
            for border in list(borders):
                b_res = self.pose.residue(border)
                # bs = set(range(1, b_res.natoms() + 1)) - set(b_res.all_bb_atoms())
                bs = list(range(1, b_res.natoms() + 1))
                for b in bs:
                    bxyz = b_res.xyz(b)
                    for o, oxyz in out_xyzs.items():
                        # is sequential bb? skip.
                        if abs(out - border) == 1 and \
                                o in out_res.all_bb_atoms() and \
                                b in b_res.all_bb_atoms():
                            continue
                        d = (bxyz - oxyz).norm()
                        if d < proximity_threshold:
                            b_atom = pyrosetta.AtomID(atomno_in=b, rsd_in=border)
                            o_atom = pyrosetta.AtomID(atomno_in=o, rsd_in=out)
                            harfun = HarmonicFunc(x0_in=d, sd_in=0.2)
                            cons.add_constraint(AtomPairConstraint(b_atom, o_atom, harfun))
        setup = pyrosetta.rosetta.protocols.constraint_movers.ConstraintSetMover()
        setup.constraint_set(cons)
        setup.apply(self.pose)
        return cons
