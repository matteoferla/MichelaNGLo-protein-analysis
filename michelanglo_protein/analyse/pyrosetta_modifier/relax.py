from typing import (Optional, Tuple, Union, Dict, List)

import pyrosetta
from Bio.SeqUtils import seq3

from .base import log
from .constraints import MutatorCon
from ...gnomad_variant import Variant

# mro: MutatorBase -> MutatorInit -> MutatorNeighbors -> MutatorCon -> MutatorRelax -> Mutator
class MutatorRelax(MutatorCon):

    def ready_relax(self, cycles: int = 1, fixed_bb: bool = False) -> pyrosetta.rosetta.protocols.moves.Mover:  # noqa
        """

        :param cycles:
        :return:
        """
        if cycles == 0:
            # relax actually accepts zero as a cycle number. It goes on for eternity.
            self.relax = pyrosetta.rosetta.protocols.moves.NullMover()  # noqa -- it's there.
            return self.relax
        self.relax = pyrosetta.rosetta.protocols.relax.FastRelax(self.scorefxn, cycles)  # noqa -- it's there.
        self.movemap = pyrosetta.MoveMap()
        if not fixed_bb:
            self.movemap.set_bb(self.neighbour_vector)
        else:
            self.movemap.set_bb(False)
        self.movemap.set_chi(self.neighbour_vector)
        # jump is false by def.
        self.relax.set_movemap(self.movemap)
        if self.scorefxn.get_weight(pyrosetta.rosetta.core.scoring.ScoreType.cart_bonded) > 0:
            # it's cartesian!
            self.relax.cartesian(True)
            self.relax.minimize_bond_angles(True)
            self.relax.minimize_bond_lengths(True)
        if hasattr(self.relax, 'set_movemap_disables_packing_of_fixed_chi_positions'):
            # this is relatively new
            self.relax.set_movemap_disables_packing_of_fixed_chi_positions(True)
        else:
            print("UPDATE YOUR PYROSETTA NOW.")
        return self.relax

    def mutate(self, aa):
        res = self.target_pdb2pose(self.target)
        if res == 0:
            raise ValueError('Residue not in structure')
        pyrosetta.toolbox.mutate_residue(self.pose, res, aa)

    def do_relax(self):
        """
        Does relax but prevents the odd case that an exploded structure is accepted.
        """
        if self.n_preventions > 3:
            log.warning(f'Relax cycle failed too many times!')
            return  # do nothing.
        if self.prevent_acceptance_of_incrementor:
            original = self.pose.clone()
            initial = self.scorefxn(self.pose)
            self.relax.apply(self.pose)
            final = self.scorefxn(self.pose)
            if initial + 0.1 < final:
                log.warning(f'Relax cycle failed {initial} < {final} for {self.target}')
                self.pose = original
                self.n_preventions += 1
                self.do_relax()
            else:
                self.n_preventions = 0
        else:
            self.relax.apply(self.pose)

    def _repack_gnomad(self, pose_idx, from_resi, to_resi) -> int:
        self.pose = self.native.clone()
        self.pose.remove_constraints()
        # local repack...
        assert from_resi in self.aminoacid_letters, f'What is a {from_resi} AA?'
        assert to_resi in self.aminoacid_letters, f'What is a {to_resi} AA?'
        pyrosetta.toolbox.mutate_residue(self.pose,
                                         mutant_position=pose_idx,
                                         mutant_aa=from_resi,  # silent mutation first
                                         pack_radius=7.0,
                                         pack_scorefxn=self.scorefxn)
        resi_sele = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(pose_idx)
        neigh_sele = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector(resi_sele, 7, True)
        self.neighbour_vector = neigh_sele.apply(self.pose)
        if self.neighbour_only_score:
            ref = self.scorefxn.get_sub_score(self.pose, self.neighbour_vector)
        else:
            ref = self.scorefxn(self.pose)
        pyrosetta.toolbox.mutate_residue(self.pose,
                                         mutant_position=pose_idx,
                                         mutant_aa=to_resi,
                                         pack_radius=7.0,
                                         pack_scorefxn=self.scorefxn)
        if self.neighbour_only_score:
            mut = self.scorefxn.get_sub_score(self.pose, self.neighbour_vector)
        else:
            mut = self.scorefxn(self.pose)
        return mut - ref

    def repack_other(self, residue_index, from_residue, to_residue):
        self.native = self.pose.clone()
        pose2pdb = self.native.pdb_info().pdb2pose
        pose_idx = pose2pdb(chain='A', res=residue_index)
        return {'ddg':         self._repack_gnomad(pose_idx, from_residue, to_residue),
                'coordinates': self.output_pdbblock(self.pose)}

    def score_gnomads(self, gnomads: List[Variant]):
        """
        This is even more sloppy than the mutant scoring. FastRelax is too slow for a big protein.
        Repacking globally returns a subpar score. Hence why each step has its own repack...

        :param gnomads: list of gnomads.
        :return:
        """
        # self.repack(self.pose)
        self.native = self.pose.clone()
        self.mark('wt')
        pose2pdb = self.native.pdb_info().pdb2pose
        ddG = {}
        for record in gnomads:  # type: Variant
            if record.type == 'nonsense':
                continue
            elif record.to_residue in ('X', '=', '*', 'fs'):
                # this is just to be clear for me reading the code that it did not come from here
                # as the next one catches it too :zany_face:
                continue
            elif record.to_residue not in self.aminoacid_letters:
                # lol.
                continue
            elif record.from_residue not in self.aminoacid_letters:
                continue
            n = pose2pdb(chain='A', res=record.residue_index)
            if n == 0:
                continue
            if record.mutation in ddG:
                # print('duplicate mutation.')
                continue
            ddG[record.mutation] = self._repack_gnomad(n, record.from_residue, record.to_residue)
        return ddG

    def make_phospho(self, ptms):
        phospho = self.pose.clone()
        MutateResidue = pyrosetta.rosetta.protocols.simple_moves.MutateResidue  # noqa -- it's there.
        pdb2pose = phospho.pdb_info().pdb2pose
        changes = 0
        resi_sele = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
        for record in ptms:
            if record['ptm'] == 'ub':
                continue  # What is a proxy for ubiquitination??
            elif record['ptm'] == 'p':
                patch = 'phosphorylated'
            elif record['ptm'] == 'ac':
                patch = 'acetylated'
            elif record['ptm'] == 'm1' and record['from_residue'].upper() == 'LYS':
                # monomethylarginine (NMM) will segfault
                patch = 'monomethylated'
            elif record['ptm'] == 'm2' and record['from_residue'].upper() == 'LYS':
                # dimethylarginine (DA2) will segfault
                patch = 'dimethylated'
            elif record['ptm'] == 'm3' and record['from_residue'].upper() == 'LYS':
                # is trimethylarginine a thing?
                patch = 'trimethylated'
            else:
                continue  # no Gal
                # raise ValueError(f'What is {record["ptm"]}?')
            new_res = f"{seq3(record['from_residue']).upper()}:{patch}"
            r = pdb2pose(res=int(record['residue_index']), chain='A')
            if r == 0:  # missing density.
                continue
            MutateResidue(target=r, new_res=new_res).apply(phospho)
            resi_sele.append_index(r)
        neigh_sele = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector(resi_sele, 7, True)
        self.neighbour_vector = neigh_sele.apply(self.pose)
        self.ready_relax(cycles=1)
        self.do_relax()
        return self.output_pdbblock(phospho)

    def repack(self, target: Optional[pyrosetta.rosetta.core.pose.Pose] = None) -> None:
        """
        This actually seems to make stuff worse on big protein.
        Not as good as ``pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn)``
        But that is slower and less good than FastRelax...

        :param target:
        :return:
        """
        if target is None:
            target = self.pose
        packer_task = pyrosetta.rosetta.core.pack.task.TaskFactory.create_packer_task(target)
        packer_task.restrict_to_repacking()
        pyrosetta.rosetta.core.pack.pack_rotamers(target, self.scorefxn, packer_task)

    def get_diff_solubility(self) -> float:
        """
        Gets the difference in solubility (fa_sol) for the protein.
        fa_sol = Gaussian exclusion implicit solvation energy between protein atoms in different residue

        :return: fa_sol kcal/mol
        """
        _, _, diff = self.get_term_scores(pyrosetta.rosetta.core.scoring.ScoreType.fa_sol)
        return diff

    def get_all_scores(self) -> Dict[str, Dict[str, Union[float, str]]]:
        data: Dict[str, Dict[str, Union[float, str]]] = {}
        for term in self.scorefxn.get_nonzero_weighted_scoretypes():
            data[term.name] = dict(zip(['native', 'mutant', 'difference'], self.get_term_scores(term)))
            data[term.name]['weight'] = self.scorefxn.get_weight(term)
            data[term.name]['meaning']: str = self.term_meanings[term.name]
        return data

    def get_term_scores(self, term: pyrosetta.rosetta.core.scoring.ScoreType) -> Tuple[float, float, float]:
        n = self.scorefxn.score_by_scoretype(self.native, term) * self.scaling_factor
        m = self.scorefxn.score_by_scoretype(self.pose, term) * self.scaling_factor
        d = m - n
        return n, m, d

    def get_diff_res_score(self) -> float:
        """
        Gets the difference in score for that residue

        :return: per_residue kcal/mol
        """
        # segfaults if score is not run globally first!
        i = self.target_pdb2pose(self.target)
        r = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(i)
        n = self.scorefxn.get_sub_score(self.native, r.apply(self.native)) * self.scaling_factor
        m = self.scorefxn.get_sub_score(self.pose, r.apply(self.native)) * self.scaling_factor
        return m - n

    def get_res_score_terms(self, pose) -> dict:
        data = pose.energies().residue_total_energies_array()  # structured numpy array
        i = self.target_pdb2pose(self.target) - 1  # pose numbering is fortran style. while python is C++
        return {data.dtype.names[j]: data[i][j] * self.scaling_factor for j in range(len(data.dtype))}

    def analyse_mutation(self, alt_resn: str, _failed=0) -> Dict:
        self.do_relax()
        self.mark('relaxed')
        self.native = self.pose.clone()
        nblock = self.output_pdbblock()
        self.mutate(alt_resn)
        self.mark('mutate')
        self.do_relax()
        self.mark('mutarelax')
        if self.scores['mutate'] < self.scores['mutarelax']:  # increment!
            if _failed < 2:
                log.debug(f"Incrementor issue: {self.scores['mutate']} < {self.scores['mutarelax']}")
                self.pose = self.native.clone()
                self.preminimize(3)
                return self.analyse_mutation(alt_resn, _failed + 1)
            else:
                self.scores['mutarelax'] = self.scores['mutate']
        return {'ddG':                  self.scores['mutarelax'] - self.scores['relaxed'],
                'scores':               self.scores,
                'native':               nblock,
                'mutant':               self.output_pdbblock(),  # pdbb
                'rmsd':                 pyrosetta.rosetta.core.scoring.CA_rmsd(self.native, self.pose),
                'dsol':                 self.get_diff_solubility(),
                'score_fxn':            self.scorefxn.get_name(),
                'ddG_residue':          self.get_diff_res_score(),
                'native_residue_terms': self.get_res_score_terms(self.native),
                'mutant_residue_terms': self.get_res_score_terms(self.pose),
                'terms':                self.get_all_scores(),
                'neighbours':           self.get_pdb_neighbours(),
                'cycles':               self.cycles,
                'radius':               self.radius,
                'neighbouring_ligand':  self.is_ligand_in_sele(),
                'n_constraints':        len(self.pose.constraint_set().get_all_constraints())
                }

    def preminimize(self, expansion: int):
        log.info('Emergency preminimisation')
        # alter vector
        self.neighbour_vector = self.calculate_neighbours_in_pyrosetta(self.radius + expansion)
        self.ready_relax(5, fixed_bb=True)
        self.do_relax()
        # reset — recalculated.
        self.neighbour_vector = self.calculate_neighbours_in_pyrosetta(self.radius)
