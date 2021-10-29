__doc__ = """
This file does all the pyrosetta operations. The energetical or ddG variable in the API of VENUS.

It is called by ``ProteinAnalyser.analyse_FF``.
To avoid segmentation faults it is run on a separate process byt this.

Pyrosetta will throw a segmentation fault if anything is done incorrectly. Such as editing a non-existent atom.
As a result ProteinAnalyser.analyse_FF uses multiprocessing to do the job on a different core.

The PyMOL neighbours could be using

cc_sele = pyrosetta.rosetta.core.select.residue_selector.CloseContactResidueSelector()
cc_sele.central_residue_group_selector(resi_sele)
cc_sele.threshold(3)

But for debug this was done...
"""

# NB. Do not change the spelling of Neighbours to Neighbors as the front end uses it too.
# I did not realise that the Americans spell it without a ``u`` until it was too embedded.
# ``colour`` is correctly spelt ``color`` throughout.

import pyrosetta, pymol2, re, os
from typing import *
from collections import namedtuple, defaultdict
from Bio.SeqUtils import seq3
from ..gnomad_variant import Variant
import logging

log = logging.getLogger()

Target = namedtuple('target', ['resi', 'chain'])

default_params_folder = os.path.join(os.path.split(__file__)[0], 'params')

pyrosetta.init(silent=True, options='-mute core basic protocols -ignore_unrecognized_res true')


class Mutator:
    """
    Relaxes around a residue on init and mutates.

    * ``.target`` mutation see Target namedtuple (resi, chain)
    * ``.neighbours`` list of neighbours
    * ``.pdbblock`` str pdb block
    * ``.pose`` pyrosetta.Pose
    * ``._pdb2pose`` points to ``self.pose.pdb_info().pdb2pose``, while target_pdb2pose accepts Target and gives back int
    """
    scaling_factor = 1  # multiplied scale

    term_meanings = defaultdict(str, {
        "fa_atr":               "Lennard-Jones attractive between atoms in different residues (r^6 term, London dispersion forces).",
        "fa_rep":               "Lennard-Jones repulsive between atoms in different residues (r^12 term, Pauli repulsion forces).",
        "fa_sol":               "Lazaridis-Karplus solvation energy.",
        "fa_intra_rep":         "Lennard-Jones repulsive between atoms in the same residue.",
        "fa_elec":              "Coulombic electrostatic potential with a distance-dependent dielectric.",
        "pro_close":            "Proline ring closure energy and energy of psi angle of preceding residue.",
        "hbond_sr_bb":          "Backbone-backbone hbonds close in primary sequence.",
        "hbond_lr_bb":          "Backbone-backbone hbonds distant in primary sequence.",
        "hbond_bb_sc":          "Sidechain-backbone hydrogen bond energy.",
        "hbond_sc":             "Sidechain-sidechain hydrogen bond energy.",
        "dslf_fa13":            "Disulfide geometry potential.",
        "rama":                 "Ramachandran preferences.",
        "omega":                "Omega dihedral in the backbone. A Harmonic constraint on planarity with standard deviation of ~6 deg.",
        "fa_dun":               "Internal energy of sidechain rotamers as derived from Dunbrack's statistics (2010 Rotamer Library used in Talaris2013).",
        "fa_dun_semi":          "Internal energy of sidechain semi-rotamers as derived from Dunbrack's statistics (2010 Rotamer Library used in Talaris2013).",
        "p_aa_pp":              "Probability of amino acid at Φ/Ψ.",
        "ref":                  "Reference energy for each amino acid. Balances internal energy of amino acid terms.  Plays role in design.",
        "METHOD_WEIGHTS":       "Not an energy term itself, but the parameters for each amino acid used by the ref energy term.",
        "lk_ball":              "Anisotropic contribution to the solvation.",
        "lk_ball_iso":          "Same as fa_sol; see below.",
        "lk_ball_wtd":          "weighted sum of lk_ball & lk_ball_iso (w1*lk_ball + w2*lk_ball_iso); w2 is negative so that anisotropic contribution(lk_ball) replaces some portion of isotropic contribution (fa_sol=lk_ball_iso).",
        "lk_ball_bridge":       "Bonus to solvation coming from bridging waters, measured by overlap of the 'balls' from two interacting polar atoms.",
        "lk_ball_bridge_uncpl": "Same as lk_ball_bridge, but the value is uncoupled with dGfree (i.e. constant bonus, whereas lk_ball_bridge is proportional to dGfree values).",
        "fa_intra_atr_xover4":  "Intra-residue LJ attraction, counted for the atom-pairs beyond torsion-relationship.",
        "fa_intra_rep_xover4":  "Intra-residue LJ repulsion, counted for the atom-pairs beyond torsion-relationship.",
        "fa_intra_sol_xover4":  "Intra-residue LK solvation, counted for the atom-pairs beyond torsion-relationship.",
        "fa_intra_elec":        "Intra-residue Coulombic interaction, counted for the atom-pairs beyond torsion-relationship.",
        "rama_prepro":          "Backbone torsion preference term that takes into account of whether preceding amono acid is Proline or not.",
        "hxl_tors":             "Sidechain hydroxyl group torsion preference for Ser/Thr/Tyr, supersedes yhh_planarity (that covers L- and D-Tyr only).",
        "yhh_planarity":        "Sidechain hydroxyl group torsion preference for Tyr, superseded by hxl_tors"
    })

    default_params = [os.path.join(default_params_folder, filename) for filename in os.listdir(default_params_folder)
                      if os.path.splitext(filename)[1] == '.params']

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
        if 'ref' in scorefxn_name:
            pass
        elif 'beta_july15' in scorefxn_name or 'beta_nov15' in scorefxn_name:
            pyrosetta.rosetta.basic.options.set_boolean_option('corrections:beta_july15', True)
        elif 'beta_nov16' in scorefxn_name:
            pyrosetta.rosetta.basic.options.set_boolean_option('corrections:beta_nov16', True)
        elif 'genpot' in scorefxn_name:
            pyrosetta.rosetta.basic.options.set_boolean_option('corrections:gen_potential', True)
        elif 'talaris' in scorefxn_name:
            pyrosetta.rosetta.basic.options.set_boolean_option(f'corrections:restore_talaris_behavior', True)
        else:
            log.warning(f'No correction applied for {scorefxn_name}!')
        # there are a few other fixes. Such as franklin2019 and spades.
        self.scorefxn = pyrosetta.create_score_function(scorefxn_name)
        ap_st = pyrosetta.rosetta.core.scoring.ScoreType.atom_pair_constraint
        self.scorefxn.set_weight(ap_st, 10)
        # correct for split as per https://www.rosettacommons.org/node/11245
        weights = self.scorefxn.weights()  # Create the EnergyMap
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
        self.ready_relax(self.cycles)
        self.n_preventions = 0
        self.prevent_acceptance_of_incrementor = prevent_acceptance_of_incrementor
        if outer_constrained:
            self.explicit_cons()

    def target_pdb2pose(self, target: Target) -> int:
        return self._pdb2pose(chain=target.chain, res=target.resi)

    @staticmethod
    def reinit(verbose: bool = False):
        if verbose:
            pyrosetta.init(options='-ignore_unrecognized_res true')
        else:
            pyrosetta.init(silent=True, options='-mute core basic protocols -ignore_unrecognized_res true')

    def load_pose(self, single_chain=False, remove_ligands=False) -> pyrosetta.Pose:
        """
        Loading from str is a bit messy. this simply does that and returns a Pose

        :return: self.pose
        """
        pose = pyrosetta.Pose()
        if self.params_filenames:
            params_paths = pyrosetta.rosetta.utility.vector1_string()
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

    @property
    def _pdb2pose(self):
        return self.pose.pdb_info().pdb2pose

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
        with pymol2.PyMOL() as pymol:
            pymol.cmd.read_pdbstr(self.pdbblock, 'blockprotein')
            sele = f"name CA and (byres chain {self.target.chain} and resi {self.target.resi} around {radius})"
            atoms = pymol.cmd.get_model(sele)
            neighbours = []
            for atom in atoms.atom:
                res = int(re.match(r'\d+', atom.resi).group())
                neighbours.append(Target(resi=res, chain=atom.chain))
        return neighbours

    def targets2vector(self, targets: List[Target]) -> pyrosetta.rosetta.utility.vector1_bool:
        neighbours = pyrosetta.rosetta.utility.vector1_bool(self.pose.total_residue())
        for target in targets:
            r = self.target_pdb2pose(target)
            neighbours[r] = True
        return neighbours

    def calculate_neighbours_in_pyrosetta(self, radius: int = 12) -> pyrosetta.rosetta.utility.vector1_bool:
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

    def ready_relax(self, cycles: int = 1, fixed_bb: bool = False) -> pyrosetta.rosetta.protocols.moves.Mover:
        """

        :param cycles:
        :return:
        """
        self.relax = pyrosetta.rosetta.protocols.relax.FastRelax(self.scorefxn, cycles)
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

    def output_pdbblock(self, pose: Optional[pyrosetta.Pose] = None) -> str:
        """
        This is weird. I did not find the equivalent to ``pose_from_pdbstring``.
        But using buffer works.

        :return: PDBBlock
        """
        if pose is None:
            pose = self.pose
        buffer = pyrosetta.rosetta.std.stringbuf()
        pose.dump_pdb(pyrosetta.rosetta.std.ostream(buffer))
        return buffer.str()

    def get_diff_solubility(self) -> float:
        """
        Gets the difference in solubility (fa_sol) for the protein.
        fa_sol = Gaussian exclusion implicit solvation energy between protein atoms in different residue

        :return: fa_sol kcal/mol
        """
        _, _, diff = self.get_term_scores(pyrosetta.rosetta.core.scoring.ScoreType.fa_sol)
        return diff

    def get_all_scores(self) -> Dict[str, Dict[str, Union[float, str]]]:
        data = {}
        for term in self.scorefxn.get_nonzero_weighted_scoretypes():
            data[term.name] = dict(zip(['native', 'mutant', 'difference'], self.get_term_scores(term)))
            data[term.name]['weight'] = self.scorefxn.get_weight(term)
            data[term.name]['meaning'] = self.term_meanings[term.name]
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

    def make_phospho(self, ptms):
        phospho = self.pose.clone()
        MutateResidue = pyrosetta.rosetta.protocols.simple_moves.MutateResidue
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

    def _repack_gnomad(self, pose_idx, from_resi, to_resi) -> int:
        self.pose = self.native.clone()
        self.pose.remove_constraints()
        # local repack...
        pyrosetta.toolbox.mutate_residue(self.pose,
                                         mutant_position=pose_idx,
                                         mutant_aa=from_resi,
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
        for record in gnomads:
            if record.type == 'nonsense':
                continue
            n = pose2pdb(chain='A', res=record.x)
            if n == 0:
                continue
            rex = re.match(r'(\w)(\d+)([\w])', record.description)
            if rex is None:
                continue
            if rex.group(0) in ddG:
                # print('duplicate mutation.')
                continue
            ddG[rex.group(0)] = self._repack_gnomad(n, rex.group(1), rex.group(3))
        return ddG


############################################################

def test():
    # these are tests.
    import requests, time
    # 1SFT/A/A/HIS`166
    pdbblock = requests.get('https://files.rcsb.org/download/1SFT.pdb').text
    tick = time.time()
    m = Mutator(pdbblock=pdbblock, target_resi=166, target_chain='A', cycles=1, radius=3)
    tock = time.time()
    print('LOAD', tock - tick)
    m.do_relax()
    m.mark('relaxed')
    tack = time.time()
    print('RELAX', tack - tock)
    native = m.pose.clone()
    m.mutate('P')
    m.mark('mutate')
    m.do_relax()
    m.mark('mutarelax')
    teck = time.time()
    print('MutaRelax', teck - tack)
    print(m.scores)
    muta = m.target_pdb2pose(m.target)
    print(pyrosetta.rosetta.core.scoring.CA_rmsd(native, m.pose, muta - 1, muta + 1))
    m.output()


def paratest():
    import requests, time
    from multiprocessing import Pipe, Process
    # 1SFT/A/A/HIS`166
    pdbblock = requests.get('https://files.rcsb.org/download/1SFT.pdb').text
    kwargs = dict(pdbblock=pdbblock, target_resi=166, target_chain='A', cycles=1, radius=3)

    def subpro(child_conn, **kwargs):  # Pipe <- Union[dict, None]:
        try:
            print('started child')
            Mutator.reinit()
            mut = Mutator(**kwargs)
            data = mut.analyse_mutation('W')  # {ddG: float, scores: Dict[str, float], native:str, mutant:str, rmsd:int}
            print('done', len(data))
            child_conn.send(data)
            print('completed child')
        except BaseException as error:
            print('error child')
            child_conn.send({'error': f'{error.__class__.__name__}:{error}'})

    parent_conn, child_conn = Pipe()
    p = Process(target=subpro, args=((child_conn),), kwargs=kwargs, name='pyrosetta')
    p.start()

    while 1:
        time.sleep(5)
        print(parent_conn.poll())
        if parent_conn.poll():
            # p.terminate()
            break
        elif not p.is_alive():
            child_conn.send({'error': 'segmentation fault'})
            break
    msg = parent_conn.recv()
    print('DONE!')


if __name__ == '__main__':
    test()
    # paratest()
