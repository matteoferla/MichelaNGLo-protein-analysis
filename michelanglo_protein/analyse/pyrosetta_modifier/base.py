from .params import get_default_params_filenames
import pyrosetta
import logging
from collections import namedtuple, defaultdict

log = logging.getLogger()
Target = namedtuple('target', ['resi', 'chain'])

# mro: MutatorBase -> MutatorInit -> MutatorNeighbors -> MutatorCon -> MutatorRelax -> MutatorDescribe -> Mutator
class MutatorBase:
    scaling_factor = 0.239  # multiplied scale. This is also set in the app (to the same value, just because)

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
        "p_aa_pp":              "Probability of amino acid at &phi;/&phi;.",
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

    default_params = get_default_params_filenames()
    aminoacid_letters = ('A', 'C', 'D', 'E', 'F', 'G',
                         'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q',
                         'R', 'S', 'T', 'V', 'W', 'Y')

    def get_scorefunction(self, scorefxn_name: str) -> pyrosetta.ScoreFunction:
        """
        Preps and gets the scorefunction

        :param scorefxn_name: name
        :return:
        """
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
        return pyrosetta.create_score_function(scorefxn_name)

    @staticmethod
    def reinit(verbose: bool = False):
        if verbose:
            pyrosetta.init(options='-ignore_unrecognized_res true')
        else:
            pyrosetta.init(silent=True, options='-mute core basic protocols -ignore_unrecognized_res true')
