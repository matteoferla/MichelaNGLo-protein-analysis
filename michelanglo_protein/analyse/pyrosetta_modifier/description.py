from .relax import MutatorRelax

import pyrosetta_help as ph
import pyrosetta
from typing import List

# mro: MutatorBase -> MutatorInit -> MutatorNeighbors -> MutatorCon -> MutatorRelax -> MutatorDescribe -> Mutator
class MutatorDescribe(MutatorRelax):
    def describe_neighbors(self) -> List[dict]:
        # pose_idx, pdb_idx pdb_chain resn is_protein omega ss betaturn hbonds
        descriptions: List[dict] = []
        cis = ph.get_cis_residues(self.pose)
        beta = ph.get_betaturns(self.pose)
        ph.get_ss(self.pose)  # fill secstruct
        hbonds = ph.get_hbond_dicts(self.pose)
        pdb_info = self.pose.pdb_info()
        # self.neighbour_vector is a pyrosetta.rosetta.utility.vector1_bool
        # whereas I want a pyrosetta.rosetta.core.select.residue_selector.ResidueVector,
        # which is a pyrosetta.rosetta.utility.vector1_unsigned_long,
        for neigh_i in pyrosetta.rosetta.core.select.residue_selector.ResidueVector( self.neighbour_vector ):
            residue = self.pose.residue(neigh_i)
            resi_hbonds = hbonds[neigh_i] if neigh_i in hbonds else []  # technically it should be there
            # get_hbond_dicts does not clean split by own and other:
            for hbond in resi_hbonds:
                if hbond['acc_resi'] == neigh_i:
                    hbond['direction'] = 'acceptor'
                    hbond['other_pdb_idx'] = pdb_info.number(hbond['don_resi'])
                    hbond['other_pdb_chain'] = pdb_info.chain(hbond['don_resi'])
                    hbond['other_atm_name'] = hbond['don_atm_name']
                    hbond['own_atm_name'] = hbond['acc_atm_name']
                else:
                    hbond['direction'] = 'donor'
                    hbond['other_pdb_idx'] = pdb_info.number(hbond['acc_resi'])
                    hbond['other_pdb_chain'] = pdb_info.chain(hbond['acc_resi'])
                    hbond['other_atm_name'] = hbond['acc_atm_name']
                    hbond['own_atm_name'] = hbond['don_atm_name']
            descriptions.append(
                dict(pose_idx=neigh_i,
                     pdb_idx=pdb_info.number(neigh_i),
                     pdb_chain=pdb_info.chain(neigh_i),
                     resn=residue.name3(),
                     is_protein=residue.is_protein(),
                     omega='cis' if neigh_i in cis else 'trans',
                     # unlike SS string, which skips non-peptides this will give L for them:
                     ss=self.pose.secstruct(neigh_i),
                     betaturn=beta[neigh_i] if neigh_i in beta else '',
                     hbonds=resi_hbonds,
                     )
            )
        return descriptions
