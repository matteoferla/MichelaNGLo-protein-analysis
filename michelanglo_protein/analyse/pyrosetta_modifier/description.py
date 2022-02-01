from .relax import MutatorRelax

import pyrosetta_help as ph
import pyrosetta
from typing import List

class ExtendedBondDataType(ph.BondDataType):
    direction: str  # acceptor/donor
    other_pdb_idx: int
    other_pdb_chain: str
    other_atm_name: str
    own_atm_name: str

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
        for neigh_i in pyrosetta.rosetta.core.select.residue_selector.ResidueVector( self.neighbour_vector ):  # noqa
            residue = self.pose.residue(neigh_i)
            xlink = {}
            if ph.is_xlinked(residue):
                xlink = ph.get_xlink_details(neigh_i, self.pose)
                # pose2pdb:
                xlink['other_pdb_chain'] = pdb_info.chain(xlink['other_idx'])
                xlink['other_pdb_idx'] = pdb_info.number(xlink['other_idx'])
            # get_hbond_dicts does not clean split by own and other:
            raw_resi_hbonds = hbonds[neigh_i] if neigh_i in hbonds else []  # technically it should be there
            neigh_hbonds: List[ExtendedBondDataType] = []
            for hbond in raw_resi_hbonds:
                if hbond['acc_resi'] == neigh_i:
                    # acceptor
                    neigh_hbonds.append(ExtendedBondDataType(**hbond,
                                                             direction='acceptor',
                                                             other_pdb_idx=pdb_info.number(hbond['don_resi']),
                                                             other_pdb_chain=pdb_info.chain(hbond['don_resi']),
                                                             other_atm_name=hbond['don_atm_name'],
                                                             own_atm_name=hbond['acc_atm_name']
                                                             )
                                        )
                else:
                    # donor
                    neigh_hbonds.append(ExtendedBondDataType(**hbond,
                                                             direction='donor',
                                                             other_pdb_idx=pdb_info.number(hbond['acc_resi']),
                                                             other_pdb_chain=pdb_info.chain(hbond['acc_resi']),
                                                             other_atm_name=hbond['acc_atm_name'],
                                                             own_atm_name=hbond['don_atm_name']
                                                             )
                                        )
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
                     hbonds=neigh_hbonds,
                     xlink=xlink
                     )
            )
        return descriptions
