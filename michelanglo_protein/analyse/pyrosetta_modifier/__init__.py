from .params import get_default_params_filenames
from .base import MutatorBase, log, Target
from .description import MutatorDescribe

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

import pyrosetta

# todo should this line be here or in the app initialisation?
# it is also very simplistic. But pyrosetta_help params generation on the fly is too slow though.
pyrosetta.init(silent=True, options='-ex1 -ex2 -mute core basic protocols -ignore_unrecognized_res true')

# mro: MutatorBase -> MutatorInit -> MutatorNeighbors -> MutatorCon -> MutatorRelax -> MutatorDescribe -> Mutator
class Mutator(MutatorDescribe):
    """
    Relaxes around a residue on init and mutates.

    * ``.target`` mutation see Target namedtuple (resi, chain)
    * ``.neighbours`` list of neighbours
    * ``.pdbblock`` str pdb block
    * ``.pose`` pyrosetta.Pose
    * ``._pdb2pose`` points to ``self.pose.pdb_info().pdb2pose``, while target_pdb2pose accepts Target and gives back int
    """


############################################################

if __name__ == '__main__':
    from .test import test, paratest

    test()
    # paratest()
