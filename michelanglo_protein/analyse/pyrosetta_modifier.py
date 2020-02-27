__doc__ = """
This file does all the pyrosetta operations. The energetical or ddG variable in the API of VENUS.

It is called by ``ProteinAnalyser.analyse_FF``.
To avoid segmentation faults it is run on a separate process byt this.

Pyrosetta will throw a segmentation fault if anything is done incorrectly. Such as editing a non-existent atom.
As a result ProteinAnalyser.analyse_FF uses multiprocessing to do the job on a different core.
"""


import pyrosetta, pymol2, re, os
from typing import List, Dict #, TypedDict
from collections import namedtuple



pyrosetta.init(silent=True, options='-mute core basic protocols -ignore_unrecognized_res true')

Target = namedtuple('target', ['resi', 'chain'])


class Mutator:
    """
    Relaxes around a residue on init and mutates.

    * ``.target`` mutation see Target namedtuple (resi, chain)
    * ``.neighbours`` list of neighbours
    * ``.pdbblock`` str pdb block
    * ``.pose`` pyrosetta.Pose
    * ``._pdb2pose`` points to ``self.pose.pdb_info().pdb2pose``, while target_pdb2pose accepts Target and gives back int
    """

    def __init__(self, pdbblock: str, target_resi: int, target_chain: str = 'A', cycles: int = 1, radius: int = 4):
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
        """
        self.scorefxn = pyrosetta.get_fa_scorefxn()
        self.scores = {} #gets filled by .mark()
        ## Load
        self.target = Target(target_resi, target_chain)
        self.pdbblock = pdbblock
        self.load_pose()  # fills self.pose
        self.mark('raw')

        ## Find neighbourhood
        self.calculate_neighbours(radius)  # fills self.neighbours

        ## Read relax
        self.ready_relax(cycles)

    def target_pdb2pose(self, target:Target) -> int:
        return self._pdb2pose(chain=target.chain, res=target.resi)

    @staticmethod
    def reinit(verbose:bool =False):
        if verbose:
            pyrosetta.init(options='-ignore_unrecognized_res true')
        else:
            pyrosetta.init(silent=True, options='-mute core basic protocols -ignore_unrecognized_res true')

    def load_pose(self) -> pyrosetta.Pose:
        """
        Loading from str is a bit messy. this simply does that and returns a Pose

        :return: self.pose
        """
        self.pose = pyrosetta.Pose()
        pyrosetta.rosetta.core.import_pose.pose_from_pdbstring(self.pose, self.pdbblock)
        self._pdb2pose = self.pose.pdb_info().pdb2pose
        return self.pose

    def calculate_neighbours(self, radius: int = 4) -> List[Target]:
        """
        Gets the residues within the radius of target. THis method uses PyMOL!
        It fills self.neighbours and returns it.

        :return: self.neighbours
        :rtype: List[Target]
        """
        with pymol2.PyMOL() as pymol:
            pymol.cmd.read_pdbstr(self.pdbblock, 'blockprotein')
            sele = f"name CA and (byres chain {self.target.chain} and resi {self.target.resi} around {radius})"
            atoms = pymol.cmd.get_model(sele)
            self.neighbours = []
            for atom in atoms.atom:
                res = int(re.match('\d+', atom.resi).group())
                self.neighbours.append(Target(resi=res, chain=atom.chain))
        return self.neighbours

    def ready_relax(self, cycles: int = 1) -> pyrosetta.rosetta.protocols.moves.Mover:
        """

        :param cycles:
        :return:
        """
        self.relax = pyrosetta.rosetta.protocols.relax.FastRelax(self.scorefxn, cycles)
        self.movemap = pyrosetta.MoveMap()
        for neigh in self.neighbours:
            n = self.target_pdb2pose(neigh)
            self.movemap.set_bb(n, True)
            self.movemap.set_chi(n, True)
        self.relax.set_movemap(self.movemap)
        return self.relax

    def mark(self, label: str) -> Dict:
        """
        Save the score to ``.scores``

        :param label: scores is Dict. label is key.
        :return: {ddG: float, scores: Dict[str, float], native:str, mutant:str, rmsd:int}
        """
        self.scores[label] = self.scorefxn(self.pose)
        return self.scores

    def mutate(self, aa):
        pyrosetta.toolbox.mutate_residue(self.pose, self.target_pdb2pose(self.target), aa)

    def do_relax(self):
        self.relax.apply(self.pose)

    def output_pdbblock(self, tmpfile='temp.pdb'):
        """
        This is weird. I did not find the equivalent to ``pose_from_pdbstring``.
        Even though the std::osstream is one of the options for dump_pdb (``p.dump_pdb(stream)``), I cannot make one.
        Namely, ``stream = pyrosetta.rosetta.std.ostream()`` requires a ``pyrosetta.rosetta.std.streambuf`` object.
        So I am doing it the stupid way.

        :param tmpfile: filename
        :return:
        """
        self.pose.dump_pdb(tmpfile)
        block = open(tmpfile,'r').read()
        if os.path.exists(tmpfile):
            os.remove(tmpfile)
        return block

    def get_diff_solubility(self) -> float:
        """
        Gets the difference in solubility (fa_sol) for the protein.
        fa_sol = Gaussian exclusion implicit solvation energy between protein atoms in different residue
        :return: fa_sol kcal/mol
        """
        n = self.scorefxn.score_by_scoretype(self.native, pyrosetta.rosetta.core.scoring.ScoreType.fa_sol)
        m = self.scorefxn.score_by_scoretype(self.pose, pyrosetta.rosetta.core.scoring.ScoreType.fa_sol)
        return m - n

    def get_diff_res_score(self) -> float:
        """
        Gets the difference in score for that residue

        :return: per_residue kcal/mol
        """
        #segfaults if score is not run globally first!
        i = self.target_pdb2pose(self.target)
        r = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(i)
        n = self.scorefxn.get_sub_score(self.native, r.apply(self.native))
        m = self.scorefxn.get_sub_score(self.pose, r.apply(self.native))
        return m - n

    def get_res_score_terms(self, pose) -> dict:
        data = pose.energies().residue_total_energies_array() #structured numpy array
        i = self.target_pdb2pose(self.target) + 1
        return {data.dtype.names[j]: data[i][j] for j in range(len(data.dtype))}

    def analyse_mutation(self, alt_resn:str) -> Dict:
        self.do_relax()
        self.mark('relaxed')
        self.native = self.pose.clone()
        nblock = self.output_pdbblock()
        self.mutate(alt_resn)
        self.mark('mutate')
        self.do_relax()
        self.mark('mutarelax')
        return {'ddG': self.scores['mutarelax'] - self.scores['relaxed'],
                'scores': self.scores,
                'native': nblock,
                'mutant': self.output_pdbblock(), ## pdbb
                'rmsd': pyrosetta.rosetta.core.scoring.CA_rmsd(self.native, self.pose),
                'dsol': self.get_diff_solubility(),
                'score_fxn': self.scorefxn.get_name(),
                'ddG_residue': self.get_diff_res_score(),
                'native_residue_terms': self.get_res_score_terms(self.native),
                'mutant_residue_terms': self.get_res_score_terms(self.pose)
                }


#######################################################################################################################

def test():
    ## these are tests.
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
            #p.terminate()
            break
        elif not p.is_alive():
            child_conn.send({'error': 'segmentation fault'})
            break
    msg = parent_conn.recv()
    print('DONE!')

if __name__ == '__main__':
    #test()
    paratest()



