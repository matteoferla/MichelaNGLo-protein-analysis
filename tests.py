import unittest


# from michelanglo_app.transplier import PyMolTranspiler
#
# transpiler = PyMolTranspiler(file='michelanglo_app/demo/1gfl.pse')
# #print(transpiler.get_view())
# #print(transpiler.get_reps())
# print(transpiler.get_loadfun_js(tag_wrapped=True, viewport='viewport', funname='loadfun'))


from michelanglo_protein import ProteinAnalyser, Structure, Mutation, global_settings
from michelanglo_protein.analyse import Mutator

import unittest, requests, json, time, pyrosetta

class MutatorTest(unittest.TestCase):
    def setUp(self):
        # 1SFT is alanine racemase
        self.pdbblock = requests.get('https://files.rcsb.org/download/1SFT.pdb').text

    def tearDown(self):
        pass

    def test_basic(self):
        # formerly from michelanglo_protein.analyse.pyrosetta_modifier import test
        # 1SFT/A/A/HIS`166
        tick = time.time()
        m = Mutator(pdbblock=self.pdbblock, target_resi=166, target_chain='A', cycles=1, radius=3)
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

    def test_main_mutator(self):
        import requests
        pa = ProteinAnalyser()
        pa.structural = Structure('test', 'test w alaR',
                                  1, 999, 'test',
                                  type='custom',
                                  coordinates=self.pdbblock
                                  )
        pa.mutation = Mutation('A166W')
        pa.analyse_FF(spit_process=False)

    def test_sub_mutator(self):
        import requests
        pa = ProteinAnalyser()
        pa.structural = Structure('test', 'test w alaR',
                                  1, 999, 'test',
                                  type='custom',
                                  coordinates=self.pdbblock
                                  )
        pa.mutation = Mutation('A166W')
        pa.analyse_FF(spit_process=True)


class FullTest(unittest.TestCase):
    def enable_logging(self):
        import logging
        from michelanglo_protein.analyse.pyrosetta_modifier import log
        handler = logging.StreamHandler()
        handler.setLevel(logging.DEBUG)
        log.setLevel(logging.DEBUG)
        log.addHandler(handler)
        return log

    def get_protein(self, uniprot_id:str, mutation:str) -> ProteinAnalyser:
        protein = ProteinAnalyser(uniprot=uniprot_id, taxid=9606)
        protein.load()
        protein.add_alphafold2()
        protein.mutation = Mutation(mutation)
        protein.predict_effect()
        return protein

    def setUp(self):
        self.log = self.enable_logging()
        global_settings.startup('../protein-data')

    def tearDown(self):
        pass

    def test_legacy(self):
        protein:ProteinAnalyser = self.get_protein('Q8WXF1', 'A175D')
        protein.predict_effect(full=True)

    def test_full(self):
        protein:ProteinAnalyser = self.get_protein('Q8WXF1', 'A175D')
        # protein.retrieve_structures_from_swissmodel()
        structure = protein.get_best_model()
        protein.analyse_structure(structure)
        out = protein.analyse_FF(spit_process=True)  # calls analyse_mutation
        self.log.info(out)
        self.assertIsInstance(out['ddG'], float)
        self.assertLess(out['ddG'], 5)

    def test_gnomad(self):
        protein:ProteinAnalyser = self.get_protein('Q8WXF1', 'A175D')
        protein.analyse_structure()
        out = protein.analyse_gnomad_FF(spit_process=True)
        self.log.info(out)


    def test_alpha(self):
        protein:ProteinAnalyser = self.get_protein('O75911', 'L273D')
        protein.add_alphafold2()
        protein.pdbs =[]
        protein.swissmodel = []
        protein.predict_effect()
        # protein.retrieve_structures_from_swissmodel()
        # protein.add_alphafold2()
        protein.analyse_structure()
        out = protein.analyse_FF(spit_process=True)  # calls analyse_mutation
        self.log.info(out)
        self.assertIsInstance(out['ddG'], float)
        self.assertLess(out['ddG'], 5)

    def untest_adhoc(self):
        # test to try with sefault settings
        # add to `.analyse_FF` this:
        # with open('/Users/matteo/Desktop/test.json', 'w') as f:
        #     json.dump(init_settings, f)
        with open('/Users/matteo/Desktop/test.json', 'r') as f:
            settings = json.load(f)
        Mutator(**settings).analyse_mutation('D')

if __name__ == '__main__':
    unittest.main()
