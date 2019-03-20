import unittest
from protein import Protein


class TestProtein(unittest.TestCase):

    def test_warn(self):
        print('Two userwarnings coming up')
        Protein.settings.verbose = True
        Protein.settings.missing_attribute_tolerant = True # default for now but may change
        irak = Protein()
        with self.assertWarns(UserWarning) as cm:
            foo = irak.foo
        with self.assertWarns(UserWarning) as cm:
            irak.other['bar'] = 'bar'
            self.assertEqual(irak.bar, 'bar')

    def test_error_tollerant(self):
        print('testing failsafe decorator')
        Protein.settings.error_tolerant = False
        mock = Protein()
        try:
            self.assertRaises(ValueError, mock._test_failsafe)
        except:
            pass
        Protein.settings.error_tolerant = True
        self.assertIsNone(mock._test_failsafe())

    def test_parse(self):
        print('testing Uniprot parsing')
        Protein.settings.error_tolerant = False
        Protein.settings.tolerate_missing_attributes = False
        irak = Protein()
        irak.uniprot = 'Q9NWZ3'
        irak.gene_name = 'IRAK4'
        irak.parse_uniprot()
        self.assertEqual(irak.sequence[0], 'M')


    def test_full_parse(self):
        print('testing serial parsing')
        Protein.settings.error_tollerant = False
        Protein.settings.tollerate_missing_attributes = False
        Protein.settings.verbose = True
        irak = Protein(uniprot = 'Q9NWZ3', gene_name = 'IRAK4')
        irak.parse_all(mode='serial')


if __name__ == '__main__':
    print('*****Test********')

    unittest.main()
    #irak.parse_uniprot()
