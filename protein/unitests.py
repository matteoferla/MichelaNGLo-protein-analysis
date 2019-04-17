import unittest
from . import ProteinCore


class TestProteinCore(unittest.TestCase):

    def test_warn(self):
        print('Two userwarnings coming up')
        ProteinCore.settings.verbose = True
        ProteinCore.settings.missing_attribute_tolerant = True # default for now but may change
        irak = ProteinCore()
        with self.assertWarns(UserWarning) as cm:
            foo = irak.foo
        with self.assertWarns(UserWarning) as cm:
            irak.other['bar'] = 'bar'
            self.assertEqual(irak.bar, 'bar')

    def test_error_tollerant(self):
        print('testing failsafe decorator')
        ProteinCore.settings.error_tolerant = False
        mock = ProteinCore()
        try:
            self.assertRaises(ValueError, mock._test_failsafe)
        except:
            pass
        ProteinCore.settings.error_tolerant = True
        self.assertIsNone(mock._test_failsafe())

    def test_parse(self):
        print('testing Uniprot parsing')
        ProteinCore.settings.error_tolerant = False
        ProteinCore.settings.tolerate_missing_attributes = False
        irak = ProteinCore()
        irak.uniprot = 'Q9NWZ3'
        irak.gene_name = 'IRAK4'
        irak.parse_uniprot()
        self.assertEqual(irak.sequence[0], 'M')


    def test_full_parse(self):
        print('testing serial parsing')
        ProteinCore.settings.error_tollerant = False
        ProteinCore.settings.tollerate_missing_attributes = False
        ProteinCore.settings.verbose = True
        irak = ProteinCore(uniprot = 'Q9NWZ3', gene_name = 'IRAK4')
        irak.parse_all(mode='serial')


if __name__ == '__main__':
    print('*****Test********')

    unittest.main()
    #irak.parse_uniprot()
