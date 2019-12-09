__doc__ = """
Generates the data to be used.

However, it repeated itself. Hence the switch to uniprot_master_parser.UniprotReader.
"""

from pprint import PrettyPrinter
pprint = PrettyPrinter().pprint

import os, json
from warnings import warn
from ._protein_gatherer import ProteinGatherer
from .uniprot_master_parser import UniprotReader
from .PDB_blast import Blaster
import random


def announce(*msgs):
    print('*' * 50)
    print(*msgs)
    print('*' * 50)

class ProteomeGatherer:
    """
    Gets everything ready and parses!
    >>> ProteomeGatherer(skip=True)
    Will assume all is prepared.
    """

    def __init__(self, skip=False, remake_pickles=False):
        """
        Calls the variaous parts that get things ready.
        :param skip: boolean, if true it runs parse_proteome
        """
        ################ FETCH ALL RAW FILES #################################################
        ProteinGatherer.prosettings.verbose = True
        announce('Retrieving references')
        if not skip:
            ProteinGatherer.settings.retrieve_references(ask=False)

        ################ PARSE UNIPROT #######################################################
        self.master_file = os.path.join(ProteinGatherer.settings.temp_folder, 'uniprot_sprot.xml')
        # This class parses the uniprot FTP file and can do various things. such as making a small one that is only human.
        # But mainly the `UniprotReader.convert('uniprot_sprot.xml')` method whcih generates the JSON files required.
        # first_n_protein is for testing.
        announce('Convert Master Uniprot file')
        if not skip:
            UniprotReader.convert(uniprot_master_file = self.master_file, first_n_protein=0)

        ################ BLAST PDB #######################################################
        # uncompresses the pdbaa
        announce('Extracting blast db')
        if not skip:
            Blaster.extract_db()
        # blasts the human.fa agains the newly extracted pdbaa
        announce('Blasting')
        if not skip:
            Blaster.pdb_blaster()
        # converts the files to something reasonable.
        announce('Parsing blast output')
        if not skip:
            Blaster.parse('blastpdb','blastpdb2')
        ################ ASSEMBLE #######################################################
        announce('Assembing proteome')
        ProteinGatherer.settings.fetch = False
        ProteinGatherer.settings.error_tolerant = False
        ProteinGatherer.settings.missing_attribute_tolerant = False
        ProteinGatherer.settings.verbose = True
        uniprot_ids = list(set(json.load(open(os.path.join(ProteinGatherer.settings.data_folder,'human_prot_namedex.json'))).values()))
        random.shuffle(uniprot_ids)   ## debug...
        for acc in uniprot_ids:
            prot = ProteinGatherer(uniprot=acc)
            if not remake_pickles:
                prot.gload()
            else:
                print(prot.uniprot)
                prot.parse_uniprot()
                prot.parse_all(mode='parallel')
                prot.complete()
                warn(f'lines silenced for now at {__name__} line 77')
                prot.get_percent_modelled()
                #prot.parse_ExAC_type()
                prot.gdump()

