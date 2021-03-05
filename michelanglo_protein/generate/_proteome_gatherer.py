__doc__ = """
Generates the data to be used.

However, it repeated itself. Hence the switch to uniprot_master_parser.UniprotMasterReader.
"""

from pprint import PrettyPrinter
pprint = PrettyPrinter().pprint

import os, json
from warnings import warn
from ._protein_gatherer import ProteinGatherer
from .uniprot_master_parser import UniprotMasterReader
from .split_gnomAD import gnomAD
from .PDB_blast import Blaster
from ..settings_handler import global_settings #the instance not the class.
import random


def announce(*msgs):
    print('*' * 50)
    print(*msgs)
    print('*' * 50)

class ProteomeGathererOLD:
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
        warn('This script is depractated in favour of uniprot_master_parser.UniprotMasterReader', category=DeprecationWarning)
        ######## FETCH ALL RAW FILES #########################
        ProteinGatherer.settings.verbose = True
        announce('Retrieving references')
        if not skip:
            ProteinGatherer.settings.retrieve_references(ask=False)

        ######## PARSE UNIPROT ############################
        self.master_file = os.path.join(ProteinGatherer.settings.temp_folder, 'uniprot_sprot.xml')
        # This class parses the uniprot FTP file and can do various things. such as making a small one that is only human.
        # But mainly the `UniprotMasterReader.convert('uniprot_sprot.xml')` method whcih generates the JSON files required.
        # first_n_protein is for testing.
        announce('Convert Master Uniprot file')
        if not skip:
            #UniprotMasterReader.convert(uniprot_master_file = self.master_file, first_n_protein=0)
            UniprotMasterReader()

        ######## BLAST PDB ############################
        # uncompresses the pdbaa
        announce('Extracting blast db')
        if not skip:
            Blaster.extract_db()
        # blasts the human.fa against the newly extracted pdbaa
        announce('Blasting')
        if not skip:
            Blaster.pdb_blaster()
        # converts the files to something reasonable.
        announce('Parsing blast output')
        if not skip:
            Blaster.parse('blastpdb','blastpdb2')
        ######## ASSEMBLE ############################
        announce('Assembing proteome')
        ProteinGatherer.settings.fetch = False
        ProteinGatherer.settings.error_tolerant = False
        ProteinGatherer.settings.missing_attribute_tolerant = False
        ProteinGatherer.settings.verbose = True
        uniprot_ids = list(set(json.load(open(os.path.join(ProteinGatherer.settings.data_folder,'human_prot_namedex.json'))).values()))
        random.shuffle(uniprot_ids)   # debug...
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

class ProteomeGatherer:
    settings = global_settings

    def __init__(self, data_folder: str):
        ######## FETCH ALL RAW FILES #########################
        self.settings.verbose = True
        announce('Retrieving references')
        self.settings.verbose = True  # False
        self.settings.startup(data_folder=data_folder)
        self.settings.retrieve_references(ask=False, refresh=False)
        self.settings.error_tolerant = True
        announce('Parsing Uniprot')
        UniprotMasterReader()
        announce('Splitting gnomAD files')
        gnomAD(genomasterfile=os.path.join(self.settings.reference_folder,
                                           'gnomAD.genomes.r2.1.1.exome_calling_intervals.sites.vcf.bgz'),
               exomasterfile=os.path.join(self.settings.reference_folder, 'gnomAD.exomes.r2.1.1.sites.vcf.bgz'),
               namedexfile=os.path.join(self.settings.dictionary_folder, 'taxid9606-names2uniprot.json'),
               folder=os.path.join(self.settings.temp_folder, 'gnomAD')
               ).split()
        announce('Adding gnomAD files')
        path = os.path.join(global_settings.pickle_folder, f'taxid9606')
        for pf in os.listdir(path):
            try:
                protein = ProteinGatherer().load(file=os.path.join(path, pf))
                protein.gnomAD = []
                protein.parse_gnomAD()
                protein.get_PTM()
                protein.dump()
            except:
                pass
