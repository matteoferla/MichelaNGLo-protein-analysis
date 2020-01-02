"""Download and parse all the protein data (use this, not proteome gatherer)"""

from michelanglo_protein import ProteinAnalyser, ProteinCore, Mutation, Structure
from michelanglo_protein.settings_handler import global_settings
from michelanglo_protein.generate import ProteinGatherer, ProteomeGatherer
from michelanglo_protein.generate.split_gnomAD import gnomAD
from michelanglo_protein.generate.uniprot_master_parser import UniprotMasterReader
from michelanglo_protein.generate.split_phosphosite import Phosphosite
import os, json

if __name__ == '__main__':
    global_settings.verbose = True #False
    global_settings.startup(data_folder='../protein-data')
    #global_settings.retrieve_references(ask=False, refresh=False)
    global_settings.error_tolerant = True
    ## Phosphosite
    #Phosphosite().split().write()
    ## Uniprot
    UniprotMasterReader(first_n_protein=10)
    # gnomAD data needs to be split up after that the dictionaries are made.
    taxid=9606 #that's humans
    gnomAD(genomasterfile=os.path.join(global_settings.reference_folder,'gnomad.genomes.r2.1.1.exome_calling_intervals.sites.vcf.bgz'),
           exomasterfile=os.path.join(global_settings.reference_folder, 'gnomad.exomes.r2.1.1.sites.vcf.bgz'),
           namedexfile=os.path.join(global_settings.dictionary_folder, 'taxid9606-names2uniprot.json'),
           folder=os.path.join(global_settings.temp_folder, 'gnomAD')
           ).split()
    path = os.path.join(global_settings.pickle_folder, f'taxid{taxid}')
    for pf in os.listdir(path):
        try:
            protein = ProteinGatherer().load(file=os.path.join(path, pf))
            protein.gnomAD = []
            protein.parse_gnomAD()
            protein.get_PTM()
            protein.dump()
        except:
            pass