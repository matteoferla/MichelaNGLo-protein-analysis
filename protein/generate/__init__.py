"""
Generates the data to be used.
"""

from pprint import PrettyPrinter
pprint = PrettyPrinter().pprint

import os, json
from protein import Protein
from .uniprot_to_jsons import UniprotReader
from .PDB_blast import Blaster

def generate():
    ################ FETCH ALL RAW FILES #################################################
    Protein.settings.verbose = True
    #Protein.settings.retrieve_references(ask=False)

    ################ PARSE UNIPROT #######################################################
    master_file = os.path.join(Protein.settings.reference_folder, 'uniprot_sprot.xml')
    # This class parses the uniprot FTP file and can do various things. such as making a small one that is only human.
    # But mainly the `UniprotReader.convert('uniprot_sprot.xml')` method whcih generates the JSON files required.
    # first_n_protein is for testing.
    #UniprotReader.convert(uniprot_master_file = master_file, first_n_protein=0)

    ################ BLAST PDB #######################################################
    # uncompresses the pdbaa
    #Blaster.extract_db()
    # blasts the human.fa agains the newly extracted pdbaa
    #Blaster.pdb_blaster()
    # converts the files to something reasonable.
    # Blaster.parse('blastpdb','blastpdb2')
    parse_proteome()

def parse_proteome():
    Protein.settings.fetch = False
    Protein.settings.error_tolerant = False
    Protein.settings.missing_attribute_tolerant = False
    Protein.settings.verbose = True
    uniprot_ids = set(json.load(open(os.path.join(Protein.settings.data_folder,'human_prot_namedex.json'))).values())
    for acc in uniprot_ids:
        prot = Protein(uniprot=acc)
        prot.parse_uniprot()
        prot.parse_pLI()
        prot.parse_swissmodel()
        prot.parse_gNOMAD()
        prot.get_percent_modelled()
        prot.parse_ExAC_type()
        prot.parse_pdb_blast()
        prot.fetch_binders()
        pprint({k: len(prot.partners[k]) for k in prot.partners})
        prot.dump()



refs=(
     'ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/ExAC.r1.sites.vep.vcf.gz',
     'http://geneontology.org/gene-associations/goa_human.gaf.gz',
     'http://purl.obolibrary.org/obo/go.obo',
     'http://interactome.baderlab.org/data/Raul-Vidal(Nature_2005).psi',
     'http://slorth.biochem.sussex.ac.uk/download/h.sapiens_ssl_predictions.csv',
     'https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-3.5.166/BIOGRID-ALL-3.5.166.mitab.zip',
     'https://stringdb-static.org/download/protein.links.v10.5/9606.protein.links.v10.5.txt.gz',
     'http://www.ensembl.org/biomart/martview/436b8b2b06f64bbee960592afda10817?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id|hsapiens_gene_ensembl.default.feature_page.ensembl_transcript_id|hsapiens_gene_ensembl.default.feature_page.external_gene_name|hsapiens_gene_ensembl.default.feature_page.uniprotswissprot|hsapiens_gene_ensembl.default.feature_page.uniprot_gn|hsapiens_gene_ensembl.default.feature_page.ensembl_peptide_id&FILTERS=&VISIBLEPANEL=resultspanel',
     'ftp://ftp.nextprot.org/pub/current_release/mapping/nextprot_refseq.txt',
     'https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.exome_calling_intervals.sites.vcf.bgz')
