import os
from protein import Protein
from .uniprot_to_jsons import UniprotReader

def generate():
    ################ FETCH ALL RAW FILES #################################################
    Protein.settings.verbose = True
    Protein.settings.retrieve_references(ask=False)

    ################ PARSE UNIPROT #######################################################
    master_file = os.path.join(Protein.settings.temp_folder, 'uniprot_sprot.xml')
    # This class parses the uniprot FTP file and can do various things. such as making a small one that is only human.
    # But mainly the `UniprotReader.convert('uniprot_sprot.xml')` method whcih generates the JSON files required.
    # first_n_protein is for testing.
    #UniprotReader.convert(uniprot_master_file = master_file, first_n_protein=0)

refs=(
     'ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/ExAC.r1.sites.vep.vcf.gz',
     'http://geneontology.org/gene-associations/goa_human.gaf.gz',
     'http://purl.obolibrary.org/obo/go.obo',
     'http://interactome.baderlab.org/data/Raul-Vidal(Nature_2005).psi',
     'http://slorth.biochem.sussex.ac.uk/download/h.sapiens_ssl_predictions.csv',
     'https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-3.5.166/BIOGRID-ALL-3.5.166.mitab.zip',
     'https://stringdb-static.org/download/protein.links.v10.5/9606.protein.links.v10.5.txt.gz',
     'http://www.ensembl.org/biomart/martview/436b8b2b06f64bbee960592afda10817?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id|hsapiens_gene_ensembl.default.feature_page.ensembl_transcript_id|hsapiens_gene_ensembl.default.feature_page.external_gene_name|hsapiens_gene_ensembl.default.feature_page.uniprotswissprot|hsapiens_gene_ensembl.default.feature_page.uniprot_gn|hsapiens_gene_ensembl.default.feature_page.ensembl_peptide_id&FILTERS=&VISIBLEPANEL=resultspanel',
     'ftp://ftp.nextprot.org/pub/current_release/mapping/nextprot_refseq.txt')
