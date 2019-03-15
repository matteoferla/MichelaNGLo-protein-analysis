__description__ = """
This is the handler for the settings to control where to save stuff, etc.
It allows customisation of output if the script is not running on a server.
"""
################## Environment ###########################

import os
from pprint import PrettyPrinter

import requests

pprint = PrettyPrinter().pprint
from warnings import warn

class GlobalSettings:
    """
    This class is container for the paths, which are used by both Variant and Tracker classes.
    Hence why in these two is the attribute .settings
    """
    verbose = False
    subdirectory_names = ('manual', 'transcript', 'protein', 'uniprot', 'pfam', 'pdb', 'ELM', 'ELM_variant', 'pdb_pre_allele', 'pdb_post_allele', 'ExAC', 'pdb_blast', 'pickle', 'references', 'go',
                          'binders')
    fetch = True
    missing_attribute_tolerant = True
    error_tolerant = False

    # getter of data_folder
    def _get_datafolder(self):
        return self._datafolder

    # setter of data_folder
    def _set_datafolder(self, new_folder):
        self.data_subdirectories = []
        self._datafolder = new_folder
        if not os.path.isdir(new_folder):
            os.mkdir(new_folder)
        for directory in self.subdirectory_names:
            if new_folder:
                path = os.path.join(new_folder, directory)
            else:
                path = directory
                warn('Setting the data directory to the base directory is a stupid idea.')
            self.data_subdirectories.append(path)
            setattr(self, directory + '_folder', path)
            if not os.path.isdir(path):
                os.mkdir(path)

    data_folder = property(_get_datafolder, _set_datafolder)

    def __init__(self, data_folder='data', page_folder='pages', home_url='/'):
        self.data_folder = data_folder
        self.page_folder = page_folder
        self.home_url = home_url
        self._obodict={}

    def get_folder_of(self, name):
        return getattr(self, name + '_folder')

    @classmethod
    def degunk(cls, verbose=False):
        for dir in cls.data_subdirectories:
            for file in os.listdir(dir):
                if os.stat(os.path.join(dir, file)).st_size < 100 and not os.path.isdir(
                        os.path.join(dir, file)):
                    if verbose: print('Removing file {}'.format(file))
                    os.remove(os.path.join(dir, file))
        if verbose: print('clean-up complete')

    @classmethod
    def wipe_html(cls):
        for file in os.listdir(cls.page_folder):
            if '.htm' in file or '.pdb' in file:
                os.remove(os.path.join(cls.page_folder, file))

    def retrieve_references(self):
        print('*' * 20)
        print('CORE REFERENCE DATA IS MISSING')
        print('There are two options, you have never ever run this script before or the folder {0} is not corrent'.format(self.reference_folder))
        print('this is super experimental (i.e. I\'ve never bother')
        i = input('Continue y/[n] _')
        if not i or i in ('N', 'n'):
            print('Exiting...')
            exit()
        addresses = ('ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/functional_gene_constraint/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt',
                     'ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/ExAC.r1.sites.vep.vcf.gz',
                     'http://geneontology.org/gene-associations/goa_human.gaf.gz',
                     'http://purl.obolibrary.org/obo/go.obo',
                     'http://interactome.baderlab.org/data/Raul-Vidal(Nature_2005).psi',
                     'http://slorth.biochem.sussex.ac.uk/download/h.sapiens_ssl_predictions.csv',
                     'https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-3.5.166/BIOGRID-ALL-3.5.166.mitab.zip',
                     'https://stringdb-static.org/download/protein.links.v10.5/9606.protein.links.v10.5.txt.gz',
                     'http://www.ensembl.org/biomart/martview/436b8b2b06f64bbee960592afda10817?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id|hsapiens_gene_ensembl.default.feature_page.ensembl_transcript_id|hsapiens_gene_ensembl.default.feature_page.external_gene_name|hsapiens_gene_ensembl.default.feature_page.uniprotswissprot|hsapiens_gene_ensembl.default.feature_page.uniprot_gn|hsapiens_gene_ensembl.default.feature_page.ensembl_peptide_id&FILTERS=&VISIBLEPANEL=resultspanel',
                     'ftp://ftp.nextprot.org/pub/current_release/mapping/nextprot_refseq.txt')
        for url in addresses:
            file = os.path.split(url)[1]
            if os.path.isfile(os.path.join(self.reference_folder, file)):
                continue
            else:
                req = requests.get(url)  # headers={"Accept": "application/xml"}
                if req.status_code != 200:
                    raise ConnectionError('Could not retrieve data: ' + req.text)
                data = req.text
                open(file, 'w').write(data)
                if os.path.splitext(url) == 'gz':
                    raise Exception('Okay. this will not work on windows... Can you unzip this file?? cd {0}; tar -x {1}; cd ../..'.format(self.reference_folder, file))
            raise NotImplementedError('Due to crappy windows 8 computer... this part is manual in VM: cat *.psi > cat.psi where psi files are from http://interactome.baderlab.org/data/')

    def _open_reference(self, file):
        fullfile = os.path.join(self.references_folder, file)
        if not os.path.isfile(fullfile):
            self.retrieve_references()
        return open(fullfile)

    def open(self, kind):
        kdex = {'ExAC_pLI': 'fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt',
                'ExAC_vep': 'ExAC.r1.sites.vep.vcf',
                'ID_mapping': 'HUMAN_9606_idmapping_selected.tab',
                'ssl': 'h.sapiens_ssl_predictions.csv',
                'go': 'go.obo',
                'go_human': 'goa_human.gaf',
                'huri':'cat.psi',
                'biogrid':'BIOGRID-ALL-3.5.166.mitab.txt',
                'string':'9606.protein.links.v10.5.txt',
                'ensembl':'ensemb.txt',
                'nextprot':'nextprot_refseq.txt'}
        assert kind in kdex, 'This is weird. unknown kind, should be: {0}'.format(list(kdex.keys()))
        return self._open_reference(kdex[kind])

global_settings = GlobalSettings()










