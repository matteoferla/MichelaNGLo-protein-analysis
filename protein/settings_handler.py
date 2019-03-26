__description__ = """
This is the handler for the settings to control where to save stuff, etc.
It allows customisation of output if the script is not running on a server.
The key parts are:




Note that the folder pages (.pages_folder) was for when it was not for a server. .wipe_html() clears them.
"""
################## Environment ###########################

import os
from pprint import PrettyPrinter

#these are needed for reference file retrieval
import urllib, gzip, shutil, tarfile

pprint = PrettyPrinter().pprint
from warnings import warn

class GlobalSettings:
    """
    This class is container for the paths, which are used by both Variant and Tracker classes.
    Hence why in these two is the attribute .settings
    """
    verbose = False
    subdirectory_names = ('reference', 'temp', 'uniprot','pdbblast', 'pickle', 'binders')

                          #'manual', 'transcript', 'protein', 'uniprot', 'pfam', 'pdb', 'ELM', 'ELM_variant', 'pdb_pre_allele', 'pdb_post_allele', 'ExAC', 'pdb_blast', 'pickle', 'references', 'go',
                          #'binders')
    fetch = True
    missing_attribute_tolerant = True
    error_tolerant = False
    addresses = ('ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz',
                 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/pdbaa.tar.gz',
                 'ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/functional_gene_constraint/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt')

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

    def degunk(self):
        """
        Removes the zero sized files that may ahve arised from error or keyboard interuption.
        :return:
        """
        for dir in self.data_subdirectories:
            for file in os.listdir(dir):
                if os.stat(os.path.join(dir, file)).st_size < 100 and not os.path.isdir(
                        os.path.join(dir, file)):
                    if self.verbose: print('Removing file {}'.format(file))
                    os.remove(os.path.join(dir, file))
        if self.verbose: print('clean-up complete')

    def wipe_html(self):
        """
        No longer needed.
        :return:
        """
        for file in os.listdir(self.page_folder):
            if '.htm' in file or '.pdb' in file:
                os.remove(os.path.join(self.page_folder, file))

    def retrieve_references(self, ask = True, issue = ''):
        if ask:
            print('*' * 20)
            print('CORE reference DATA IS MISSING '+issue)
            print('CORE reference DATA IS MISSING')
            print('There are two options, you have never ever run this script before or the folder {0} is not corrent'.format(self.reference_folder))
            print('this is super experimental (i.e. I\'ve never bother)')
            i = input('Continue y/[n] _')
            if not i or i in ('N', 'n'):
                print('Exiting...')
                exit()
        for url in self.addresses:
            file = os.path.join(self.reference_folder, os.path.split(url)[1])
            if os.path.isfile(file):
                if self.verbose:
                    print('{0} file is present already'.format(file))
            else:
                if self.verbose:
                    print('{0} file is being downloaded'.format(file))
                self._get_url(url, file)
            #raise NotImplementedError('Due to crappy windows 8 computer... this part is manual in VM: cat *.psi > cat.psi where psi files are from http://interactome.baderlab.org/data/')

    def _get_url(self, url, file):
        req = urllib.request.Request(url)
        response = urllib.request.urlopen(req)
        data = response.read()
        with open(file, 'wb') as w:
            w.write(data)

    def _open_reference(self, file):
        fullfile = os.path.join(self.reference_folder, file)
        if not os.path.isfile(fullfile):
            self.retrieve_references(issue = fullfile)
        ## handle compression
        unfile = os.path.join(self.temp_folder, file.replace('.gz','').replace('.tar',''))
        if '.tar.gz' in file and not os.path.exists(unfile):
            tar = tarfile.open(file)
            tar.extractall()
            tar.close()
        elif '.gz' in file and not os.path.isfile(unfile):
            if self.verbose:
                print('{0} file is being temporily extracted to {1}'.format(file, unfile))
                exit(69)
            with open(unfile, 'wb') as f_out:
                with gzip.open(file, 'rb') as f_in:
                    shutil.copyfileobj(f_in, f_out)
        elif '.gz' in file:
            if self.verbose: print('{0} file is already decompressed'.format(file))
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
                'nextprot':'nextprot_refseq.txt',
                'swissmodel':'swissmodel_index.json'}
        assert kind in kdex, 'This is weird. unknown kind, should be: {0}'.format(list(kdex.keys()))
        return self._open_reference(kdex[kind])

global_settings = GlobalSettings()










