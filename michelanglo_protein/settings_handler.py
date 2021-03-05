"""
The Singleton GlobalSettings is the handler for the settings to control where to save stuff, etc.
It allows customisation of output if the script is not running on a server.
The key parts are:

- startup. It is initialised when the module is imported. However, it is not ready as the files need to be configured with startup.
- retrieve_references. download all the bits. Phosphosite need manuall download. See `licence_note` in `michelanglo_protein.generate.split_phosphosite`.


Note that the folder pages (.pages_folder) was for when it was not for a server. say .wipe_html() clears them.
This is old code.
"""
######### Environment ##############

import os, json
import zipfile
from pprint import PrettyPrinter

#these are needed for reference file retrieval
import requests, gzip, shutil, tarfile
from urllib.request import urlopen

pprint = PrettyPrinter().pprint
from warnings import warn

class Singleton(type): #https://stackoverflow.com/questions/6760685/creating-a-singleton-in-python
    """There can only be one setting."""
    _instances = {}
    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        else:
            warn('Attempt to initialise another instance of a singleton. Returning original.')
        return cls._instances[cls]

## Some reference files have to be fetched manually...
## this is for humans/will need to be updated below.
_refs = (
    'ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/ExAC.r1.sites.vep.vcf.gz',
    'http://geneontology.org/gene-associations/goa_human.gaf.gz',
    'http://purl.obolibrary.org/obo/go.obo',
    'http://interactome.baderlab.org/data/Raul-Vidal(Nature_2005).psi',
    'http://slorth.biochem.sussex.ac.uk/download/h.sapiens_ssl_predictions.csv',
    'https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-3.5.166/BIOGRID-ALL-3.5.166.mitab.zip',
    'https://stringdb-static.org/download/protein.links.v10.5/9606.protein.links.v10.5.txt.gz',
    'http://www.ensembl.org/biomart/martview/436b8b2b06f64bbee960592afda10817?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id|hsapiens_gene_ensembl.default.feature_page.ensembl_transcript_id|hsapiens_gene_ensembl.default.feature_page.external_gene_name|hsapiens_gene_ensembl.default.feature_page.uniprotswissprot|hsapiens_gene_ensembl.default.feature_page.uniprot_gn|hsapiens_gene_ensembl.default.feature_page.ensembl_peptide_id&FILTERS=&VISIBLEPANEL=resultspanel',
    'ftp://ftp.nextprot.org/pub/current_release/mapping/nextprot_refseq.txt',
    'https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.exome_calling_intervals.sites.vcf.bgz'
    'ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz')

# http://www.sbg.bio.ic.ac.uk/~missense3d/download/all_dataset.xlsx
# http://www.sbg.bio.ic.ac.uk/~missense3d/download/1052_hom_str.zip


class GlobalSettings(metaclass=Singleton):
    """
    This class is container for the paths, which are used by both Variant and Tracker classes.
    Hence why in these two is the attribute .settings
    """
    verbose = False #:verbose boolean controls the verbosity of the whole module.
    subdirectory_names = ('reference', 'temp', 'uniprot', 'pdbblast', 'pickle', 'binders', 'dictionary')

                          #'manual', 'transcript', 'protein', 'uniprot', 'pfam', 'pdb', 'ELM', 'ELM_variant', 'pdb_pre_allele', 'pdb_post_allele', 'ExAC', 'pdb_blast', 'pickle', 'references', 'go',
                          #'binders')
    fetch = True #: boolean for whether to download data from the interwebs.
    missing_attribute_tolerant = True
    error_tolerant = False
    addresses = ['ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz',
                 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/pdbaa.tar.gz',
                 'ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/functional_gene_constraint/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt',
                 'ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz',
                 'ftp://ftp.wwpdb.org/pub/pdb/derived_data/index/resolu.idx',
                 'https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.exome_calling_intervals.sites.vcf.bgz',
                 'https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz',
                 'http://www.sbg.bio.ic.ac.uk/~missense3d/download/1052_hom_str.zip'] #from http://www.sbg.bio.ic.ac.uk/~missense3d/dataset.html
    for _s in (9606, 3702, 6239, 7227, 10090, 36329, 83332, 83333, 93061, 190650, 208964, 284812, 559292):
        addresses.append(f'https://swissmodel.expasy.org/repository/download/core_species/{_s}_meta.tar.gz')
    manual_task_note = """'# Manual TASKS\nRemember that user has manually downloaded _site_dataset.gz files from https://www.phosphosite.org/staticDownloads at phosphosite.'"""

    # getter of data_folder
    def _get_datafolder(self):
        if not self._initialised:
            self.init()
        return self._datafolder

    # setter of data_folder
    def _set_datafolder(self, new_folder):
        self.data_subdirectories = []
        self._datafolder = new_folder
        if not os.path.isdir(new_folder):
            os.mkdir(new_folder)
            warn(f'The folder {new_folder} was created. This means that you will have to run `global_settings.retrieve_references(ask=False, refresh=True)`')
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

    def __init__(self, home_url='/'):
        self._obodict={}
        self.home_url = home_url
        self._initialised = False

    def startup(self, data_folder='data'):
        if self._initialised:
            raise Exception('The module is already initialised.')
        self._initialised = True
        self.data_folder = data_folder
        #self.page_folder = page_folder  # does nothing.
        print(f'Folder path set to {self.data_folder}')
        return self


    def get_folder_of(self, name):
        if not self._initialised:
            self.init()
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

    def retrieve_references(self, ask = True, refresh=False, issue = ''):
        if not self._initialised:
            raise ValueError('You have not initialised the settings (set the folder) >>> run xxx.settings.startup()')
        if ask:
            print('*' * 20)
            print('CORE reference DATA IS MISSING --trigger by '+issue)
            print('There are two options, you have never ever run this script before or the folder {0} is not corrent'.format(self.reference_folder))
            print('this is super experimental (i.e. I\'ve never bother)')
            i = input('Continue y/[n] _')
            if not i or i in ('N', 'n'):
                print('Exiting...')
                exit()
        for url in self.addresses:
            self._deal_w_url(url, refresh)
        # convert dodgy ones.
        self.create_json_from_idx('resolu.idx', 'resolution.json')
        print(self.manual_task_note)

        #implement cat *.psi > cat.psi where psi files are from http://interactome.baderlab.org/data/')

    def _deal_w_url(self, url, refresh=False) -> str:
        """
        If the file does not exist or refresh is true, it downloads (calling``self._get_url(url, file)``) and unzips the page (``self._unzip_file(file)``).
        :param url:
        :param refresh:
        :return: the file name (full)
        """
        file = os.path.join(self.reference_folder, os.path.split(url)[1])
        unfile = file.replace('.gz', '').replace('.tar', '').replace('.zip', '')
        if os.path.isfile(unfile) and not refresh:
            if self.verbose:
                print('{0} unzipped file is present already'.format(unfile))
        elif os.path.isfile(file) and not refresh:
            if self.verbose:
                print('{0} zipped file is present already, but not unzipped'.format(file))
            self._unzip_file(file)
        else:
            if os.path.isfile(unfile):
                os.remove(unfile)
            if self.verbose:
                print('{0} file is being downloaded'.format(file))
            self._get_url(url, file)
            self._unzip_file(file)
        return file

    def _get_url(self, url, file):
        if 'ftp://' in url:
            req = urlopen(url)
            with open(file, 'wb') as fp:
                shutil.copyfileobj(req, fp)
        else:
            with requests.get(url, stream=True) as r:
                with open(file, 'wb') as f:
                    shutil.copyfileobj(r.raw, f)

    def _unzip_file(self, file):
        unfile = file.replace('.gz', '').replace('.tar', '').replace('.zip', '')
        if os.path.exists(unfile):
            if self.verbose:
                print('{0} file has already been extracted to {1}'.format(file, unfile))
            return self
        elif '.tar.gz' in file:
                os.mkdir(unfile)
                tar = tarfile.open(file)
                tar.extractall(path=unfile)
                tar.close()
        elif '.gz' in file:  #ignore the .bgz of gnomAD. it is too big.
                with open(unfile, 'wb') as f_out:
                    with gzip.open(file, 'rb') as f_in:
                        shutil.copyfileobj(f_in, f_out)
                #if self.verbose: print('{0} file is already decompressed'.format(file))
        elif '.zip' in file:
            with zipfile.ZipFile(file, 'r') as zip_ref:
                zip_ref.extractall(self.reference_folder)
        else:
            pass #not a compressed file
        return self

    def _open_reference(self, file, mode='r'):
        fullfile = os.path.join(self.reference_folder, file)
        if mode == 'w':
            return open(fullfile, 'w')
        elif not os.path.isfile(fullfile):
            self.retrieve_references(issue = fullfile)
        # handle compression
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
                'string':'9606.michelanglo_protein.links.v10.5.txt',
                'ensembl':'ensemb.txt',
                'nextprot':'nextprot_refseq.txt',
                'pdb_chain_uniprot': 'pdb_chain_uniprot.tsv',
                'elm':'elm_classes.tsv',
                'resolution': 'resolution.json',
                'ensembl-uniprot': 'Homo_sapiens.GRCh38.99.uniprot.tsv'}
        if 'swissmodel' in kind:
            taxid = kind.replace('swissmodel','')
            if taxid == '':
                taxid = '9606' #legacy.
            return self._open_reference(f'{taxid}_meta/SWISS-MODEL_Repository/INDEX.json')
        else:
            assert kind in kdex, 'This is weird. unknown kind, should be: {0}'.format(list(kdex.keys()))
            return self._open_reference(kdex[kind])

    def create_json_from_idx(self, infile, outfile):
        # resolu.idx is in the weirdest format.
        fh = self._open_reference(infile)
        for row in fh:
            if not row.strip():
                break
        header = next(fh).split()
        next(fh) #dashes
        parts = [dict(zip(header, [f.strip() for f in row.split(';')])) for row in fh if row.strip()]
        json.dump(parts, self. _open_reference(outfile, mode='w'))

global_settings = GlobalSettings()