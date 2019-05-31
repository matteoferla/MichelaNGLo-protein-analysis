import pickle, os, re
from datetime import datetime
from .settings_handler import global_settings #the instance not the class.
from collections import namedtuple
import gzip

Variant = namedtuple('Variant', ['id', 'x', 'y', 'impact', 'description'])
Variant.__doc__="""
Stores the gNOMAD data for easy use by FeatureViewer and co. Can be converted to Mutation.
"""

class Structure: #a C++ coder would hate this... Sturcture as in protein structure
    """
    No longer a namedtuple.
    Stores the structural data for easy use by FeatureViewer and co. Can be converted to StructureAnalyser
    type = rcsb | swissmodel | homologue
    """
    __slots__ = ['id', 'description', 'x', 'y', 'url','type','chain','offset', 'coordinates', 'extra']
    def __init__(self, id, description, x:int, y:int, url, type='rcsb',chain='A',offset:int=0, coordinates=None, extra=None):
        """
        Stores the structural data for easy use by FeatureViewer and co. Can be converted to StructureAnalyser
        type = rcsb | swissmodel | homologue
        """
        self.id = id
        self.description = description
        self.x = int(x)
        self.y = int(y)
        self.url = url
        self.type = type
        self.chain = chain
        self.offset = int(offset)
        self.extra = extra
        self.coordinates = coordinates


class ProteinCore:
    """
    This is a lightweight version of Protein that is intended to run of pre parsed pickles.
    It forms the base of Protein. This does no protein analyses.
    It has IO powers! .dump/.gdump saves an instance .load/.gload loads and can work as a class method if the filename is provided as an argument.
    The gzipped forms (.gdump and .gload) are about 1/3 the size. 50 KB.
    """
    settings = global_settings
    version = 1.0 #this is for pickled file migration/maintenance.


    ############## elm
    elmdata = []
    with open(os.path.join(settings.data_folder, 'elm_classes.tsv')) as fh:
        header = ("Accession", "ELMIdentifier", "FunctionalSiteName", "Description", "Regex", "Probability", "#Instances", "#Instances_in_PDB")
        for line in fh:
            if line[0] == '#':
                continue
            if "Accession" in line:
                continue
            elmdata.append(dict(zip(header, line.replace('"','').split('\t'))))

    def __init__(self, gene_name='', uniprot = '', uniprot_name = '', sequence='', **other):
        ### predeclaration (and cheatsheet)
        self.gene_name = gene_name
        self.uniprot_name = uniprot_name ## S39AD_HUMAN
        #### uniprot derivved
        self.uniprot = uniprot ## uniprot accession
        self.alt_gene_name_list = []
        self.accession_list = [] ## Q96H72 etc.
        self.sequence = sequence  ###called seq in early version causing eror.rs
        self.recommended_name = '' #Zinc transporter ZIP13
        self.alternative_fullname_list = []
        self.alternative_shortname_list = []
        self.features={}  #see _parse_protein_feature. Dictionary of key: type of feature, value = list of dict with the FeatureViewer format (x,y, id, description)
        self.partners ={'interactant': [],  #from uniprot
                        'BioGRID': [],  #from biogrid downlaoad
                        'SSL': [],  #Slorth data
                        'HuRI': [],
                        'stringDB highest': [],  # score >900
                        'stringDB high': [],  #900 > score > 700
                        'stringDB medium': [], #400 > score > 400
                        'stringDB low': [] #score < 400
                        } # lists not sets as it gave a pickle issue.
        self.diseases=[] # 'description', 'name', 'id', 'MIM'
        self.pdbs = []  # {'description': elem.attrib['id'], 'id': elem.attrib['id'], 'x': loca[0], 'y': loca[1]}
        self.ENSP = ''
        self.ENST = ''
        self.ENSG = ''
        ### ExAC
        self.gNOMAD = [] #formerlly alleles
        #self.ExAC_type (property= 'Unparsed' # Dominant | Recessive | None | Unknown (=???)
        self.pLI = -1
        self.pRec = -1
        self.pNull = -1
        ### pdb
        self.pdb_matches =[] #{'match': align.title[0:50], 'match_score': hsp.score, 'match_start': hsp.query_start, 'match_length': hsp.align_length, 'match_identity': hsp.identities / hsp.align_length}
        self.swissmodel = []
        self.percent_modelled = -1
        ### junk
        self.other = other ### this is a garbage bin. But a handy one.
        self.logbook = [] # debug purposes only. See self.log()
        self._threads = {}
        #not needed for ProteinLite
        self.xml = None

    ############################## property objects

    @property
    def ExAC_type(self):
        if self.pLI < 0:  # error.
            return 'Unknown'
        elif self.pLI > max(self.pRec, self.pNull):
            return 'Dominant'
        elif self.pRec > max(self.pLI, self.pNull):
            return 'Recessive'
        elif self.pNull > max(self.pLI, self.pRec):
            return 'None'
        else:
            return 'Unknown'

    ############################# IO #############################
    def dump(self, file=None):
        if not file:
            file = os.path.join(self.settings.pickle_folder, '{0}.p'.format(self.uniprot))
        self.complete()  # wait complete.
        pickle.dump(self.__dict__, open(file, 'wb'))
        self.log('Data saved to {} as pickled dictionary'.format(file))

    def gdump(self, file=None):
        if not file:
            file = os.path.join(self.settings.pickle_folder, f'{self.uniprot}.pgz')
        self.complete()  # wait complete.
        with gzip.GzipFile(file, 'w') as f:
            pickle.dump(self.__dict__, f)
        self.log('Data saved to {} as gzipped pickled dictionary'.format(file))

    #decorator /fake @classmethod
    def _ready_load(fun):
        """
        Prepare loading for both load and gload. Deals with the fact that load can be used as a class method or a bound method.
        :return:
        """
        def loader(clelf,file=None): # clelf made-up portmanteau of cls and self, non PEP
            if isinstance(clelf, type):  # clelf is a class
                self = clelf.__new__(clelf)
                assert file, 'file is mandatory for `Protein.load()` as a class method. optional as bound method.'
            else:
                self = clelf
                if not file:
                    if fun.__name__ == 'load':
                        extension = '.p'
                    else:
                        extension = '.pgz'
                    file = os.path.join(self.settings.pickle_folder, self.uniprot+extension)
            fun(self, file)
            return self
        return loader

    @_ready_load
    def load(self, file):
        self.__dict__ = pickle.load(open(file, 'rb'))
        self.log('Data from the pickled dictionary {}'.format(file))
        return self

    @_ready_load
    def gload(self, file):
        with gzip.GzipFile(file, 'r') as f:
            self.__dict__ = pickle.load(f)
        self.log('Data from the gzipped pickled dictionary {}'.format(file))
        return self

    ####################### Misc Magic methods ##################
    def __len__(self):  ## sequence lenght
        return len(self.sequence)

    def log(self, text):
        """
        Logging is primarily for protein_full
        :param text:
        :return:
        """
        msg = '[{}]\t'.format(str(datetime.now())) + text
        self.logbook.append(msg)
        if self.settings.verbose:
            print(msg)
        return self

    def __str__(self):
        if len(self.gene_name):
            return self.gene_name
        else:
            return self.uniprot

    def complete(self):
        """
        Make sure that all subthreads are complete. Not used for Lite!
        """
        for k in self._threads:
            if self._threads[k] and self._threads[k].is_alive():
                self._threads[k].join()
        self._threads = {}
        return self
