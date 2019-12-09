import pickle, os, re, json
from datetime import datetime
from .settings_handler import global_settings #the instance not the class.
from collections import namedtuple
import gzip, requests

from warnings import warn
from .metadata_from_PDBe import PDBMeta

Variant = namedtuple('Variant', ['id', 'x', 'y', 'impact', 'description', 'homozygous'], defaults=(None, None, None, None, None, None))
Variant.__doc__="""
Stores the gnomAD data for easy use by FeatureViewer and co. Can be converted to Mutation.
"""

class Structure:
    #lolz. a C++ coder would hate this name. Sturcture as in "protein structure"
    #that is not funny. Why I did I think it was?
    #Why am I talking to my past self?!
    """
    No longer a namedtuple.
    Stores the structural data for easy use by FeatureViewer and co. Can be converted to StructureAnalyser
    type = rcsb | swissmodel | homologue
    """
    settings = global_settings

    #__slots__ = ['id', 'description', 'x', 'y', 'url','type','chain','offset', 'coordinates', 'extra']
    def __init__(self, id, description, x:int, y:int, code, type='rcsb',chain='A',offset:int=0, coordinates=None, extra=None, url=''):
        """
        Stores the structural data for easy use by FeatureViewer and co. Can be converted to StructureAnalyser
        type = rcsb | swissmodel | homologue | www | local
        """
        self.id = id
        self.description = description
        self.x = int(x)  # resi in the whole uniprot protein
        self.y = int(y)  # end resi in the whole uniprot protein
        self.offset = int(offset) # offset is the number *subtracted* from the PDB index to make it match the position in Uniprot.
        self.pdb_start = None  # no longer used. TO be deleted.
        self.pdb_end = None   # ditto.
        self.resolution = 0
        self.code = code
        self.chain_definitions = None #filled by SIFTS
        self.type = type.lower()
        self.chain = chain
        if extra is None:
            self.extra = {}
        else:
            self.extra = extra
        self.coordinates = coordinates
        self.url = url  ## for type = www or local or swissmodel
        # https://files.rcsb.org/download/{self.code}.pdb does not work (often) while the url is something odd.

    def to_dict(self):
        return {'x': self.x, 'y': self.y, 'id': self.id, 'description': self.description}

    def __str__(self):
        return str(self.to_dict())

    def get_coordinates(self):
        if self.type == 'rcsb':
            r = requests.get(f'https://files.rcsb.org/download/{self.code}.pdb')
        elif self.type == 'swissmodel':
            r = requests.get(self.url)
        elif self.type == 'www':
            r = requests.get(self.url)
        elif self.type == 'local':
            self.coordinates = open(self.url).read()
            return self.coordinates
        else:
            warn(f'Model type {self.type}  for {self.id} could not be recognised.')
            return None
        if r.status_code == 200:
            self.coordinates = r.text
        else:
            warn(f'Model {self.code} failed.')
        return self.coordinates

    def includes(self, position, offset=0):
        """
        Generally there should not be an offset as x and y are from Uniprot data so they are already fixed!
        :param position:
        :param offset:
        :return:
        """
        if self.x + offset > position:
            return False
        elif self.y + offset < position:
            return False
        else:
            return True


    def lookup_sifts(self):
        """
        SIFTS data. for PDBe query see elsewhere.
        There are four start/stop pairs that need to be compared to get a good idea of a protein.
        For a lengthy discussion see https://blog.matteoferla.com/2019/09/pdb-numbering-rollercoaster.html
        Also for a good list of corner case models see https://proteopedia.org/wiki/index.php/Unusual_sequence_numbering
        :return: self
        """
        def get_offset(detail):
            if detail['PDB_BEG'] == 'None':
                # assuming 1 is the start, which is pretty likely.
                b = int(detail['RES_BEG'])
                if b != 1:
                    warn('SP_BEG is not 1, yet PDB_BEG is without a crystallised start')
            else:
                r = re.search('(-?\d+)', detail['PDB_BEG'])
                if r is None:
                    return self
                b = int(r.group(1))
            return int(detail['SP_BEG']) - b

        if self.type != 'rcsb':
            return self
        details = self._get_sifts()
        ## get matching chain.
        try:
            detail = next(filter(lambda x: self.chain == x['CHAIN'], details))
        except StopIteration:
            warn(f'{self.code} {self.chain} not in {details}')
            return self
        self.offset = get_offset(detail)
        self.chain_definitions = [{'chain': d['CHAIN'],
                                   'uniprot': d['SP_PRIMARY'],
                                   'x': int(d["SP_BEG"]),
                                   'y': int(d["SP_END"]),
                                   'offset': get_offset(d),
                                   'range': f'{d["SP_BEG"]}-{d["SP_END"]}',
                                   'name': None,
                                   'description': None} for d in details]
        return self

    def _get_sifts(self, all_chains=True): #formerly called .lookup_pdb_chain_uniprot
        details = []
        headers = 'PDB     CHAIN   SP_PRIMARY      RES_BEG RES_END PDB_BEG PDB_END SP_BEG  SP_END'.split()
        with self.settings.open('pdb_chain_uniprot') as fh:
            for row in fh:
                if self.code.lower() == row[0:4]:
                    entry = dict(zip(headers, row.split()))
                    if self.chain == entry['CHAIN'] or all_chains:
                        details.append(entry)
        return details

    def lookup_resolution(self):
        if self.type != 'rcsb':
            return self
        with self.settings.open('resolution') as fh:
            resolution = json.load(fh)
            for entry in resolution:
                if entry['IDCODE'] == self.code:
                    if entry['RESOLUTION'].strip():
                        self.resolution = float(entry['RESOLUTION'])
                    break
            else:
                warn(f'No resolution info for {self.code}')
        return self

    def lookup_ligand(self):
        warn('TEMP! Returns the data... not self')
        return PDBMeta(self.code+'_'+self.chain).data


class ProteinCore:
    """
    This is a lightweight version of Protein that is intended to run off pre parsed pickles.
    It forms the base of Protein. This does no protein analyses.
    It has IO powers though .dump/.gdump saves an instance .load/.gload loads and can work as a class method if the filename is provided as an argument.
    The gzipped forms (.gdump and .gload) are about 1/3 the size. 50 KB.
    """
    settings = global_settings
    version = 1.0 #this is for pickled file migration/maintenance.

    def __init__(self, gene_name='', uniprot = '', uniprot_name = '', sequence='', organism = None, taxid=None, **other):
        ### predeclaration (and cheatsheet)
        if organism: # dictionary with keys common scientific and NCBI Taxonomy
            self.organism = organism
        else:
            self.organism = {'common': 'NA', 'scientific': 'NA', 'NCBI Taxonomy': 'NA', 'other': 'NA'} ##obs? ignore for human purposes.
        if taxid:
            self.organism['NCBI Taxonomy'] = taxid
        self.gene_name = gene_name
        self.uniprot_name = uniprot_name.strip() ## S39AD_HUMAN
        #### uniprot derivved
        self.uniprot = uniprot.strip() ## uniprot accession
        self.uniprot_dataset = '' ## Swiss-Prot good, TrEMBL bad.
        self.alt_gene_name_list = []
        self.accession_list = [] ## Q96H72 etc.
        self.sequence = sequence  ###called seq in early version causing eror.rs
        self.recommended_name = '' #Zinc transporter ZIP13
        self.alternative_fullname_list = []
        self.alternative_shortname_list = []
        self.properties={}
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
        self.gnomAD = [] #formerlly alleles
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
        self.timestamp = datetime.now()
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
    def _get_species_folder(self):
        if self.organism['NCBI Taxonomy'] == 'NA':
            species = f'taxid{self.get_species_for_uniprot()}'
        else:
            species = f'taxid{self.organism["NCBI Taxonomy"]}'
        path = os.path.join(self.settings.pickle_folder, species)
        if not os.path.exists(path):
            os.mkdir(path)
        return path


    def exists(self, file=None):
        """
        Method to check if file exists already.
        Actually loads it sneakily!
        :return:
        """
        if file is not None:
            return os.path.exists(file)
        try:
            path = self._get_species_folder()
        except ValueError:
            return False
        for extension, loader in (('.p', self.load), ('.pgz', self.gload)):
            file = os.path.join(path, self.uniprot + extension)
            if os.path.exists(file):
                loader()
                return True
        return False

    def dump(self, file=None):
        if not file:
            path = self._get_species_folder()
            file = os.path.join(path, '{0}.p'.format(self.uniprot))
        self.complete()  # wait complete.
        pickle.dump(self.__dict__, open(file, 'wb'))
        self.log('Data saved to {} as pickled dictionary'.format(file))

    def gdump(self, file=None):
        if not file:
            path = self._get_species_folder()
            file = os.path.join(path,  f'{self.uniprot}.pgz')
        self.complete()  # wait complete.
        if not os.path.exists(path):
            os.mkdir(path)
        with gzip.GzipFile(file, 'w') as f:
            pickle.dump(self.__dict__, f)
        self.log('Data saved to {} as gzipped pickled dictionary'.format(file))

    def get_species_for_uniprot(self):
        warn('You have triggered a fallback. If you know your filepath to load use it.')
        uniprot2species = json.load(open(os.path.join(self.settings.dictionary_folder, 'uniprot2species.json')))
        if self.uniprot in uniprot2species.keys():
            return uniprot2species[self.uniprot]
        else:
            raise ValueError('Cannot figure out species of uniprot to load it. Best bet is to fetch it.')

    #decorator /fake @classmethod
    def _ready_load(fun):
        """
        Prepare loading for both load and gload.
        Formerly allowed it to run as a class method, code not fixed.
        :return:
        """
        def loader(self, file=None):
            if not file:
                path = self._get_species_folder()
                if fun.__name__ == 'load':
                    extension = '.p'
                else:
                    extension = '.pgz'
                file = os.path.join(path, self.uniprot+extension)
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
        Make sure that all subthreads are complete. Not used for Core!
        """
        for k in self._threads:
            if self._threads[k] and self._threads[k].is_alive():
                self._threads[k].join()
        self._threads = {}
        return self
