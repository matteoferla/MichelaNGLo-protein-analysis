import pickle, os
from datetime import datetime
from protein.settings_handler import global_settings #the instance not the class.

class ProteinLite:
    """
    This is a lightweight version of Protein that is intended to run of pre parsed pickles.
    It forms the base of Protein.
    """
    settings = global_settings

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
        self.ExAC_type = 'Unparsed' # Dominant | Recessive | None | Unknown (=???)
        self.pLI = -1
        self.pRec = -1
        self.pNull = -1
        ### pdb
        self.pdb_matches =[] #{'match': align.title[0:50], 'match_score': hsp.score, 'match_start': hsp.query_start, 'match_length': hsp.align_length, 'match_identity': hsp.identities / hsp.align_length}
        self.swissmodel = []
        self.percent_modelled = 0
        ### other ###
        self.user_text = ''
        ### mutation ###
        self.mutation = None
        ### junk
        self.other = other ### this is a garbage bin. But a handy one.
        self.logbook = [] # debug purposes only. See self.log()
        self._threads = {}
        #not needed for ProteinLite
        self.xml = None

    ############################# IO #############################
    def dump(self, file=None):
        if not file:
            file = os.path.join(self.settings.pickle_folder, '{0}.p'.format(self.uniprot))
        self.complete()  # wait complete.
        pickle.dump(self.__dict__, open(file, 'wb'))
        self.log('Data saved to {} as pickled dictionary'.format(file))

    #@classmethod
    def load(clelf, file=None):  # clelf made-up portmanteau of cls and self, non PEP
        if isinstance(clelf, type):  # clelf is a class
            self = clelf.__new__(clelf)
            assert file, 'file is mandatory for `Protein.load()` as a class method. optional as bound method.'
        else:
            self = clelf
            if not file:
                file = os.path.join(self.settings.pickle_folder,self.uniprot+'.p')
        self.__dict__ = pickle.load(open(file, 'rb'))
        self.log('Data from the pickled dictionary {}'.format(file))
        return self

    ####################### Misc Magic methods ##################
    def __len__(self):  ## sequence lenght
        return len(self.sequence)

    def log(self, text):
        msg = '[{}]\t'.format(str(datetime.now())) + text
        self.logbook.append(msg)
        if self.settings.verbose:
            print(msg)

    def __str__(self):
        if self.gene_name:
            return self.gene_name
        else:
            return self.uniprot

    def complete(self):
        """
        Make sure that all subthreads are complete. Not used for Lite.
        """
        for k in self._threads:
            if self._threads[k] and self._threads[k].is_alive():
                self._threads[k].join()
        self._threads = {}
        return self


    ################### mutant related
    def predict_effect(self):
        pass
        ##TOdo write>!

    def check_mutation(self, mutation):
        if len(self.sequence) > mutation.residue_index and self.sequence[mutation.residue_index - 1] == mutation.from_residue:
            return True
        else:
            return False # call mutation_discrepancy to see why.

    def mutation_discrepancy(self, mutation):
        # returns a string explaining the `check_mutation` discrepancy error
        neighbours=''
        if len(self.sequence) < mutation.residue_index:
            return 'Uniprot {g} is {l} amino acids long, while user claimed a mutation at {i}.'.format(
                g=self.uniprot,
                i=mutation.residue_index,
                l=len(self.sequence)
                )
        else:
            if mutation.residue_index < 5:
                neighbours = '{pre}*{i}*{post}'.format(pre=self.sequence[:mutation.residue_index-1],
                                                       i=self.sequence[mutation.residue_index-1],
                                                       post=self.sequence[mutation.residue_index:mutation.residue_index+5])
            elif mutation.residue_index > 5 and len(self.sequence) > mutation.residue_index + 5:
                neighbours = '{pre}*{i}*{post}'.format(pre=self.sequence[mutation.residue_index - 6:mutation.residue_index - 1],
                                                       i=self.sequence[mutation.residue_index - 1],
                                                       post=self.sequence[mutation.residue_index:mutation.residue_index + 5])
            elif len(self.sequence) < mutation.residue_index + 5:
                neighbours = '{pre}*{i}*{post}'.format(pre=self.sequence[mutation.residue_index - 6:mutation.residue_index - 1],
                                                       i=self.sequence[mutation.residue_index - 1],
                                                       post=self.sequence[mutation.residue_index:])
            else:
                print('impossible?!')
                neighbours = 'ERROR.'
            return 'Residue {i} is {n} in Uniprot {g}, while user claimed it was {f}. (neighbouring residues: {s}'.format(
                                                                            i=mutation.residue_index,
                                                                            n=self.sequence[mutation.residue_index - 1],
                                                                            g=self.uniprot,
                                                                            f=mutation.from_residue,
                                                                            s=neighbours
                                                                            )

    def get_features_at_position(self, position):
        return self.get_features_near_position(position,wobble=0)

    def get_features_near_position(self, position, wobble=10):
        ## dodgy position?
        if isinstance(position,str):
            position = int(position)
        elif not isinstance(position,int):
            position = position.residue_index
        else:
            position = position
        ## deal with it.
        valid = [{**f,'type': g} for g in self.features for f in self.features[g] if f['x'] - wobble < position < f['y'] + wobble]
        return valid

