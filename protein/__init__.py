__docs__ = """
Protein was formerly called Variant.
"""

import os
import pickle
import re
import csv
import json
import time
import markdown
import xmlschema
import threading
from collections import defaultdict
from datetime import datetime
from warnings import warn
from shutil import copyfile
import requests  # for xml fetcher.

from protein.ET_monkeypatch import ET #monkeypatched version
from protein.settings_handler import global_settings
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

from protein._protein_uniprot_mixin import _UniprotMixin
# Protein inherits _UniprotMixin, which in turn inherits _BaseMixin
# `.settings` class attribute is global_settings from settings_handler.py and is added by _BaseMixin.
# _BaseMixin is inherited by _UniprotMixin contains _failsafe decorator, __getattr__ and settings

#######################
class Protein(_UniprotMixin):
    """
    This class handles each protein entry from Uniprot. See __init__ for the list of variables stored.
    It fills them from various other sources.

        >>> Protein()
        >>> Protein(xml_entry) # equivalent to Protein()._parse_uniprot_xml(xml_entry)
        >>> Protein.load('filename')

    NB. The ET.Element has to be monkeypatched. See `help(ElementalExpansion)`
    """

    # these older commands should be made redundant
    #    fetch = True
    #    croak = True
    #    tollerate_no_SNV = True

    # decorator
    def _failsafe(func):
        def wrapper(self, *args, **kargs):
            # the call happned after chekcing if it should croak on error so to make the traceback cleaner.
            if self.settings.error_tollerant:
                try:
                    return func(self, *args, **kargs)
                except Exception as error:
                    print('Error caught in method `Protein().{n}`: {e}'.format(n=func.__name__, e=error))
                    return None
            else:
                return func(self, *args, **kargs)

        return wrapper

    ############################# INIT #############################

    def __init__(self, entry=None, gene_name='', uniprot = '', uniprot_name = '', sequence='', **other):
        ### predeclaration (and cheatsheet)
        self.xml = entry
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
                        'stringDB highest': [],  # score >900
                        'stringDB high': [],  #900 > score > 700
                        'stringDB medium': [], #400 > score > 400
                        'stringDB low': [] #score < 400
                        } # lists not sets as it gave a pickle issue.
        self.diseases=[] # 'description', 'name', 'id', 'MIM'
        self.pdbs = []  # {'description': elem.attrib['id'], 'id': elem.attrib['id'], 'x': loca[0], 'y': loca[1]}
        ### ExAC
        self.ExAC_type = 'Unparsed' # Dominant | Recessive | None | Unknown (=???)
        self.pLI = -1
        self.pRec = -1
        self.pNull = -1
        ### pdb
        self.pdb_matches =[] #{'match': align.title[0:50], 'match_score': hsp.score, 'match_start': hsp.query_start, 'match_length': hsp.align_length, 'match_identity': hsp.identities / hsp.align_length}
        ### other ###
        self.user_text = ''
        ### mutation ###
        self.mutation = None
        ### junk
        self.other = other ### this is a garbage bin. But a handy one.
        self.logbook = [] # debug purposes only. See self.log()
        self._threads = {}
        ### fill
        if entry:
            self._parse_uniprot_xml(entry)

    def complete(self):
        for k in self._threads:
            if self._threads[k] and self._threads[k].is_alive():
                self._threads[k].join()
        self._threads = {}
        return self

    def __len__(self): ## sequence lenght
        return len(self.sequence)


    ############################# IO #############################
    def dump(self, file=None):
        if not file:
            file = os.path.join(self.settings.pickle_folder, '{0}.p'.format(self.uniprot_name))
        self.complete() # wait complete.
        pickle.dump(self.__dict__, open(file, 'wb'))
        self.log('Data saved to {} as pickled dictionary'.format(file))

    @classmethod
    def load(cls, file):
        self = cls.__new__(cls)
        self.__dict__ = pickle.load(open(file, 'rb'))
        self.log('Data from the pickled dictionary {}'.format(file))
        return self

    def write_uniprot(self, file=None):
        if not file:
            file=self.uniprot_name+'.xml'
        with open(file,'w') as w:
            w.write(
                '<?xml version="1.0" encoding="UTF-8"?><uniprot xmlns="http://uniprot.org/uniprot" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" ' +
                'xsi:schemaLocation="http://uniprot.org/uniprot http://www.uniprot.org/docs/uniprot.xsd">')
            if isinstance(self.xml,str):
                w.write(self.xml)
            else:
                w.write(ET.tostring(self.xml).decode())
            w.write('</uniprot>')


    def log(self, text):
        msg = '[{}]\t'.format(str(datetime.now())) + text
        self.logbook.append(msg)
        if self.settings.verbose:
            print(msg)

    ############################# Data gathering #############################
    def _assert_fetchable(self, text):
        if not self.settings.fetch:
            raise FileNotFoundError(
                text + ' failed previously (or never run previously) and `.settings.fetch` is enabled')


    def xml_fetcher(self, mode):  # mode = uniprot or pfam
        file = os.path.join(self.settings.get_folder_of(mode), self.uniprot + '_' + mode + '.xml')
        if os.path.isfile(file):
            with open(file, 'r') as w:
                xml = w.read()
            self.log('{0} read from file {1}'.format(mode, file))
        else:
            warn('Missing file {0}'.format(file))
            self._assert_fetchable(mode)
            if mode == 'uniprot':
                requestURL = 'https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=100&accession={acc}&taxid=9606'.format(
                    acc=self.uniprot)
            elif mode == 'pfam':
                requestURL = 'https://pfam.xfam.org/protein?output=xml&acc={acc}'.format(acc=self.uniprot)
            else:
                raise ValueError('Only options for mode are uniprot or pfam')
            req = requests.get(requestURL, headers={"Accept": "application/xml"})
            if req.status_code != 200:
                raise ConnectionError('Could not retrieve data: ' + req.text)
            xml = req.text
            with open(file, 'w') as w:
                w.write(xml)
        return xml

    def xml_parser(self, mode, xml):
        """
        THIS METHOD WILL BE DEPRACATED ONCE PFAM PARSER IS REWRITTEN!!
        This is weird. uniprot works fine getting it parsed via its XML schema, while pfam fails even with validation = lax.
        So the pfam data uses stackoverflow.com/questions/7684333/converting-xml-to-dictionary-using-elementtree
        and not schema. I favour the latter as it ain't a just a snippet from SO, but the latter seems okay.
        I looked at how https://github.com/prody/ProDy/blob/master/prody/database/pfam.py does it and they
        simply returns the ET root yet say it returns a dictionary. :-S
        :param mode:
        :param xml:
        :return:
        """

        def deep_clean(element, damn):  # the values have the annoying {http} field this removes them.
            try:
                if isinstance(element, dict):
                    return {k.replace(damn, '').replace('@', '').replace('#', ''): deep_clean(element[k], damn) for k in
                            element}
                elif isinstance(element, list):
                    return [deep_clean(e, damn) for e in element]
                elif isinstance(element, str):
                    if re.fullmatch('\d*', element):
                        return int(element)
                    elif re.fullmatch('-?\+?\d*\.?\d*e?\d*', element):
                        return float(element)
                    else:
                        return element.replace(damn, '')  # there should not be a damn
                else:
                    return element
            except ValueError:
                return element

        # from https://stackoverflow.com/questions/7684333/converting-xml-to-dictionary-using-elementtree
        def etree_to_dict(t):
            d = {t.tag: {} if t.attrib else None}
            children = list(t)
            if children:
                dd = defaultdict(list)
                for dc in map(etree_to_dict, children):
                    for k, v in dc.items():
                        dd[k].append(v)
                d = {t.tag: {k: v[0] if len(v) == 1 else v
                             for k, v in dd.items()}}
            if t.attrib:
                d[t.tag].update(('@' + k, v)
                                for k, v in t.attrib.items())
            if t.text:
                text = t.text.strip()
                if children or t.attrib:
                    if text:
                        d[t.tag]['#text'] = text
                else:
                    d[t.tag] = text
            return d

        sch_URL = {'uniprot': 'https://www.uniprot.org/docs/uniprot.xsd',
                   'pfam': 'https://pfam.xfam.org/static/documents/schemas/protein.xsd'}[mode]

        etx = ET.XML(xml)
        if mode == 'uniprot':
            schema = xmlschema.XMLSchema(sch_URL)
            entry_dict = schema.to_dict(etx)['{http://uniprot.org/uniprot}entry'][0]
        elif mode == 'pfam':
            protoentry_dict = etree_to_dict(etx)
            if '{https://pfam.xfam.org/}matches' in protoentry_dict['{https://pfam.xfam.org/}pfam'][
                '{https://pfam.xfam.org/}entry']:
                entry_dict = protoentry_dict['{https://pfam.xfam.org/}pfam']['{https://pfam.xfam.org/}entry'][
                    '{https://pfam.xfam.org/}matches']['{https://pfam.xfam.org/}match']
            else:
                entry_dict = []
        else:
            raise ValueError
        damn = {'uniprot': '{http://uniprot.org/uniprot}', 'pfam': '{https://pfam.xfam.org/}'}[mode]
        return deep_clean(entry_dict, damn)

    def xml_fetch_n_parser(self, mode):  # mode = Unprot or pfam
        assert mode in ['uniprot', 'pfam'], 'Unknown mode: ' + mode
        xml = self.xml_fetcher(mode)
        return self.xml_parser(mode, xml)


    def parse_uniprot(self):
        """
        This is an rewritten version that does not use ET -> dict.
        """
        xml = self.xml_fetcher('uniprot')
        self.xml = ET.fromstring(xml)[0]
        self._parse_uniprot_xml(self.xml)
        return self

    @_failsafe
    def parse_pfam(self):
        # https://pfam.xfam.org/help#tabview=tab11
        xml = self.xml_fetch_n_parser('pfam')
        if isinstance(xml, list):
            self.pfam = xml
        elif isinstance(xml, dict):
            self.pfam = [xml]
        else:
            raise ValueError('not list or dict pfam data: ' + str(xml))

    @_failsafe
    def fetch_binders(self):
        file = os.path.join(self.settings.binders_folder, self.uniprot + '.json')
        if os.path.isfile(file):
            with open(file) as f:
                self.partners = json.load(f)
        else:
            for row in self.settings.open('ssl'):
                # CHEK1	MTOR	ENSG00000149554	ENSG00000198793	H. sapiens	BioGRID	2208318, 2342099, 2342170	3
                if self.gene_name in row:
                    protein_set = set(row.split('\t')[:2])
                    protein_set.discard(self.gene_name)
                    if protein_set == 1:
                        self.partners['SSL'].append(protein_set.pop())
                    else:  # most likely a partial match.
                        pass
                        # warn('Impossible SSL '+row)
            for row in self.settings.open('huri'):
                # Unique identifier for interactor A	Unique identifier for interactor B	Alternative identifier for interactor A	Alternative identifier for interactor B	Aliases for A	Aliases for B	Interaction detection methods	First author	Identifier of the publication	NCBI Taxonomy identifier for interactor A	NCBI Taxonomy identifier for interactor B	Interaction types	Source databases	Interaction identifier(s)	Confidence score	Complex expansion	Biological role A	Biological role B	Experimental role A	Experimental role B	Interactor type A	Interactor type B	Xref for interactor A	Xref for interactor B	Xref for the interaction	Annotations for interactor A	Annotations for interactor B	Annotations for the interaction	NCBI Taxonomy identifier for the host organism	Parameters of the interaction	Creation date	Update date	Checksum for interactor A	Checksum for interactor B	Checksum for interaction	negative	Feature(s) for interactor A	Feature(s) for interactor B	Stoichiometry for interactor A	Stoichiometry for interactor B	Participant identification method for interactor A	Participant identification method for interactor B
                # -	uniprotkb:Q6P1W5-2	ensembl:ENST00000507897.5|ensembl:ENSP00000426769.1|ensembl:ENSG00000213204.8	ensembl:ENST00000373374.7|ensembl:ENSP00000362472.3|ensembl:ENSG00000142698.14	human orfeome collection:2362(author assigned name)	human orfeome collection:5315(author assigned name)	"psi-mi:""MI:1112""(two hybrid prey pooling approach)"	Yu et al.(2011)	pubmed:21516116	taxid:9606(Homo Sapiens)	taxid:9606(Homo Sapiens)	"psi-mi:""MI:0407""(direct interaction)"	-	-	-	-	-	-	"psi-mi:""MI:0496""(bait)"	"psi-mi:""MI:0498""(prey)"	"psi-mi:""MI:0326""(protein)"	"psi-mi:""MI:0326""(protein)"	-	-	-	"comment:""vector name: pDEST-DB""|comment:""centromeric vector""|comment:""yeast strain: Y8930"""	"comment:""vector name: pDEST-AD""|comment:""centromeric vector""|comment:""yeast strain: Y8800"""	"comment:""Found in screens 1."""	taxid:4932(Saccharomyces cerevisiae)	-	6/30/2017	-	-	-	-	-	DB domain (n-terminal): gal4 dna binding domain:n-n	AD domain (n-terminal): gal4 activation domain:n-n	-	-	"psi-mi:""MI1180""(partial DNA sequence identification)"	"psi-mi:""MI1180""(partial DNA sequence identification)"
                if ':' + self.gene_name + '(' in row:  #
                    protein_set = re.findall('\:(\w+)\(gene name\)', row)
                    if len(protein_set) == 2:
                        protein_set = set(protein_set)
                        protein_set.discard(self.gene_name)
                        if protein_set == 1:
                            self.partners['HuRI'].append(protein_set.pop())
                    # t = set(row.split('\t')[:2])
                    # print(t)
                    # if t:
                    #     for p in t:
                    #         if 'uniprot' in p:
                    #             # match uniprot id to gene with go_human or similar.
                    #             match=[r for r in self.settings.open('go_human') if p.replace('uniprotkb:','') in r]
                    #             if match:
                    #                 self.partners['HuRI'].append(match[0].split('\t')[2]) #uniprotKB	A0A024RBG1	NUDT4B		GO:0003723
                    #             else:
                    #                 warn('Unmatched',p)
                    #         else:
                    #             warn('Unmatched', p)
                    else:
                        warn('Impossible HuRI ' + row)
            for row in self.settings.open('biogrid'):
                # ID Interactor A	ID Interactor B	Alt IDs Interactor A	Alt IDs Interactor B	Aliases Interactor A	Aliases Interactor B	Interaction Detection Method	Publication 1st Author	Publication Identifiers	Taxid Interactor A	Taxid Interactor B	Interaction Types	Source Database	Interaction Identifiers	Confidence Values
                # entrez gene/locuslink:6416	entrez gene/locuslink:2318	biogrid:112315|entrez gene/locuslink:MAP2K4	biogrid:108607|entrez gene/locuslink:FLNC	entrez gene/locuslink:JNKK(gene name synonym)|entrez gene/locuslink:JNKK1(gene name synonym)|entrez gene/locuslink:MAPKK4(gene name synonym)|entrez gene/locuslink:MEK4(gene name synonym)|entrez gene/locuslink:MKK4(gene name synonym)|entrez gene/locuslink:PRKMK4(gene name synonym)|entrez gene/locuslink:SAPKK-1(gene name synonym)|entrez gene/locuslink:SAPKK1(gene name synonym)|entrez gene/locuslink:SEK1(gene name synonym)|entrez gene/locuslink:SERK1(gene name synonym)|entrez gene/locuslink:SKK1(gene name synonym)	entrez gene/locuslink:ABP-280(gene name synonym)|entrez gene/locuslink:ABP280A(gene name synonym)|entrez gene/locuslink:ABPA(gene name synonym)|entrez gene/locuslink:ABPL(gene name synonym)|entrez gene/locuslink:FLN2(gene name synonym)|entrez gene/locuslink:MFM5(gene name synonym)|entrez gene/locuslink:MPD4(gene name synonym)	psi-mi:"MI:0018"(two hybrid)	"Marti A (1997)"	pubmed:9006895	taxid:9606	taxid:9606	psi-mi:"MI:0407"(direct interaction)	psi-mi:"MI:0463"(biogrid)	biogrid:103	-
                if self.gene_name in row:
                    protein_set = set([re.search('locuslink:([\w\-\.]+)\|?', e.replace('\n', '')).group(1) for e in row.split('\t')[2:4]])
                    protein_set.discard(self.gene_name)
                    if len(protein_set) == 1:
                        matched_protein = protein_set.pop()
                        self.partners['BioGRID'].append(matched_protein)
            if len(self.ENSP) > 10:
                with self.settings.open('string') as ref:
                    for row in ref:
                        if self.ENSP in row:
                            protein_set = set(row.split())
                            protein_set.discard('9606.' + self.ENSP)
                            score = 0
                            converted_gene = ''
                            for matched_protein in protein_set:
                                if matched_protein.isdigit():
                                    score = int(matched_protein)
                                else:
                                    matched_protein = matched_protein.replace('9606.', '')
                                    with self.settings.open('ensembl') as ref:
                                        for gene in ref:
                                            if matched_protein in gene:
                                                converted_gene = gene.split('\t')[2]
                                                break
                                        else:
                                            converted_gene = matched_protein  # a lie
                            if score > 900:  # highest confidence
                                self.partners['stringDB highest'].append(converted_gene)
                            elif score > 700:  # high confidence
                                self.partners['stringDB high'].append(converted_gene)
                            elif score > 400:  # medium confidence
                                self.partners['stringDB medium'].append(converted_gene)
            with open(file, 'w') as f:
                json.dump({db: list(self.partners[db]) for db in self.partners}, f)  # makes no difference downstream
        return self

    @_failsafe
    def query_ELM(self):
        assert self.uniprot, 'No uniprot entry for ELM to parse...'
        # Whole gene.
        file = os.path.join(self.settings.ELM_folder, self.uniprot + '_ELM.tsv')
        if os.path.isfile(file):
            self.log('Reading ELM data from file')
        else:
            self._assert_fetchable('ELM')
            req = requests.get('http://elm.eu.org/start_search/{id}.tsv'.format(id=self.uniprot))
            self.log('Retrieving ELM data from web')
            if req.status_code == 429:
                warn('Too many ELM requests..')
                time.sleep(60)
                self.query_ELM()
                return
            elif req.status_code != 200:
                raise ConnectionError(
                    'Could not retrieve uniprot ID data from ELM (code {0}) for gene {1}'.format(req.status_code,
                                                                                                 self.uniprot))
            response = req.text
            open(file, 'w').write(response)
        data = list(csv.DictReader(open(file, 'r'), delimiter='\t'))
        self.ELM = [(entry['start'],entry['stop'],entry['elm_identifier']) for entry in data if
                    entry['is_filtered'] == 'FALSE' or entry['is_filtered'] == 'False']

    def query_ELM_for_mutant(self):
        raise NotImplementedError
        #todo implement better.
        # Variant
        mfile = os.path.join(self.settings.ELM_variant_folder, self.uniprot + '_' + self.file_friendly_mutation + '_ELM_variant.tsv')
        if os.path.isfile(mfile):
            self.log('Reading ELM variant data from file')
        else:
            self._assert_fetchable('ELM (variant)')
            mseq = self.seq[0:self.resi - 1] + self.to_resn + self.seq[self.resi:]
            mseq_trim = mseq[max(self.resi - 30, 0):min(self.resi + 30, len(self.seq))]
            assert mseq_trim, 'Something went wrong in parsing {0}'.format(mseq)
            req = requests.get('http://elm.eu.org/start_search/{}'.format(mseq_trim))
            if req.status_code == 429:
                time.sleep(60)
                self.query_ELM()
                return
            elif req.status_code != 200:
                raise ConnectionError('Could not retrieve data: ' + req.text)
            response = req.text
            self.log('Retrieving ELM variant data from web')
            open(mfile, 'w').write(response)
        mdata = list(csv.DictReader(open(mfile, 'r'), delimiter='\t'))
        lbound = max(self.resi - 30, 0)
        ubound = min(self.resi + 30, len(self.seq))
        for entry in data:
            if int(entry['start']) < self.resi < int(entry['stop']):
                for m in range(len(mdata)):
                    mentry = mdata[m]
                    if entry['elm_identifier'] == mentry['elm_identifier'] and int(entry['start']) == int(
                            mentry['start']) + lbound:
                        del mdata[m]
                        self.mutational_effect.append(
                            'Possible linear motif <a target="_blank" href="http://elm.eu.org/elms/{0}">{0} <i class="fas fa-external-link-square"></i></a> spanning mutation: preserved'.format(
                                entry['elm_identifier']))
                        break
                else:
                    self.mutational_effect.append(
                        'Possible linear motif <a target="_blank" href="http://elm.eu.org/elms/{0}">{0} <i class="fas fa-external-link-square"></i></a> spanning mutation: lost'.format(
                            entry['elm_identifier']))
        if mdata:
            for entry in mdata:
                if int(entry['start']) < 30 and int(entry['stop']) > 30:
                    self.mutational_effect.append(
                        'Possible new linear motif <a target="_blank" href="http://elm.eu.org/elms/{0}">{0} <i class="fas fa-external-link-square"></i></a> spanning mutation'.format(
                            entry['elm_identifier']))

    def fetch_ENSP(self):
        """EMBL ids should vome from the Unirpto entry. However, in some cases too many it is absent."""
        for row in self.settings.open('ensembl'):
            if self.gene_name in row:
                if row.count(' ') > 5:
                    self.ENSP = row.split()[5]
                break
        else:
            warn('Unknown Ensembl protein id for ' + self.gene_name)
        return self


    @_failsafe
    def query_ELM(self):
        assert self.uniprot, 'No uniprot entry for ELM to parse...'
        # Whole gene.
        file = os.path.join(self.settings.ELM_folder, self.uniprot + '_ELM.tsv')
        if os.path.isfile(file):
            self.log('Reading ELM data from file')
        else:
            self._assert_fetchable('ELM')
            req = requests.get('http://elm.eu.org/start_search/{id}.tsv'.format(id=self.uniprot))
            self.log('Retrieving ELM data from web')
            if req.status_code == 429:
                warn('Too many ELM requests..')
                time.sleep(60)
                self.query_ELM()
                return
            elif req.status_code != 200:
                raise ConnectionError(
                    'Could not retrieve uniprot ID data from ELM (code {0}) for gene {1}'.format(req.status_code,
                                                                                                 self.uniprot))
            response = req.text
            open(file, 'w').write(response)
        with open(file, 'r') as fh:
            data = list(csv.DictReader(fh, delimiter='\t'))
            self.ELM = [(entry['start'],entry['stop'],entry['elm_identifier']) for entry in data if
                        entry['is_filtered'] == 'FALSE' or entry['is_filtered'] == 'False']
        return self

    @_failsafe
    def parse_pLI(self):
        for line in csv.DictReader(self.settings.open('ExAC_pLI'), delimiter='\t'):
            # transcript	gene	chr	n_exons	cds_start	cds_end	bp	mu_syn	mu_mis	mu_lof	n_syn	n_mis	n_lof	exp_syn	exp_mis	exp_lof	syn_z	mis_z	lof_z	pLI	pRec	pNull
            if self.gene_name == line['gene']:
                self.pLI = float(line['pLI'])  # intolerant of a single loss-of-function variant (like haploinsufficient genes, observed ~ 0.1*expected)
                self.pRec = float(line['pRec'])  # intolerant of two loss-of-function variants (like recessive genes, observed ~ 0.5*expected)
                self.pNull = float(line['pNull'])  # completely tolerant of loss-of-function variation (observed = expected)
                break
        else:
            warn('Gene {} not found in ExAC table.'.format(self.gene_name))
        return self

    @_failsafe
    def get_ExAC(self):

        def parse(line):
            data = dict(zip(['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'], line.split('\t')))
            # unpack info
            info = {x.split('=')[0]: x.split('=')[1] for x in data['INFO'].split(';') if '=' in x}
            del data['INFO']
            # definitions were taken from the .vep.vcf
            csq = dict(zip(
                'Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|GMAF|AFR_MAF|AMR_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|ExAC_MAF|ExAC_Adj_MAF|ExAC_AFR_MAF|ExAC_AMR_MAF|ExAC_EAS_MAF|ExAC_FIN_MAF|ExAC_NFE_MAF|ExAC_OTH_MAF|ExAC_SAS_MAF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info|context|ancestral'.split(
                    '|'), info['CSQ'].split('|')))
            del info['CSQ']
            return {**data, **info, **csq}

        if self.pLI == -1:  # No ExAC name so don't bother.
            return
        file = os.path.join(self.settings.ExAC_folder, self.uniprot + '_ExAC.vcf')
        matches = []
        if os.path.isfile(file):
            self.log('Reading from cached _ExAC.vcf file')
            matches = list(open(file))
        else:
            self.log('Parsing ExAC.r1.sites.vep.vcf')
            # out.write('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO')
            gmatch = '|{}|'.format(self.gene_name)
            with self.settings.open('ExAC_vep') as fh:
                for line in fh:
                    if line and line[0] != '#' and gmatch in line:
                        matches.append(line)
                        # data = parse(line)
            with open(file, 'w') as out:
                out.writelines(matches)
            self.log('... {} Matches found'.format(len(matches)))
        # parse entries. Oddly some are incorrect if the gene happens to be potentially ambiguous.
        parsed = [parse(line) for line in matches]
        self.alleles = [x for x in parsed if x['SYMBOL'] == self.gene_name]
        if self.alleles:
            self.ENSG = self.alleles[0]['Gene']
        verdict_list = [self.verify_allele(a) for a in self.iter_allele()]
        verdict_sum = sum(verdict_list)
        verdict_len = len(verdict_list)
        if verdict_len and verdict_sum / verdict_len > 0.8:
            self.log('ExAC data usable. There are {s} valid out of {l}'.format(s=verdict_sum, l=verdict_len))
        else:
            self.log('ExAC data cannot be integrated. There are {s} valid out of {l}'.format(s=verdict_sum, l=verdict_len))
            self.alleles = []  # no ExAC.

    def iter_allele(self, filter=True, consequence=None, split=True):
        """
        Generator that gives the allele...
        :param filter: Bool. True (default) means only PASS records, False all.
        :param consequence: Str. Optional filter. 'intron_variant', 'missense_variant', '5_prime_UTR_variant', 'synonymous_variant', 'splice_region_variant&intron_variant', '3_prime_UTR_variant', 'missense_variant&splice_region_variant', 'stop_gained', 'splice_region_variant&synonymous_variant', 'inframe_deletion', 'frameshift_variant', 'splice_region_variant&5_prime_UTR_variant', 'splice_acceptor_variant', 'splice_donor_variant'
        :return: yields (resi as str but maybe range, X/Y)
        """
        for allele in self.alleles:
            if allele['HGVSp']:  # ['HGVSp']
                if filter:
                    if allele['FILTER'] != 'PASS':
                        continue
                if consequence:
                    if allele['Consequence'] != consequence and allele['Consequence'] not in consequence:
                        continue  # 'Protein_position': '', 'Amino_acids'
                elif allele['Consequence'] in ('synonymous_variant', 'splice_region_variant&synonymous_variant'):
                    continue
                else:
                    pass  # It is fine.
                if split:
                    yield (allele['Protein_position'], allele['Amino_acids'])
                else:
                    yield allele['HGVSp'].split(':')[-1]

    def verify_allele(self, allele):
        s = allele[0].split('-')[0]
        if not s.isdigit():
            return False
        i = int(s)
        # There is the problem that some ExAC entries do not match.
        if i > len(self.seq):
            return False
        elif self.seq[i - 1] == allele[1].split('/')[0]:
            return True
        else:
            return False

    # @_failsafe
    def add_manual_data(self):
        man = {x['Gene']: x for x in csv.DictReader(open(os.path.join(self.settings.manual_folder, 'manual.csv')))}
        ### PDB
        pdb_candidate = os.path.join(self.settings.manual_folder, self.uniprot + '.pdb')
        if os.path.isfile(pdb_candidate):
            pdb_man_file = pdb_candidate
        elif self.uniprot in man and man[self.uniprot]['PDB']:
            pdb_candidate = os.path.join(self.settings.manual_folder, man[self.uniprot]['PDB'].lstrip().rstrip())
            if os.path.isfile(pdb_candidate):
                pdb_man_file = pdb_candidate
            else:  # csv is wrong
                pdb_man_file = None
                warn('Cannot find structure file, claimed by manual.csv: '+pdb_candidate)
        else:
            pdb_man_file = None
        if pdb_man_file:
            self.pdb_file = os.path.join(self.settings.page_folder, self.uniprot + '.pdb')
            self.pdb_resi = 0  # self.resi - int(man[self.gene]['PDB_offset'])
            if not os.path.isfile(self.pdb_file):  # make the html copy
                copyfile(pdb_man_file, self.pdb_file)
        ### MD
        txt_candidate = os.path.join(self.settings.manual_folder, self.uniprot + '.md')
        if os.path.isfile(txt_candidate):
            txt = txt_candidate
        elif self.gene_name in man and man[self.uniprot]['Text']:
            txt_candidate = os.path.join(self.settings.manual_folder, man[self.uniprot]['Text'])
            if os.path.isfile(txt_candidate):
                txt = txt_candidate
            else:
                warn('Cannot find text file, ' + txt_candidate + ' claimed by manual.csv.')
                txt = None
        else:
            txt = None
        if txt:
            self.user_text = markdown.markdown(open(txt).read())
        # log
        if self.user_text and self.pdb_file:
            self.log('Manual data added.')
        elif self.user_text and not self.pdb_file:
            self.log('Manual data. No structure.')
        else:  # this is unnedded except for legacy pickled files.
            self.user_text = ''
            self.pdb_file = ''
        # done
        return self

    @_failsafe
    def match_pdb(self):
        # variable self.matched_pdbs not used...
        file = os.path.join(self.settings.pdb_blast_folder, self.uniprot + '_blastPDB.xml')
        if os.path.isfile(file):
            pass
        else:
            raise Exception('This is impossible. Preparsed.') #todo make a settings flag for this
            self._assert_fetchable(self.gene_name + ' PDB match')
            self.log('Blasting against PDB dataset of NCBI')
            result_handle = NCBIWWW.qblast("blastp", "pdb", self.seq)  # the whole sequence.
            open(file, "w").write(result_handle.read())
        blast_record = NCBIXML.read(open(file))
        # Parse!
        matches = []
        for align in blast_record.alignments:
            for hsp in align.hsps:
                if hsp.score > 100:
                    d = {'match': align.title[0:50], 'match_score': hsp.score, 'match_start': hsp.query_start, 'match_length': hsp.align_length, 'match_identity': hsp.identities / hsp.align_length}
                    self.pdb_matches.append(d)
                    #if hsp.query_start < self.resi < hsp.align_length + hsp.query_start:

                    # self.pdb_resi=self.resi+(hsp.sbjct_start-hsp.query_start) #this is a massive pickle as the NCBI pdb data has different numbering.
        if matches:
            self.log('GENE MUTANT WITHIN CRYSTAL STRUCTURE!!')
            # mostly based on identity, but if a long one exists and it is only marginally worse the go for that.
            best = sorted(matches, key=lambda m: m['match_length'] * m['match_identity'] ** 3, reverse=True)[0]
            for k in best:
                setattr(self, k, best[k]) #forgot what the hell this is!
        else:
            self.log('UNCRYSTALLISED MUTATION')
        return self

    def parse_all(self, mode='parallel'):
        """
        Gets all the data. if running in parallel see self._threads list.
        :param mode: parallel | background (=parallel but not complete) | serial (or anything else)
        :return:
        """
        ### Sanity check
        assert self.gene_name, 'There is no gene name'
        if not self.uniprot:
            warn('There is no uniprot value in this entry!')
            self.uniprot = json.load(open(os.path.join(self.settings.data_folder,'human_prot_namedex.json')))[self.gene_name]
        ### fetch!
        tasks = {'Uniprot': self.parse_uniprot, #done.
                 'PFam': self.parse_pfam, #done.
                 'ELM': self.query_ELM,
                 'pLI': self.parse_pLI,
                 'ExAC': self.get_ExAC,
                 'manual': self.add_manual_data,
                 #'GO terms': self.fetch_go,
                 'Binding partners': self.fetch_binders}
        #tasks = {'Uniprot': self.parse_uniprot}
        if mode in ('parallel','background'):
            threads = {}
            for k, fn in tasks.items():
                t = threading.Thread(target=fn)
                t.start()
                threads[k] = t
            if mode == 'parallel':
                for tn, t in threads.items():
                    t.join()
                return self
            else:
                self._threads = threads
                return self
        else:
            for task_fn in tasks.values():
                task_fn()
        return self

    def predict_effect(self):
        pass
        ##TOdo write>!

    #to do fix
    def parse_ExAC_type(self):
        if self.pLI < 0:  # error.
            self.ExAC_type='Unknown'
        elif self.pLI > max(self.pRec, self.pNull):
            self.ExAC_type ='Dominant'
        elif self.pRec > max(self.pLI, self.pNull):
            self.ExAC_type ='Recessive'
        elif self.pNull > max(self.pLI, self.pRec):
            self.ExAC_type ='None'
        else:
            self.ExAC_type ='Unknown'
        return self

    @_failsafe
    def parse_pLI(self):
        for line in csv.DictReader(self.settings.open('ExAC_pLI'), delimiter='\t'):
            # transcript	gene	chr	n_exons	cds_start	cds_end	bp	mu_syn	mu_mis	mu_lof	n_syn	n_mis	n_lof	exp_syn	exp_mis	exp_lof	syn_z	mis_z	lof_z	pLI	pRec	pNull
            if self.gene_name == line['gene']:
                self.pLI = float(line['pLI'])  # intolerant of a single loss-of-function variant (like haploinsufficient genes, observed ~ 0.1*expected)
                self.pRec = float(line['pRec'])  # intolerant of two loss-of-function variants (like recessive genes, observed ~ 0.5*expected)
                self.pNull = float(line['pNull'])  # completely tolerant of loss-of-function variation (observed = expected)
                self.parse_ExAC_type()
                break
        else:
            warn('Gene {} not found in ExAC table.'.format(self.gene))
        return self


    ## depraction zone.
    def write(self, file=None):
        raise Exception('DEPRACATED. use write_uniprot')

    @_failsafe
    def _test_failsafe(self):
        raise ValueError('Will failsafe catch it? ({})'.format(self.settings.error_tollerant))

    def check_mutation(self, mutation):
        if len(self.sequence) > mutation.residue_index and self.sequence[mutation.residue_index - 1] == mutation.from_residue:
            return True
        else:
            return False # call mutation_discrepancy to see why.

    def mutation_discrepancy(self, mutation):
        # returns a string explaining the check_mutation discrepancy
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


class Mutation:
    """
    Moved out of Protein.
    """
    def __init__(self, mutation=None):
        self.from_residue = ''
        self.to_residue = ''
        self.residue_index = 0
        if mutation:
            self.parse_mutation(mutation)

    def __str__(self):
        # note taht this is not file-safe
        return self.from_residue+str(self.residue_index)+self.to_residue

    #raise NotImplementedError('Under upgrade')
    def parse_mutation(self, mutation):
        ##### clean
        assert mutation.find('.c') == -1, 'Chromosome mutation not accepted. Use Protein.'
        # remove the p.
        mutation = mutation.replace('p.', '').replace('P.', '')
        for (single, triple) in (
                ('A', 'Ala'), ('B', 'Asx'), ('C', 'Cys'), ('D', 'Asp'), ('E', 'Glu'), ('F', 'Phe'), ('G', 'Gly'),
                ('H', 'His'),
                ('I', 'Ile'), ('K', 'Lys'), ('L', 'Leu'), ('M', 'Met'), ('N', 'Asn'), ('P', 'Pro'), ('Q', 'Gln'),
                ('R', 'Arg'),
                ('S', 'Ser'), ('T', 'Thr'), ('V', 'Val'), ('W', 'Trp'), ('X', 'Xaa'), ('Y', 'Tyr'), ('Z', 'Glx')):
            if mutation.find(triple) != -1 or mutation.find(triple.lower()) != -1 or mutation.find(triple.upper()):
                mutation = mutation.replace(triple, single).replace(triple.upper(), single).replace(triple.lower(), single)
        self.mutation = mutation
        ###### split
        if self.mutation.find('fs') != -1 or self.mutation.find('*') != -1:
            rex = re.match('(\w)(\d+)', self.mutation)
            if rex:
                self.from_residue = rex.group(1)
                self.residue_index = int(rex.group(2))
                self.to_residue = '*'
                self.file_friendly_mutation = self.from_residue + str(self.residue_index) + 'stop'
                self.clean_mutation = self.from_residue + str(self.residue_index) + '*'
            else:
                raise ValueError('odd mutation of type fs' + self.mutation)
        elif self.mutation.find('del') != -1 or self.mutation.find('\N{GREEK CAPITAL LETTER DELTA}') != -1:
            self.mutation = self.mutation.replace('\N{GREEK CAPITAL LETTER DELTA}',
                                                  'del')  # would love the other way round but unicode in 2018 still crashes stuff...
            rex = re.match('del(\w)(\d+)', self.mutation)
            if not rex:
                rex = re.match('(\w)(\d+)del', self.mutation)
            if rex:
                self.from_residue = rex.group(1)
                self.residue_index = int(rex.group(2))
                self.to_residue = rex.group(1)  # technically not but shush...
                self.file_friendly_mutation = 'del' + self.from_residue + str(self.residue_index)
                self.clean_mutation = 'del' + self.from_residue + str(self.residue_index)
                warn('Mutation parsing: Deletions are not handled correctly atm...')
                self.log('Residue deletion is not handled correctly at the moment...')
            else:
                raise ValueError('odd mutation of type deletion' + self.mutation)
        else:
            rex = re.match('(\w)(\d+)(\w)', self.mutation)
            if rex:
                assert rex.group(1) in 'QWERTYIPASDFGHKLCVNM*', 'The from mutant is not a valid amino acid'
                assert rex.group(3) in 'QWERTYIPASDFGHKLCVNM*', 'The to mutant is not a valid amino acid'
                self.from_residue = rex.group(1)
                self.residue_index = int(rex.group(2))
                self.to_residue = rex.group(3)
                self.file_friendly_mutation = self.from_residue + str(self.residue_index) + self.to_residue
                self.clean_mutation = self.from_residue + str(self.residue_index) + self.to_residue
                #self.blossum_score = self.blossum[self.from_residue][self.to_residue]
            else:
                raise ValueError(self.mutation + ' is an odd_mutation')
        return self
