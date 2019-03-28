# Protein inherits _UniprotMixin, which in turn inherits _BaseMixin
# `.settings` class attribute is global_settings from settings_handler.py and is added by _BaseMixin.
# _BaseMixin is inherited by _UniprotMixin contains _failsafe decorator, __getattr__ and settings

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
from warnings import warn
from shutil import copyfile
import requests  # for xml fetcher.

from .ET_monkeypatch import ET  # monkeypatched version

from ._protein_uniprot_mixin import _UniprotMixin
from ._protein_base_mixin import _BaseMixin
from ._protein_disused_mixin import _DisusedMixin
from ._protein_lite import ProteinLite

class Protein(ProteinLite, _BaseMixin, _DisusedMixin, _UniprotMixin):
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
            if self.settings.error_tolerant:
                try:
                    return func(self, *args, **kargs)
                except Exception as error:
                    warn('Error caught in method `Protein().{n}`: {e}'.format(n=func.__name__, e=error))
                    return None
            else:
                return func(self, *args, **kargs)

        return wrapper

    ############################# INIT #############################
    @classmethod
    def from_uniprot(cls, entry):
        self=cls()
        self._parse_uniprot_xml(entry)

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
    def fetch_binders(self):
        file = os.path.join(self.settings.binders_folder, self.uniprot + '.json')
        if os.path.isfile(file):
            with open(file) as f:
                self.partners = json.load(f)
        else:
            ###################### SSL http://slorth.biochem.sussex.ac.uk/download/h.sapiens_ssl_predictions.csv
            self.log('parsing SSL...')
            for row in self.settings.open('ssl'):
                # CHEK1	MTOR	ENSG00000149554	ENSG00000198793	H. sapiens	BioGRID	2208318, 2342099, 2342170	3
                if self.gene_name in row:
                    protein_set = set(row.split('\t')[:2])
                    protein_set.discard(self.gene_name)
                    if len(protein_set) == 1:
                        self.partners['SSL'].append(protein_set.pop())
                    else:  # most likely a partial match. maybe???
                        warn('Impossible SSL '+row)
                        warn('{0} after discarding {1} erroneously became {2}'.format(protein_set, self.gene_name, set(row.split('\t')[:2])))
            ###################### HURI cat.psi
            self.log('parsing HURI...')
            for row in self.settings.open('huri'):
                # Unique identifier for interactor A	Unique identifier for interactor B	Alternative identifier for interactor A	Alternative identifier for interactor B	Aliases for A	Aliases for B	Interaction detection methods	First author	Identifier of the publication	NCBI Taxonomy identifier for interactor A	NCBI Taxonomy identifier for interactor B	Interaction types	Source databases	Interaction identifier(s)	Confidence score	Complex expansion	Biological role A	Biological role B	Experimental role A	Experimental role B	Interactor type A	Interactor type B	Xref for interactor A	Xref for interactor B	Xref for the interaction	Annotations for interactor A	Annotations for interactor B	Annotations for the interaction	NCBI Taxonomy identifier for the host organism	Parameters of the interaction	Creation date	Update date	Checksum for interactor A	Checksum for interactor B	Checksum for interaction	negative	Feature(s) for interactor A	Feature(s) for interactor B	Stoichiometry for interactor A	Stoichiometry for interactor B	Participant identification method for interactor A	Participant identification method for interactor B
                # -	uniprotkb:Q6P1W5-2	ensembl:ENST00000507897.5|ensembl:ENSP00000426769.1|ensembl:ENSG00000213204.8	ensembl:ENST00000373374.7|ensembl:ENSP00000362472.3|ensembl:ENSG00000142698.14	human orfeome collection:2362(author assigned name)	human orfeome collection:5315(author assigned name)	"psi-mi:""MI:1112""(two hybrid prey pooling approach)"	Yu et al.(2011)	pubmed:21516116	taxid:9606(Homo Sapiens)	taxid:9606(Homo Sapiens)	"psi-mi:""MI:0407""(direct interaction)"	-	-	-	-	-	-	"psi-mi:""MI:0496""(bait)"	"psi-mi:""MI:0498""(prey)"	"psi-mi:""MI:0326""(protein)"	"psi-mi:""MI:0326""(protein)"	-	-	-	"comment:""vector name: pDEST-DB""|comment:""centromeric vector""|comment:""yeast strain: Y8930"""	"comment:""vector name: pDEST-AD""|comment:""centromeric vector""|comment:""yeast strain: Y8800"""	"comment:""Found in screens 1."""	taxid:4932(Saccharomyces cerevisiae)	-	6/30/2017	-	-	-	-	-	DB domain (n-terminal): gal4 dna binding domain:n-n	AD domain (n-terminal): gal4 activation domain:n-n	-	-	"psi-mi:""MI1180""(partial DNA sequence identification)"	"psi-mi:""MI1180""(partial DNA sequence identification)"
                if ':' + self.gene_name + '(' in row:  #
                    protein_set = re.findall('\:(\w+)\(gene name\)', row)
                    if len(protein_set) == 2:
                        protein_set = set(protein_set)
                        protein_set.discard(self.gene_name)
                        if len(protein_set) == 1:
                            self.partners['HuRI'].append(protein_set.pop())
                        elif len(protein_set) == 0: ### multimeric
                            pass
                        else:
                            warn('Impossible HURI ' + row)
                            warn('{0} after discarding {1} erroneously became {2}'.format(protein_set, self.gene_name, re.findall('\:(\w+)\(gene name\)', row)))

                    else:
                        warn('Impossible HURI ' + row)
            ###################### biogrid https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-3.5.166/BIOGRID-ALL-3.5.166.mitab.zip
            self.log('parsing biogrid...')
            for row in self.settings.open('biogrid'):
                # ID Interactor A	ID Interactor B	Alt IDs Interactor A	Alt IDs Interactor B	Aliases Interactor A	Aliases Interactor B	Interaction Detection Method	Publication 1st Author	Publication Identifiers	Taxid Interactor A	Taxid Interactor B	Interaction Types	Source Database	Interaction Identifiers	Confidence Values
                # entrez gene/locuslink:6416	entrez gene/locuslink:2318	biogrid:112315|entrez gene/locuslink:MAP2K4	biogrid:108607|entrez gene/locuslink:FLNC	entrez gene/locuslink:JNKK(gene name synonym)|entrez gene/locuslink:JNKK1(gene name synonym)|entrez gene/locuslink:MAPKK4(gene name synonym)|entrez gene/locuslink:MEK4(gene name synonym)|entrez gene/locuslink:MKK4(gene name synonym)|entrez gene/locuslink:PRKMK4(gene name synonym)|entrez gene/locuslink:SAPKK-1(gene name synonym)|entrez gene/locuslink:SAPKK1(gene name synonym)|entrez gene/locuslink:SEK1(gene name synonym)|entrez gene/locuslink:SERK1(gene name synonym)|entrez gene/locuslink:SKK1(gene name synonym)	entrez gene/locuslink:ABP-280(gene name synonym)|entrez gene/locuslink:ABP280A(gene name synonym)|entrez gene/locuslink:ABPA(gene name synonym)|entrez gene/locuslink:ABPL(gene name synonym)|entrez gene/locuslink:FLN2(gene name synonym)|entrez gene/locuslink:MFM5(gene name synonym)|entrez gene/locuslink:MPD4(gene name synonym)	psi-mi:"MI:0018"(two hybrid)	"Marti A (1997)"	pubmed:9006895	taxid:9606	taxid:9606	psi-mi:"MI:0407"(direct interaction)	psi-mi:"MI:0463"(biogrid)	biogrid:103	-
                if self.gene_name in row:
                    protein_set = set()
                    for e in row.split('\t')[2:4]:
                        rex = re.search('locuslink:([\w\-\.]+)\|?', e.replace('\n', ''))
                        if rex:
                            protein_set.add(rex.group(1))
                        else:
                            pass
                            #warn('locuslink not found in {}'.format(e.replace('\n', '')))
                    protein_set.discard(self.gene_name)
                    if len(protein_set) == 1:
                        matched_protein = protein_set.pop()
                        self.partners['BioGRID'].append(matched_protein)
            if len(self.ENSP) > 10:
                ###################### biogrid https://stringdb-static.org/download/protein.links.v10.5/9606.protein.links.v10.5.txt.gz
                self.log('parsing string...')
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
                            else:
                                self.partners['stringDB low'].append(converted_gene)
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
                self.log('pLI: {pLI}, pRec {pRec}, pNull {pNull}'.format(**line))
                break
        else:
            msg = 'Gene {} not found in ExAC table.'.format(self.gene_name)
            warn(msg)
            self.log(msg)
        return self

    @_failsafe
    def parse_gNOMAD(self):
        # preparsed.
        # from protein.generate.split_gNOMAD import gNOMAD
        # gNOMAD().write('gNOMAD')
        file = os.path.join(self.settings.temp_folder, 'gNOMAD', self.uniprot +'.json')
        if os.path.exists(file):
            for snp in json.load(open(file)):
                resi = int(snp['residue_index'].split('-')[0])
                self.gNOMAD.append({'id': 'gNOMAD_{x}_{x}_{id}'.format(x=resi, id=snp['id']),
                                    'description': '{from_residue}{residue_index}{to_residue} ({id})'.format(**snp),
                                    'x': resi,
                                    'y': resi,
                                    'impact': snp['impact']
                                    })
        else:
            self.gNOMAD = []
            warn('No gNOMAD data from {0} {1}?'.format(self.gene_name, self.uniprot))
        self.log('gNOMAD mutations: {0}'.format(len(self.gNOMAD)))
        return self






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

    def parse_all(self, mode='parallel'):
        """
        Gets all the data. if running in parallel see self._threads list.
        :param mode: parallel | background (=parallel but not complete) | serial (or anything else)
        :return:
        """
        ### Sanity check
        if not self.uniprot:
            warn('There is no uniprot value in this entry!')
            self.uniprot = json.load(open(os.path.join(self.settings.data_folder,'human_prot_namedex.json')))[self.gene_name]
        if not self.gene_name:
            self.parse_uniprot()
        ### fetch!
        tasks = {'Uniprot': self.parse_uniprot, #done.
                 'PFam': self.parse_pfam, #done.
                 'Swissmodel': self.parse_swissmodel,
                 'pLI': self.parse_pLI,
                 'gNOMAD': self.parse_gNOMAD,
                 'parse_pdb_blast': self.parse_pdb_blast,
                 'manual': self.add_manual_data,
                 'Binding partners': self.fetch_binders}
        #prot.get_percent_modelled()
        #prot.parse_ExAC_type()
        #prot.dump()
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
            warn('Gene {} not found in ExAC table.'.format(self.gene_name))
        return self


    def parse_swissmodel(self):
        models=json.load(self.settings.open('swissmodel'))['index']
        for model in models:
            if self.uniprot == model['uniprot_ac']:
                self.swissmodel.append({'x': model['from'],
                                        'y': model['to'],
                                        'id': model['coordinate_id'],
                                        'description': '{template} (id:{seqid:.0}%)'.format(**model),
                                        'url': model['url']})
        self.log('Swissmodel has {0} models.'.format(len(self.swissmodel)))
        return self


    def get_percent_modelled(self):
        """Returns the [0,1] fraction of the protein length that is modelled (repeats arent recounted)"""

        def clean(text):
            return int(text.lstrip().replace('-', ' ').replace(',', ' ').split(' ')[0]) - 1

        try:
            if len(self) == 0:
                return 0
            state = [False for i in range(len(self))]
            for dataset in (self.swissmodel,self.pdbs):
                for model in dataset:
                    for i in range(clean(model['x']),clean(model['y'])):
                        state[i]=True
            self.percent_modelled = sum(state)/len(self)
        except:
            self.percent_modelled = 0
        return self.percent_modelled

    def parse_pdb_blast(self):
        file = os.path.join(self.settings.temp_folder, 'blastpdb2', self.uniprot + '.json')
        if os.path.exists(file):
            self.pdb_matches = json.load(open(file))
        else:
            warn('No PDB blast data from {0} {1}?'.format(self.gene_name, self.uniprot))
        return self



    @_failsafe
    def _test_failsafe(self):
        raise ValueError('Will failsafe catch it? ({})'.format(self.settings.error_tollerant))
