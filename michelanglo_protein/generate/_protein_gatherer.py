__doc__ = """
Protein inherits _UniprotMixin, which in turn inherits _BaseMixin
`.settings` class attribute is global_settings from settings_handler.py and is added by _BaseMixin.
_BaseMixin is inherited by _UniprotMixin contains _failsafe decorator, __getattr__ and settings
"""

import os
import pickle
import re
import csv
import json
import time
import markdown
# import xmlschema
import threading
from collections import defaultdict
from warnings import warn
from shutil import copyfile
import requests  # for xml fetcher.

from .ET_monkeypatch import ET  # monkeypatched version

from ._protein_uniprot_mixin import _UniprotMixin
from ._protein_base_mixin import _BaseMixin
from ._protein_disused_mixin import _DisusedMixin
from ..core import ProteinCore, Variant, Structure
from Bio.SeqUtils import ProtParam, ProtParamData


class ProteinGatherer(ProteinCore, _BaseMixin, _DisusedMixin, _UniprotMixin):
    """
    This class handles each protein entry from Uniprot. See __init__ for the list of variables stored.
    It fills them from various other sources.

        >>> ProteinGatherer()
        >>> ProteinGatherer.from_uniprot(xml_entry)
        >>> ProteinGatherer.load('filename')

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

    ############### INIT ###############
    @classmethod
    def from_uniprot(cls, entry):
        '''This requires a xml entry'''
        self = cls()
        self.xml = entry
        self._parse_uniprot_xml(entry)
        return self

    def write_uniprot(self, file=None):
        if not file:
            file = self.uniprot_name + '.xml'
        with open(file, 'w') as w:
            w.write(
                '<?xml version="1.0" encoding="UTF-8"?><uniprot xmlns="http://uniprot.org/uniprot" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" ' +
                'xsi:schemaLocation="http://uniprot.org/uniprot http://www.uniprot.org/docs/uniprot.xsd">')
            if isinstance(self.xml, str):
                w.write(self.xml)
            else:
                w.write(ET.tostring(self.xml).decode())
            w.write('</uniprot>')

    ############### Data gathering ###############
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
                requestURL = 'https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=100&accession={acc}'.format(
                    # &taxid=9606
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
            ########### SSL http://slorth.biochem.sussex.ac.uk/download/h.sapiens_ssl_predictions.csv
            self.log('parsing SSL...')
            for row in self.settings.open('ssl'):
                # CHEK1	MTOR	ENSG00000149554	ENSG00000198793	H. sapiens	BioGRID	2208318, 2342099, 2342170	3
                if self.gene_name in row:
                    protein_set = set(row.split('\t')[:2])
                    protein_set.discard(self.gene_name)
                    if len(protein_set) == 1:
                        self.partners['SSL'].append(protein_set.pop())
                    else:  # most likely a partial match. maybe???
                        warn('Impossible SSL ' + row)
                        warn('{0} after discarding {1} erroneously became {2}'.format(protein_set, self.gene_name,
                                                                                      set(row.split('\t')[:2])))
            ########### HURI cat.psi
            self.log('parsing HURI...')
            for row in self.settings.open('huri'):
                # Unique identifier for interactor A	Unique identifier for interactor B	Alternative identifier for interactor A	Alternative identifier for interactor B	Aliases for A	Aliases for B	Interaction detection methods	First author	Identifier of the publication	NCBI Taxonomy identifier for interactor A	NCBI Taxonomy identifier for interactor B	Interaction types	Source databases	Interaction identifier(s)	Confidence score	Complex expansion	Biological role A	Biological role B	Experimental role A	Experimental role B	Interactor type A	Interactor type B	Xref for interactor A	Xref for interactor B	Xref for the interaction	Annotations for interactor A	Annotations for interactor B	Annotations for the interaction	NCBI Taxonomy identifier for the host organism	Parameters of the interaction	Creation date	Update date	Checksum for interactor A	Checksum for interactor B	Checksum for interaction	negative	Feature(s) for interactor A	Feature(s) for interactor B	Stoichiometry for interactor A	Stoichiometry for interactor B	Participant identification method for interactor A	Participant identification method for interactor B
                # -	uniprotkb:Q6P1W5-2	ensembl:ENST00000507897.5|ensembl:ENSP00000426769.1|ensembl:ENSG00000213204.8	ensembl:ENST00000373374.7|ensembl:ENSP00000362472.3|ensembl:ENSG00000142698.14	human orfeome collection:2362(author assigned name)	human orfeome collection:5315(author assigned name)	"psi-mi:""MI:1112""(two hybrid prey pooling approach)"	Yu et al.(2011)	pubmed:21516116	taxid:9606(Homo Sapiens)	taxid:9606(Homo Sapiens)	"psi-mi:""MI:0407""(direct interaction)"	-	-	-	-	-	-	"psi-mi:""MI:0496""(bait)"	"psi-mi:""MI:0498""(prey)"	"psi-mi:""MI:0326""(protein)"	"psi-mi:""MI:0326""(michelanglo_protein)"	-	-	-	"comment:""vector name: pDEST-DB""|comment:""centromeric vector""|comment:""yeast strain: Y8930"""	"comment:""vector name: pDEST-AD""|comment:""centromeric vector""|comment:""yeast strain: Y8800"""	"comment:""Found in screens 1."""	taxid:4932(Saccharomyces cerevisiae)	-	6/30/2017	-	-	-	-	-	DB domain (n-terminal): gal4 dna binding domain:n-n	AD domain (n-terminal): gal4 activation domain:n-n	-	-	"psi-mi:""MI1180""(partial DNA sequence identification)"	"psi-mi:""MI1180""(partial DNA sequence identification)"
                if ':' + self.gene_name + '(' in row:  #
                    protein_set = re.findall('\:(\w+)\(gene name\)', row)
                    if len(protein_set) == 2:
                        protein_set = set(protein_set)
                        protein_set.discard(self.gene_name)
                        if len(protein_set) == 1:
                            self.partners['HuRI'].append(protein_set.pop())
                        elif len(protein_set) == 0:  ## multimeric
                            pass
                        else:
                            warn('Impossible HURI ' + row)
                            warn('{0} after discarding {1} erroneously became {2}'.format(protein_set, self.gene_name,
                                                                                          re.findall(
                                                                                              '\:(\w+)\(gene name\)',
                                                                                              row)))

                    else:
                        warn('Impossible HURI ' + row)
            ########### biogrid https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-3.5.166/BIOGRID-ALL-3.5.166.mitab.zip
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
                            # warn('locuslink not found in {}'.format(e.replace('\n', '')))
                    protein_set.discard(self.gene_name)
                    if len(protein_set) == 1:
                        matched_protein = protein_set.pop()
                        self.partners['BioGRID'].append(matched_protein)
            if len(self.ENSP) > 10:
                ########### biogrid https://stringdb-static.org/download/protein.links.v10.5/9606.protein.links.v10.5.txt.gz
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
        self.ELM = [(entry['start'], entry['stop'], entry['elm_identifier']) for entry in data if
                    entry['is_filtered'] == 'FALSE' or entry['is_filtered'] == 'False']

    def query_ELM_for_mutant(self):
        raise NotImplementedError
        # todo implement better.
        # Variant
        mfile = os.path.join(self.settings.ELM_variant_folder,
                             self.uniprot + '_' + self.file_friendly_mutation + '_ELM_variant.tsv')
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
        """EMBL ids should have come from the Unirpto entry. However, in some cases (too many) it is absent."""
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
            self.ELM = [(entry['start'], entry['stop'], entry['elm_identifier']) for entry in data if
                        entry['is_filtered'] == 'FALSE' or entry['is_filtered'] == 'False']
        return self

    @_failsafe
    def parse_pLI(self):
        for line in csv.DictReader(self.settings.open('ExAC_pLI'), delimiter='\t'):
            # transcript	gene	chr	n_exons	cds_start	cds_end	bp	mu_syn	mu_mis	mu_lof	n_syn	n_mis	n_lof	exp_syn	exp_mis	exp_lof	syn_z	mis_z	lof_z	pLI	pRec	pNull
            if self.gene_name == line['gene']:
                self.pLI = float(line[
                                     'pLI'])  # intolerant of a single loss-of-function variant (like haploinsufficient genes, observed ~ 0.1*expected)
                self.pRec = float(line[
                                      'pRec'])  # intolerant of two loss-of-function variants (like recessive genes, observed ~ 0.5*expected)
                self.pNull = float(
                    line['pNull'])  # completely tolerant of loss-of-function variation (observed = expected)
                self.log('pLI: {pLI}, pRec {pRec}, pNull {pNull}'.format(**line))
                break
        else:
            msg = 'Gene {} not found in ExAC table.'.format(self.gene_name)
            warn(msg)
            self.log(msg)
        return self

    @_failsafe
    def parse_gnomAD(self):
        """
        preparsed. using from michelanglo_protein.generate.split_gnomAD import gnomAD
        which splits the huge gnomAD into lots
        :return:
        """
        #
        #
        # gnomAD().write('gnomAD')
        file = os.path.join(self.settings.temp_folder, 'gnomAD', self.uniprot + '.json')
        self.gnomAD = []
        if os.path.exists(file):
            for snp in json.load(open(file)):
                resi = snp['residue_index']
                if isinstance(resi, str):
                    resi = int(resi.split('-')[0])
                # print('HERE', snp)
                variant = Variant(id=f'gnomAD_{resi}_{resi}_{snp["id"]}',
                                  description='{from_residue}{residue_index}{to_residue} ({id})'.format(**snp),
                                  x=resi, y=resi,
                                  impact=snp['impact'],
                                  homozygous=snp['homozygous'])
                self.gnomAD.append(variant)
        else:
            self.gnomAD = []
            warn('No gnomAD data from {0} {1}? Have your run michelanglo_protein.generate.split_gnomAD?'.format(
                self.gene_name, self.uniprot))
        self.log('gnomAD mutations: {0}'.format(len(self.gnomAD)))
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

    def parse_all(self, mode='parallel'):
        """
        Gets all the data. if running in parallel see self._threads list.
        :param mode: parallel | background (=parallel but not complete) | serial (or anything else)
        :return:
        """
        ## Sanity check
        if not self.uniprot:
            warn('There is no uniprot value in this entry!')
            self.uniprot = json.load(open(os.path.join(self.settings.data_folder, 'human_prot_namedex.json')))[
                self.gene_name]
        if not self.gene_name:
            self.parse_uniprot()  # this runs off the web.
        ## fetch!
        tasks = {
            'Swissmodel': self.parse_swissmodel,
            'pLI': self.parse_pLI,
            'param': self.compute_params,
            'gnomAD': self.parse_gnomAD,
            'ptm': self.get_PTM
            # 'parse_pdb_blast': self.parse_pdb_blast
            # 'manual': self.add_manual_data,
            # 'Binding partners': self.fetch_binders
        }
        # prot.get_percent_modelled()
        # prot.parse_ExAC_type()
        # prot.dump()
        # tasks = {'Uniprot': self.parse_uniprot}
        if mode in ('parallel', 'background'):
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
        else:  # serial
            for task_fn in tasks.values():
                task_fn()
        return self

    @_failsafe
    def parse_pLI(self):
        for line in csv.DictReader(self.settings.open('ExAC_pLI'), delimiter='\t'):
            # transcript	gene	chr	n_exons	cds_start	cds_end	bp	mu_syn	mu_mis	mu_lof	n_syn	n_mis	n_lof	exp_syn	exp_mis	exp_lof	syn_z	mis_z	lof_z	pLI	pRec	pNull
            if self.gene_name == line['gene']:
                self.pLI = float(line[
                                     'pLI'])  # intolerant of a single loss-of-function variant (like haploinsufficient genes, observed ~ 0.1*expected)
                self.pRec = float(line[
                                      'pRec'])  # intolerant of two loss-of-function variants (like recessive genes, observed ~ 0.5*expected)
                self.pNull = float(
                    line['pNull'])  # completely tolerant of loss-of-function variation (observed = expected)
                break
        else:
            warn('Gene {} not found in ExAC table.'.format(self.gene_name))
        return self

    def parse_swissmodel(self):
        """
        Fills the object with ``.swissmodel`` data for key species.

        :return:
        """
        self.swissmodel = []
        # entry: {"uniprot_seq_length": 246, "provider": "PDB", "seqid": "", "from": 3, "uniprot_ac": "P31946",
        # "uniprot_seq_md5": "c82f2efd57f939ee3c4e571708dd31a8", "url": "https://swissmodel.expasy.org/repository/uniprot/P31946.pdb?from=3&to=232&template=6byk&provider=pdb",
        # "to": 232, "template": "6byk", "iso_id": "P31946-1", "coordinate_id": "5be4a9c602efd0e456a7ffeb"}
        # figure what animal it is. fix missing.
        if self.organism['NCBI Taxonomy'] == 'NA':
            self.organism['NCBI Taxonomy'] = self.get_species_for_uniprot()
        # Now. figure what animal it is.
        model_organisms = (9606, 3702, 6239, 7227, 10090, 36329, 83332, 83333, 93061, 190650, 208964, 284812, 559292)
        if int(self.organism['NCBI Taxonomy']) in model_organisms:
            reffile = 'swissmodel' + self.organism['NCBI Taxonomy']
        else:
            return self
        models = json.load(self.settings.open(reffile))['index']
        for model in models:
            if self.uniprot == model['uniprot_ac']:
                if model['provider'] == 'PDB':
                    continue
                if not model['seqid']:
                    warn('Odd entry')
                    model['seqid'] = 0
                self.swissmodel.append(
                    Structure(description='{template} (identity:{seqid:.0f}%)'.format(**model),
                              id=model['coordinate_id'],
                              chain='A',
                              code=model['coordinate_id'],
                              url=model['url'],  #is this wise? this url is junk
                              x=int(model['from']),
                              y=int(model['to']),
                              type='swissmodel'
                              )
                )
        self.log('Swissmodel has {0} models.'.format(len(self.swissmodel)))
        return self

    def get_percent_modelled(self):
        """Returns the [0,1] fraction of the protein length that is modelled (repeats arent recounted)"""

        def clean(text):
            if isinstance(text, int):
                return text
            else:
                return int(text.lstrip().replace('-', ' ').replace(',', ' ').split(' ')[0]) - 1

        try:
            if len(self) == 0:
                return 0
            state = [False for i in range(len(self))]
            for dataset in (self.swissmodel, self.pdbs):
                for model in dataset:
                    for i in range(clean(model.x), clean(model.y)):
                        state[i] = True
            self.percent_modelled = sum(state) / len(self)
        except Exception as err:
            warn(f'error arose in %modelled...{err}')
            self.percent_modelled = -1
        return self.percent_modelled

    def parse_pdb_blast(self):
        file = os.path.join(self.settings.temp_folder, 'blastpdb2', self.uniprot + '.json')
        if os.path.exists(file):
            self.pdb_matches = [Structure(**match) for match in json.load(open(file))]
        else:
            warn('No PDB blast data from {0} {1}?'.format(self.gene_name, self.uniprot))
        return self

    ############ model checks.
    def get_offsets(self):
        for m in self.pdbs:
            m.lookup_sifts()
        return self

    def get_resolutions(self):
        for m in self.pdbs:
            m.lookup_resolution()
        return self

    def augment_pdb_w_offset(self):
        """
        SIFTS data. for PDBe query see elsewhere.
        There are four start/stop pairs that need to be compared to get a good idea of a protein.
        For a lengthy discussion see https://blog.matteoferla.com/2019/09/pdb-numbering-rollercoaster.html

        :return:
        """
        warn('THIS METHOD IS DEPRACTED USE THE METHOD OF STRUCTURE> get_offsets')
        # namedtuple('Structure', ['id', 'description', 'x', 'y', 'url','type','chain','offset', 'extra'])
        for i in range(len(self.pdbs)):
            model = self.pdbs[i]
            details = self.lookup_pdb_chain_uniprot(model.url, model.chain)
            detail = details[
                0]  # offset should be the same for all.. check_discrepancy_in_pdb_chain_uniprot() not needed
            if detail['PDB_BEG'] != detail['SP_BEG'] or detail['PDB_END'] != detail['SP_END']:
                model.offset = int(detail['PDB_BEG']) - int(detail['SP_BEG'])
                self.pdbs[i] = model  # remnant from when it was a new namedtuple
        return self

    # pdb_chain_uniprot.tsv
    def lookup_pdb_chain_uniprot(self, pdb, chain):
        warn('THIS METHOD IS DEPRACTED USE THE METHOD OF STRUCTURE>')
        details = []
        headers = 'PDB     CHAIN   SP_PRIMARY      RES_BEG RES_END PDB_BEG PDB_END SP_BEG  SP_END'.split('\t')
        with self.settings.open('pdb_chain_uniprot') as fh:
            for row in fh:
                if pdb.lower() == row[0:4]:
                    details.append(dict(zip(headers, row.split('\t'))))
        return details

    def check_discrepancy_in_pdb_chain_uniprot(self, details):
        warn('THIS METHOD IS DEPRACTED USE THE METHOD OF STRUCTURE>')
        for detail in details:
            if detail['PDB_BEG'] != detail['SP_BEG'] or detail['PDB_END'] != detail['SP_END']:
                print('Sequence discrepancy.')
                return False
        return True

    def get_structures(self):
        for model in self.pdbs:
            model.get_coordinates()
        return self

    # figure out which is best model
    def get_best_model(self):
        self.correct_definitions()
        # So ideally making a multidomain concatenation would be best. But that is hard as it must not be overlapping spacially and sequentially.
        # figuring out what works best is key
        smodels = sorted(self.pdbs, key=lambda x: int(x['y']) - int(x['x']), reverse=True)
        for smodel in smodels:
            details = self.lookup_pdb_chain_uniprot(smodel['id'].split('_')[0], smodel['id'].split('_')[1])
            if self.check_discrepancy_in_pdb_chain_uniprot(details):
                return smodel
        smodels = sorted(self.swissmodel, key=lambda x: int(x['y']) - int(x['x']), reverse=True)
        if smodels:
            return smodels[0]
        else:
            return None

    @_failsafe
    def get_PTM(self):
        """
        To run this requires ``Phosphosite().split().write()`` to have been run at some point...
        :return:
        """
        self.features['PSP_modified_residues'] = []
        assert self.uniprot, 'Uniprot Acc. required. Kind of.'
        self.log(f'Getting PTM for {self.uniprot}')
        fp = os.path.join(self.settings.temp_folder, 'phosphosite', self.uniprot + '.json')
        if os.path.exists(fp):
            with open(fp) as fh:
                self.features['PSP_modified_residues'] = json.load(fh)
        return self

    @_failsafe
    def _test_failsafe(self):
        raise ValueError('Will failsafe catch it? ({})'.format(self.settings.error_tollerant))

    def compute_params(self):
        self.sequence = self.sequence.replace(' ', '').replace('X', '')
        p = ProtParam.ProteinAnalysis(self.sequence)
        self.properties = {}
        self.properties['kd'] = p.protein_scale(ProtParamData.kd, window=9,
                                                edge=.4)  # Kyte & Doolittle index of hydrophobicity J. Mol. Biol. 157:105-132(1982).
        self.properties['Flex'] = p.protein_scale(ProtParamData.Flex, window=9,
                                                  edge=.4)  # Flexibility Normalized flexibility parameters (B-values), average Vihinen M., Torkkila E., Riikonen P. Proteins. 19(2):141-9(1994).
        self.properties['hw'] = p.protein_scale(ProtParamData.hw, window=9,
                                                edge=.4)  # Hydrophilicity Hopp & Wood Proc. Natl. Acad. Sci. U.S.A. 78:3824-3828(1981)
        self.properties['em'] = p.protein_scale(ProtParamData.em, window=9,
                                                edge=.4)  # Surface accessibility Vergoten G & Theophanides T, Biomolecular Structure and Dynamics, pg.138 (1997).
        self.properties['ja'] = p.protein_scale(ProtParamData.ja, window=9,
                                                edge=.4)  # Janin Interior to surface transfer energy scale
        # DIWV requires a mod.
        return self
