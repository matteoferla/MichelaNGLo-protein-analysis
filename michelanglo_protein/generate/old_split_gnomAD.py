__doc__ = """
This script is actually independent of settings.
Entry point: ``gnomAD().split().write()``
"""

import os, time
import gzip
import json

from pprint import PrettyPrinter

pprint = PrettyPrinter().pprint

from collections import defaultdict, namedtuple
from typing import Union, Dict, List
from warnings import warn


class gnomADVariant:
    """This is the same as the namedtuple but with more stuff. It does not get written. to_dict does."""

    def __init__(self, symbol: str, identifier: str, from_residue: str, residue_index: Union[str, int],
                 to_residue: str, impact: str, count: int, homozygous: int):
        self.symbol = symbol
        self.id = identifier
        self.from_residue = from_residue
        self.to_residue = to_residue
        self.residue_index = int(residue_index)
        self.impact = impact
        self.homozygous = int(homozygous)
        self.count = int(count)

    @staticmethod
    def parse_line(line: str) -> List[Dict]:
        data = dict(zip(['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'], line.split('\t')))
        # unpack info
        info = {x.split('=')[0]: x.split('=')[1] for x in data['INFO'].split(';') if '=' in x}
        del data['INFO']
        # definitions were taken from the .vep.vcf
        if 'CSQ' in info:
            vep_id = 'CSQ'
        elif 'vep' in info:
            vep_id = 'vep'
        else:
            raise ValueError('Not a VEP.')
        veps = info[vep_id].split(',')
        del info[vep_id]
        get_csq = lambda entry: dict(zip(['Allele',
                                          'Consequence',
                                          'IMPACT',
                                          'SYMBOL',
                                          'Gene',
                                          'Feature_type',
                                          'Feature',
                                          'BIOTYPE',
                                          'EXON',
                                          'INTRON',
                                          'HGVSc',
                                          'HGVSp',
                                          'cDNA_position',
                                          'CDS_position',
                                          'Protein_position',
                                          'Amino_acids',
                                          'Codons',
                                          'Existing_variation',
                                          'ALLELE_NUM',
                                          'DISTANCE',
                                          'STRAND',
                                          'FLAGS',
                                          'VARIANT_CLASS',
                                          'MINIMISED',
                                          'SYMBOL_SOURCE',
                                          'HGNC_ID',
                                          'CANONICAL',
                                          'TSL',
                                          'APPRIS',
                                          'CCDS',
                                          'ENSP',
                                          'SWISSPROT',
                                          'TREMBL',
                                          'UNIPARC',
                                          'GENE_PHENO',
                                          'SIFT',
                                          'PolyPhen',
                                          'DOMAINS',
                                          'HGVS_OFFSET',
                                          'GMAF',
                                          'AFR_MAF',
                                          'AMR_MAF',
                                          'EAS_MAF',
                                          'EUR_MAF',
                                          'SAS_MAF',
                                          'AA_MAF',
                                          'EA_MAF',
                                          'ExAC_MAF',
                                          'ExAC_Adj_MAF',
                                          'ExAC_AFR_MAF',
                                          'ExAC_AMR_MAF',
                                          'ExAC_EAS_MAF',
                                          'ExAC_FIN_MAF',
                                          'ExAC_NFE_MAF',
                                          'ExAC_OTH_MAF',
                                          'ExAC_SAS_MAF',
                                          'CLIN_SIG',
                                          'SOMATIC',
                                          'PHENO',
                                          'PUBMED',
                                          'MOTIF_NAME',
                                          'MOTIF_POS',
                                          'HIGH_INF_POS',
                                          'MOTIF_SCORE_CHANGE',
                                          'LoF',
                                          'LoF_filter',
                                          'LoF_flags',
                                          'LoF_info',
                                          'context',
                                          'ancestral']
                                         , entry.split('|')))
        return [{**data, **info, **get_csq(entry)} for entry in veps]

    def to_dict(self) -> Dict:
        return self.__dict__

    @classmethod
    def from_line(cls, line: str):
        """
        This takes one millisecond per line.
        :param line:
        :return:
        """
        multidata = cls.parse_line(line)
        variant = []
        for data in multidata:
            if data['Consequence'] in ('missense_variant', 'stop_gained', 'stop_lost', 'frameshift_variant'):
                if not data['Amino_acids'].strip():
                    continue
                if data['CANONICAL'] == 'YES' and data['FILTER'] == 'PASS':
                    variant.append(cls(symbol=data['SYMBOL'],
                                       identifier=data['ID'],
                                       from_residue=data['Amino_acids'].split('/')[0],
                                       residue_index=int(data['Protein_position'].split('-')[0]),
                                       to_residue=data['Amino_acids'].split('/')[1] if '/' in data[
                                           'Amino_acids'] else 'X',
                                       impact=data['IMPACT'],
                                       homozygous=int(data['nhomalt']),
                                       count=int(data['AC'])))
            else:
                pass
        return variant


class gnomAD:

    def __init__(self,
                 masterfiles: Union[str, List],
                 namedexfile: str,
                 folder: str = 'gNOMAD',
                 store_in_memory: bool = False):
        """
        Instantiation starts the settings. But the settings can be changed.
        ``split`` splits the file into the self.data dict containing gene acc id as key and list of gnomADVariant.
        But the bound method `write` writes and the gnomADVariant as regular dictionary.

        :param masterfiles: path to ``gnomad.genomes.r2.1.1.exome_calling_intervals.sites.vcf.bgz`` or
                            ``gnomad.exomes.r2.1.1.sites.vcf.bgz`` from gnomAD
        :type masterfiles: Union[str,List]
        :param namedexfile: path to json generated by ``UniprotMasterReader`` which all the gene synomyms to Uniprot id.
        :type namedexfile: str
        :param folder: Where to save. Internally this is ``self.datafolder``.
        :type folder: str
        :param store_in_memory: If True ``self.data`` will contain all the data.
        """
        self._namedex = json.load(open(namedexfile))
        self.data = defaultdict(list)  #: the data parsed. It is cleared if store_in_memory is False
        if isinstance(masterfiles, str):
            masterfiles = [masterfiles]
        self.masterfiles = masterfiles
        if folder:
            if not os.path.exists(folder):
                os.mkdir(folder)
            self.datafolder = folder
        else:
            self.datafolder = None
            assert store_in_memory, 'If you don\'t want to save the files, store_in_memory has to be true.'
        self.store_in_memory = store_in_memory

    def save_entry(self, previous):
        path = os.path.join(self.datafolder, previous + '.json')
        if os.path.exists(path) and not self.store_in_memory:
            print('... preexistant!')
            with open(path, 'r') as r:
                older = json.load(r)
        else:
            older = []
        with open(path, 'w') as w:
            json.dump([v.to_dict() for v in self.data[previous]] + older, w)
        if self.store_in_memory == False:
            self.data[previous] = []

    def split(self):
        """
        This takes a whole 30 per gene.
        This means that the 58 GB file
        ``gnomad.exomes.r2.1.1.sites.vcf.bgz`` takes 9 days.
        :return:
        """
        for masterfile in self.masterfiles:
            print(masterfile)
            uniprot = ''
            tick = time.time()
            n = 0  # line count for debug.
            with gzip.open(masterfile, 'rt') as f:
                previous = ''
                for line in f:
                    if n % 100 == 0:
                        # print(n, line[:50])
                        pass
                    n += 1
                    if line[0] != '#':
                        for variant in gnomADVariant.from_line(line):
                            if variant and variant.symbol in self._namedex:
                                uniprot = self._namedex[variant.symbol]
                                # print(variant.symbol, '? This gene is', uniprot)
                                if self.datafolder and previous and uniprot != previous:
                                    tock = time.time()
                                    print(f'GENE CHANGE from {previous} TO {uniprot}. ' +
                                          f'store:{self.store_in_memory} time:{tock - tick}s')
                                    tick = tock
                                    self.save_entry(previous)
                                if variant.id in [v.id for v in self.data[uniprot]]:
                                    continue  # dejavu
                                else:
                                    self.data[uniprot].append(variant)
                                previous = uniprot
                            else:
                                if variant:
                                    warn(f'This line has a mystery gene {variant.symbol}!')
                    else:
                        pass
                        # print(line)
            if uniprot:
                self.save_entry(uniprot)
            else:
                warn(f'{masterfile} appears to be empty')
            print('Complete!')
        return self

    def write(self):
        """
        This saves self.data. Makes sense only if self.store_in_memory is True.
        """
        assert not self.store_in_memory, 'Store in memory is off?!'
        for gene in self.data:
            with open(os.path.join(self.datafolder, gene + '.json'), 'w') as w:
                json.dump([v.to_dict() for v in self.data[gene]], w)
        return self
