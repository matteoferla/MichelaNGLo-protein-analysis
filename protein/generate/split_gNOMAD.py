from ..settings_handler import global_settings
import os
import gzip
import json

from pprint import PrettyPrinter
pprint = PrettyPrinter().pprint

from collections import defaultdict


class gNomadVariant:
    def __init__(self, symbol, identifier, from_residue, residue_index, to_residue, impact):
        self.symbol = symbol
        self.id = identifier
        self.from_residue = from_residue
        self.to_residue = to_residue
        self.residue_index = residue_index
        self.impact = impact

    @staticmethod
    def parse_line(line):
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
        csq = dict(zip(
            'Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|GMAF|AFR_MAF|AMR_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|ExAC_MAF|ExAC_Adj_MAF|ExAC_AFR_MAF|ExAC_AMR_MAF|ExAC_EAS_MAF|ExAC_FIN_MAF|ExAC_NFE_MAF|ExAC_OTH_MAF|ExAC_SAS_MAF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info|context|ancestral'.split(
                '|'), info[vep_id].split('|')))
        del info[vep_id]
        return {**data, **info, **csq}

    def to_dict(self):
        return self.__dict__

    @classmethod
    def from_line(cls, line):
        data = cls.parse_line(line)
        if data['Consequence'] in ('missense_variant', 'stop_gained', 'stop_lost', 'frameshift_variant'):
            return cls(symbol=data['SYMBOL'],
                       identifier=data['ID'],
                       from_residue=data['Amino_acids'].split('/')[0],
                       residue_index=data['Protein_position'],
                       to_residue=data['Amino_acids'].split('/')[1],
                       impact=data['IMPACT'])
        else:
            return None




class gNOMAD:
    def __init__(self):
        self.namedex = json.load(open(os.path.join(global_settings.data_folder,'human_prot_namedex.json')))
        self.data = defaultdict(list)
        self.masterfile = os.path.join(global_settings.reference_folder,'gnomad.genomes.r2.1.1.exome_calling_intervals.sites.vcf.bgz')
        with gzip.open(self.masterfile, 'rt') as f:
            for line in f:
                if line[0] != '#':
                    variant = gNomadVariant.from_line(line)
                    if variant and variant.symbol in self.namedex:
                        uniprot = self.namedex[variant.symbol]
                        self.data[uniprot].append(variant)


    def write(self,folder):
        full=os.path.join(global_settings.temp_folder,folder)
        if not os.path.exists(full):
            os.mkdir(full)
        for gene in self.data:
            with open(os.path.join(full,gene+'.json'),'w') as w:
                json.dump([v.to_dict() for v in self.data[gene]],w)


'''
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
        self.alleles = []  # no ExAC.'''
