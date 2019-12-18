from michelanglo_protein.generate._protein_base_mixin import _BaseMixin
_failsafe = _BaseMixin._failsafe

from warnings import warn
import os
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

class _DisusedMixin:
    """
    These are methods that have been superseeded, but may be usedful if someone is running a local version without the generate.
    """

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
    def match_pdb(self):
        warn('The match_pdb method has been replaced. Use this only if you havent generated the data.')
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

    ## depraction zone.
    def write(self, file=None):
        raise Exception('DEPRACATED. use write_uniprot')

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
