from protein._protein_base_mixin import _BaseMixin
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
