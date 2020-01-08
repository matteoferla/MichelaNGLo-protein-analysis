import re, os
from warnings import warn

class Mutation:
    """
    Stores the mutation. Not to be confused with the namedtuple Variant, which stores gnomAD mutations.
    >>> Mutation('p.M12D')
    >>> Mutation('M12D')
    >>> Mutation(gnomAD_variant_instance)
    This class does not do analyses with Protein, but ProteinAnalysis do. Here however, wordy conversions happen.
    """
    # the following variable was made in apriori_effect.py
    _apriori_data = { 'A>A': 'identical',
                'A>C': 'bigger|more polar',
                'A>D': 'differently charged|bigger',
                'A>E': 'differently charged|bigger',
                'A>F': 'bigger',
                'A>G': 'more polar|smaller|more flexible',
                'A>H': 'differently charged|bigger|from a non-aromatic to an aromatic',
                'A>I': 'bigger',
                'A>K': 'differently charged|bigger',
                'A>L': 'bigger',
                'A>M': 'bigger',
                'A>N': 'bigger|more polar',
                'A>P': 'bigger|from a non-aromatic to an aromatic|more flexible|more polar',
                'A>Q': 'bigger|more polar',
                'A>R': 'differently charged|bigger',
                'A>S': 'bigger|more flexible|more polar',
                'A>T': 'bigger|more polar',
                'A>V': 'bigger',
                'A>W': 'bigger|from a non-aromatic to an aromatic',
                'A>Y': 'bigger|from a non-aromatic to an aromatic',
                'C>A': 'smaller|more hydrophobic',
                'C>C': 'identical',
                'C>D': 'differently charged|bigger',
                'C>E': 'differently charged|bigger',
                'C>F': 'bigger|more hydrophobic',
                'C>G': 'smaller|more flexible',
                'C>H': 'differently charged|bigger|from a non-aromatic to an aromatic',
                'C>I': 'bigger|more hydrophobic',
                'C>K': 'differently charged|bigger',
                'C>L': 'bigger|more hydrophobic',
                'C>M': 'bigger|more hydrophobic',
                'C>N': 'bigger',
                'C>P': 'from a non-aromatic to an aromatic|smaller|more flexible',
                'C>Q': 'bigger',
                'C>R': 'differently charged|bigger',
                'C>S': 'bigger|smaller|equally sized|more flexible',
                'C>T': 'bigger',
                'C>V': 'bigger|more hydrophobic',
                'C>W': 'bigger|from a non-aromatic to an aromatic|more hydrophobic',
                'C>Y': 'bigger|from a non-aromatic to an aromatic|more hydrophobic',
                'D>A': 'differently charged|smaller',
                'D>C': 'differently charged|smaller',
                'D>D': 'identical',
                'D>E': 'bigger',
                'D>F': 'differently charged|bigger',
                'D>G': 'differently charged|smaller|more flexible',
                'D>H': 'differently charged|from a non-aromatic to an aromatic|differently shaped',
                'D>I': 'differently charged|differently shaped',
                'D>K': 'differently charged|bigger',
                'D>L': 'differently charged|equally sized',
                'D>M': 'differently charged|differently shaped',
                'D>N': 'differently charged|equally sized',
                'D>P': 'differently charged|from a non-aromatic to an aromatic|smaller|more flexible',
                'D>Q': 'differently charged|bigger',
                'D>R': 'differently charged|bigger',
                'D>S': 'differently charged|smaller|more flexible',
                'D>T': 'differently charged|differently shaped',
                'D>V': 'differently charged|differently shaped',
                'D>W': 'differently charged|bigger|from a non-aromatic to an aromatic',
                'D>Y': 'differently charged|bigger|from a non-aromatic to an aromatic',
                'E>A': 'differently charged|smaller',
                'E>C': 'differently charged|smaller',
                'E>D': 'smaller',
                'E>E': 'identical',
                'E>F': 'differently charged|bigger',
                'E>G': 'differently charged|smaller|more flexible',
                'E>H': 'differently charged|from a non-aromatic to an aromatic|differently shaped',
                'E>I': 'differently charged|differently shaped',
                'E>K': 'differently charged|bigger',
                'E>L': 'differently charged|differently shaped',
                'E>M': 'differently charged|smaller',
                'E>N': 'differently charged|smaller',
                'E>P': 'differently charged|from a non-aromatic to an aromatic|smaller|more flexible',
                'E>Q': 'differently charged|equally sized',
                'E>R': 'differently charged|bigger',
                'E>S': 'differently charged|smaller|more flexible',
                'E>T': 'differently charged|differently shaped',
                'E>V': 'differently charged|differently shaped',
                'E>W': 'differently charged|from a non-aromatic to an aromatic|differently shaped',
                'E>Y': 'differently charged|bigger|from a non-aromatic to an aromatic',
                'F>A': 'smaller',
                'F>C': 'more polar|smaller',
                'F>D': 'differently charged|smaller',
                'F>E': 'differently charged|smaller',
                'F>F': 'identical',
                'F>G': 'more polar|smaller|more flexible',
                'F>H': 'differently charged|from a non-aromatic to an aromatic|smaller',
                'F>I': 'differently shaped',
                'F>K': 'differently charged|differently shaped',
                'F>L': 'smaller',
                'F>M': 'smaller',
                'F>N': 'more polar|differently shaped',
                'F>P': 'from a non-aromatic to an aromatic|smaller|more flexible|more polar',
                'F>Q': 'more polar|differently shaped',
                'F>R': 'differently charged|differently shaped',
                'F>S': 'more polar|smaller|more flexible',
                'F>T': 'more polar|smaller',
                'F>V': 'smaller',
                'F>W': 'from a non-aromatic to an aromatic|differently shaped',
                'F>Y': 'bigger|from a non-aromatic to an aromatic',
                'G>A': 'bigger|more rigid|more hydrophobic',
                'G>C': 'bigger|more rigid',
                'G>D': 'differently charged|bigger|more rigid',
                'G>E': 'differently charged|bigger|more rigid',
                'G>F': 'bigger|more rigid|more hydrophobic',
                'G>G': 'identical',
                'G>H': 'differently charged|bigger|from a non-aromatic to an aromatic|more rigid',
                'G>I': 'bigger|more rigid|more hydrophobic',
                'G>K': 'differently charged|bigger|more rigid',
                'G>L': 'bigger|more rigid|more hydrophobic',
                'G>M': 'bigger|more rigid|more hydrophobic',
                'G>N': 'bigger|more rigid',
                'G>P': 'bigger|from a non-aromatic to an aromatic|more flexible',
                'G>Q': 'bigger|more rigid',
                'G>R': 'differently charged|bigger|more rigid',
                'G>S': 'bigger|more rigid',
                'G>T': 'bigger|more rigid',
                'G>V': 'bigger|more rigid|more hydrophobic',
                'G>W': 'bigger|from a non-aromatic to an aromatic|more rigid|more hydrophobic',
                'G>Y': 'bigger|from a non-aromatic to an aromatic|more rigid|more hydrophobic',
                'H>A': 'differently charged|from an aromatic to a non-aromatic|smaller',
                'H>C': 'differently charged|from an aromatic to a non-aromatic|smaller',
                'H>D': 'differently charged|from an aromatic to a non-aromatic|differently shaped',
                'H>E': 'differently charged|from an aromatic to a non-aromatic|differently shaped',
                'H>F': 'differently charged|bigger|from an aromatic to a non-aromatic',
                'H>G': 'differently charged|from an aromatic to a non-aromatic|smaller|more flexible',
                'H>H': 'identical',
                'H>I': 'differently charged|from an aromatic to a non-aromatic|differently shaped',
                'H>K': 'from an aromatic to a non-aromatic|differently shaped',
                'H>L': 'differently charged|from an aromatic to a non-aromatic|differently shaped',
                'H>M': 'differently charged|from an aromatic to a non-aromatic|differently shaped',
                'H>N': 'differently charged|from an aromatic to a non-aromatic|differently shaped',
                'H>P': 'differently charged|smaller|more flexible',
                'H>Q': 'differently charged|from an aromatic to a non-aromatic|differently shaped',
                'H>R': 'from an aromatic to a non-aromatic|differently shaped',
                'H>S': 'differently charged|from an aromatic to a non-aromatic|smaller|more flexible',
                'H>T': 'differently charged|from an aromatic to a non-aromatic|differently shaped',
                'H>V': 'differently charged|from an aromatic to a non-aromatic|differently shaped',
                'H>W': 'differently charged|bigger',
                'H>Y': 'differently charged|bigger',
                'I>A': 'smaller',
                'I>C': 'more polar|smaller',
                'I>D': 'differently charged|differently shaped',
                'I>E': 'differently charged|differently shaped',
                'I>F': 'differently shaped',
                'I>G': 'more polar|smaller|more flexible',
                'I>H': 'differently charged|from a non-aromatic to an aromatic|differently shaped',
                'I>I': 'identical',
                'I>K': 'differently charged|differently shaped',
                'I>L': 'differently shaped',
                'I>M': 'differently shaped',
                'I>N': 'more polar|differently shaped',
                'I>P': 'from a non-aromatic to an aromatic|smaller|more flexible|more polar',
                'I>Q': 'more polar|differently shaped',
                'I>R': 'differently charged|differently shaped',
                'I>S': 'more polar|smaller|more flexible',
                'I>T': 'more polar|differently shaped',
                'I>V': 'smaller',
                'I>W': 'from a non-aromatic to an aromatic|differently shaped',
                'I>Y': 'from a non-aromatic to an aromatic|differently shaped',
                'K>A': 'differently charged|smaller',
                'K>C': 'differently charged|smaller',
                'K>D': 'differently charged|smaller',
                'K>E': 'differently charged|smaller',
                'K>F': 'differently charged|differently shaped',
                'K>G': 'differently charged|smaller|more flexible',
                'K>H': 'from a non-aromatic to an aromatic|differently shaped',
                'K>I': 'differently charged|differently shaped',
                'K>K': 'identical',
                'K>L': 'differently charged|differently shaped',
                'K>M': 'differently charged|smaller',
                'K>N': 'differently charged|differently shaped',
                'K>P': 'differently charged|from a non-aromatic to an aromatic|smaller|more flexible',
                'K>Q': 'differently charged|differently shaped',
                'K>R': 'bigger',
                'K>S': 'differently charged|smaller|more flexible',
                'K>T': 'differently charged|differently shaped',
                'K>V': 'differently charged|differently shaped',
                'K>W': 'differently charged|from a non-aromatic to an aromatic|differently shaped',
                'K>Y': 'differently charged|from a non-aromatic to an aromatic|differently shaped',
                'L>A': 'smaller',
                'L>C': 'more polar|smaller',
                'L>D': 'differently charged|equally sized',
                'L>E': 'differently charged|differently shaped',
                'L>F': 'bigger',
                'L>G': 'more polar|smaller|more flexible',
                'L>H': 'differently charged|from a non-aromatic to an aromatic|differently shaped',
                'L>I': 'differently shaped',
                'L>K': 'differently charged|differently shaped',
                'L>L': 'identical',
                'L>M': 'differently shaped',
                'L>N': 'more polar|equally sized',
                'L>P': 'from a non-aromatic to an aromatic|smaller|more flexible|more polar',
                'L>Q': 'more polar|differently shaped',
                'L>R': 'differently charged|differently shaped',
                'L>S': 'more polar|smaller|more flexible',
                'L>T': 'more polar|differently shaped',
                'L>V': 'smaller',
                'L>W': 'bigger|from a non-aromatic to an aromatic',
                'L>Y': 'bigger|from a non-aromatic to an aromatic',
                'M>A': 'smaller',
                'M>C': 'more polar|smaller',
                'M>D': 'differently charged|differently shaped',
                'M>E': 'differently charged|bigger',
                'M>F': 'bigger',
                'M>G': 'more polar|smaller|more flexible',
                'M>H': 'differently charged|from a non-aromatic to an aromatic|differently shaped',
                'M>I': 'differently shaped',
                'M>K': 'differently charged|bigger',
                'M>L': 'differently shaped',
                'M>M': 'identical',
                'M>N': 'more polar|differently shaped',
                'M>P': 'from a non-aromatic to an aromatic|smaller|more flexible|more polar',
                'M>Q': 'bigger|more polar',
                'M>R': 'differently charged|bigger',
                'M>S': 'more polar|smaller|more flexible',
                'M>T': 'more polar|differently shaped',
                'M>V': 'differently shaped',
                'M>W': 'from a non-aromatic to an aromatic|differently shaped',
                'M>Y': 'bigger|from a non-aromatic to an aromatic',
                'N>A': 'smaller|more hydrophobic',
                'N>C': 'smaller',
                'N>D': 'differently charged|equally sized',
                'N>E': 'differently charged|bigger',
                'N>F': 'differently shaped|more hydrophobic',
                'N>G': 'smaller|more flexible',
                'N>H': 'differently charged|from a non-aromatic to an aromatic|differently shaped',
                'N>I': 'differently shaped|more hydrophobic',
                'N>K': 'differently charged|differently shaped',
                'N>L': 'equally sized|more hydrophobic',
                'N>M': 'differently shaped|more hydrophobic',
                'N>N': 'identical',
                'N>P': 'from a non-aromatic to an aromatic|smaller|more flexible',
                'N>Q': 'bigger',
                'N>R': 'differently charged|differently shaped',
                'N>S': 'smaller|more flexible',
                'N>T': 'differently shaped',
                'N>V': 'differently shaped|more hydrophobic',
                'N>W': 'from a non-aromatic to an aromatic|differently shaped|more hydrophobic',
                'N>Y': 'from a non-aromatic to an aromatic|differently shaped|more hydrophobic',
                'P>A': 'from an aromatic to a non-aromatic|smaller|more rigid|more hydrophobic',
                'P>C': 'bigger|from an aromatic to a non-aromatic|more rigid',
                'P>D': 'differently charged|bigger|from an aromatic to a non-aromatic|more rigid',
                'P>E': 'differently charged|bigger|from an aromatic to a non-aromatic|more rigid',
                'P>F': 'bigger|from an aromatic to a non-aromatic|more rigid|more hydrophobic',
                'P>G': 'from an aromatic to a non-aromatic|smaller|more rigid',
                'P>H': 'differently charged|bigger|more rigid',
                'P>I': 'bigger|from an aromatic to a non-aromatic|more rigid|more hydrophobic',
                'P>K': 'differently charged|bigger|from an aromatic to a non-aromatic|more rigid',
                'P>L': 'bigger|from an aromatic to a non-aromatic|more rigid|more hydrophobic',
                'P>M': 'bigger|from an aromatic to a non-aromatic|more rigid|more hydrophobic',
                'P>N': 'bigger|from an aromatic to a non-aromatic|more rigid',
                'P>P': 'identical',
                'P>Q': 'bigger|from an aromatic to a non-aromatic|more rigid',
                'P>R': 'differently charged|bigger|from an aromatic to a non-aromatic|more rigid',
                'P>S': 'bigger|from an aromatic to a non-aromatic|more rigid',
                'P>T': 'bigger|from an aromatic to a non-aromatic|more rigid',
                'P>V': 'bigger|from an aromatic to a non-aromatic|more rigid|more hydrophobic',
                'P>W': 'bigger|more rigid|more hydrophobic',
                'P>Y': 'bigger|more rigid|more hydrophobic',
                'Q>A': 'smaller|more hydrophobic',
                'Q>C': 'smaller',
                'Q>D': 'differently charged|smaller',
                'Q>E': 'differently charged|equally sized',
                'Q>F': 'differently shaped|more hydrophobic',
                'Q>G': 'smaller|more flexible',
                'Q>H': 'differently charged|from a non-aromatic to an aromatic|differently shaped',
                'Q>I': 'differently shaped|more hydrophobic',
                'Q>K': 'differently charged|differently shaped',
                'Q>L': 'differently shaped|more hydrophobic',
                'Q>M': 'smaller|more hydrophobic',
                'Q>N': 'smaller',
                'Q>P': 'from a non-aromatic to an aromatic|smaller|more flexible',
                'Q>Q': 'identical',
                'Q>R': 'differently charged|differently shaped',
                'Q>S': 'smaller|more flexible',
                'Q>T': 'differently shaped',
                'Q>V': 'differently shaped|more hydrophobic',
                'Q>W': 'from a non-aromatic to an aromatic|differently shaped|more hydrophobic',
                'Q>Y': 'from a non-aromatic to an aromatic|differently shaped|more hydrophobic',
                'R>A': 'differently charged|smaller',
                'R>C': 'differently charged|smaller',
                'R>D': 'differently charged|smaller',
                'R>E': 'differently charged|smaller',
                'R>F': 'differently charged|differently shaped',
                'R>G': 'differently charged|smaller|more flexible',
                'R>H': 'from a non-aromatic to an aromatic|differently shaped',
                'R>I': 'differently charged|differently shaped',
                'R>K': 'smaller',
                'R>L': 'differently charged|differently shaped',
                'R>M': 'differently charged|smaller',
                'R>N': 'differently charged|differently shaped',
                'R>P': 'differently charged|from a non-aromatic to an aromatic|smaller|more flexible',
                'R>Q': 'differently charged|differently shaped',
                'R>R': 'identical',
                'R>S': 'differently charged|smaller|more flexible',
                'R>T': 'differently charged|differently shaped',
                'R>V': 'differently charged|differently shaped',
                'R>W': 'differently charged|from a non-aromatic to an aromatic|differently shaped',
                'R>Y': 'differently charged|from a non-aromatic to an aromatic|differently shaped',
                'S>A': 'smaller|more rigid|more hydrophobic',
                'S>C': 'bigger|smaller|equally sized|more rigid',
                'S>D': 'differently charged|bigger|more rigid',
                'S>E': 'differently charged|bigger|more rigid',
                'S>F': 'bigger|more rigid|more hydrophobic',
                'S>G': 'smaller|more flexible',
                'S>H': 'differently charged|bigger|from a non-aromatic to an aromatic|more rigid',
                'S>I': 'bigger|more rigid|more hydrophobic',
                'S>K': 'differently charged|bigger|more rigid',
                'S>L': 'bigger|more rigid|more hydrophobic',
                'S>M': 'bigger|more rigid|more hydrophobic',
                'S>N': 'bigger|more rigid',
                'S>P': 'from a non-aromatic to an aromatic|smaller|more flexible',
                'S>Q': 'bigger|more rigid',
                'S>R': 'differently charged|bigger|more rigid',
                'S>S': 'identical',
                'S>T': 'bigger|more rigid',
                'S>V': 'bigger|more rigid|more hydrophobic',
                'S>W': 'bigger|from a non-aromatic to an aromatic|more rigid|more hydrophobic',
                'S>Y': 'bigger|from a non-aromatic to an aromatic|more rigid|more hydrophobic',
                'T>A': 'smaller|more hydrophobic',
                'T>C': 'smaller',
                'T>D': 'differently charged|differently shaped',
                'T>E': 'differently charged|differently shaped',
                'T>F': 'bigger|more hydrophobic',
                'T>G': 'smaller|more flexible',
                'T>H': 'differently charged|from a non-aromatic to an aromatic|differently shaped',
                'T>I': 'differently shaped|more hydrophobic',
                'T>K': 'differently charged|differently shaped',
                'T>L': 'differently shaped|more hydrophobic',
                'T>M': 'differently shaped|more hydrophobic',
                'T>N': 'differently shaped',
                'T>P': 'from a non-aromatic to an aromatic|smaller|more flexible',
                'T>Q': 'differently shaped',
                'T>R': 'differently charged|differently shaped',
                'T>S': 'smaller|more flexible',
                'T>T': 'identical',
                'T>V': 'equally sized|more hydrophobic',
                'T>W': 'bigger|from a non-aromatic to an aromatic|more hydrophobic',
                'T>Y': 'bigger|from a non-aromatic to an aromatic|more hydrophobic',
                'V>A': 'smaller',
                'V>C': 'more polar|smaller',
                'V>D': 'differently charged|differently shaped',
                'V>E': 'differently charged|differently shaped',
                'V>F': 'bigger',
                'V>G': 'more polar|smaller|more flexible',
                'V>H': 'differently charged|from a non-aromatic to an aromatic|differently shaped',
                'V>I': 'bigger',
                'V>K': 'differently charged|differently shaped',
                'V>L': 'bigger',
                'V>M': 'differently shaped',
                'V>N': 'more polar|differently shaped',
                'V>P': 'from a non-aromatic to an aromatic|smaller|more flexible|more polar',
                'V>Q': 'more polar|differently shaped',
                'V>R': 'differently charged|differently shaped',
                'V>S': 'more polar|smaller|more flexible',
                'V>T': 'more polar|equally sized',
                'V>V': 'identical',
                'V>W': 'bigger|from a non-aromatic to an aromatic',
                'V>Y': 'bigger|from a non-aromatic to an aromatic',
                'W>A': 'from an aromatic to a non-aromatic|smaller',
                'W>C': 'more polar|from an aromatic to a non-aromatic|smaller',
                'W>D': 'differently charged|from an aromatic to a non-aromatic|smaller',
                'W>E': 'differently charged|from an aromatic to a non-aromatic|differently shaped',
                'W>F': 'from an aromatic to a non-aromatic|differently shaped',
                'W>G': 'more polar|from an aromatic to a non-aromatic|smaller|more flexible',
                'W>H': 'differently charged|smaller',
                'W>I': 'from an aromatic to a non-aromatic|differently shaped',
                'W>K': 'differently charged|from an aromatic to a non-aromatic|differently shaped',
                'W>L': 'from an aromatic to a non-aromatic|smaller',
                'W>M': 'from an aromatic to a non-aromatic|differently shaped',
                'W>N': 'more polar|from an aromatic to a non-aromatic|differently shaped',
                'W>P': 'more polar|smaller|more flexible',
                'W>Q': 'more polar|from an aromatic to a non-aromatic|differently shaped',
                'W>R': 'differently charged|from an aromatic to a non-aromatic|differently shaped',
                'W>S': 'more polar|from an aromatic to a non-aromatic|smaller|more flexible',
                'W>T': 'more polar|from an aromatic to a non-aromatic|smaller',
                'W>V': 'from an aromatic to a non-aromatic|smaller',
                'W>W': 'identical',
                'W>Y': 'differently shaped',
                'Y>A': 'from an aromatic to a non-aromatic|smaller',
                'Y>C': 'more polar|from an aromatic to a non-aromatic|smaller',
                'Y>D': 'differently charged|from an aromatic to a non-aromatic|smaller',
                'Y>E': 'differently charged|from an aromatic to a non-aromatic|smaller',
                'Y>F': 'from an aromatic to a non-aromatic|smaller',
                'Y>G': 'more polar|from an aromatic to a non-aromatic|smaller|more flexible',
                'Y>H': 'differently charged|smaller',
                'Y>I': 'from an aromatic to a non-aromatic|differently shaped',
                'Y>K': 'differently charged|from an aromatic to a non-aromatic|differently shaped',
                'Y>L': 'from an aromatic to a non-aromatic|smaller',
                'Y>M': 'from an aromatic to a non-aromatic|smaller',
                'Y>N': 'more polar|from an aromatic to a non-aromatic|differently shaped',
                'Y>P': 'more polar|smaller|more flexible',
                'Y>Q': 'more polar|from an aromatic to a non-aromatic|differently shaped',
                'Y>R': 'differently charged|from an aromatic to a non-aromatic|differently shaped',
                'Y>S': 'more polar|from an aromatic to a non-aromatic|smaller|more flexible',
                'Y>T': 'more polar|from an aromatic to a non-aromatic|smaller',
                'Y>V': 'from an aromatic to a non-aromatic|smaller',
                'Y>W': 'differently shaped',
                'Y>Y': 'identical'}

    names =(('A', 'Ala', 'Alanine'),
            ('B', 'Asx', 'Aspartate/asparagine'),
            ('C', 'Cys', 'Cysteine'),
            ('D', 'Asp', 'Aspartate'),
            ('E', 'Glu', 'Glutamate'),
            ('F', 'Phe', 'Phenylanine'),
            ('G', 'Gly', 'Glycine'),
            ('H', 'His', 'Histidine'),
            ('I', 'Ile', 'Isoleucine'),
            ('K', 'Lys', 'Lysine'),
            ('L', 'Leu', 'Leucine'),
            ('M', 'Met', 'Methionine'),
            ('N', 'Asn', 'Asparagine'),
            ('P', 'Pro', 'Proline'),
            ('Q', 'Gln', 'Glutamine'),
            ('R', 'Arg', 'Arginine'),
            ('S', 'Ser', 'Serine'),
            ('T', 'Thr', 'Threonine'),
            ('U', 'Sel', 'Selenocysteine'),
            ('V', 'Val', 'Valine'),
            ('W', 'Trp', 'Tryptophan'),
            ('X', 'Xaa', 'Any'),
            ('Y', 'Tyr', 'Tyrosine'),
            ('Z', 'Glx', 'Glutamate/glutamine'),
            ('*','Stop', 'Stop'))

    aa_list = tuple('QWERTYIPASDFGHKLCVNM*')


    def __init__(self, mutation=None):
        self.from_residue = ''
        self.to_residue = ''
        self.residue_index = 0
        self.apriori_effect = 'TBD'
        self.surface_expose = ''
        self.clean_mutation = None
        #self.exposure_effect see property.getter exposure_effect
        self.elm = []  # michelanglo_protein.check_elm(mutation) fills it.
        if mutation:
            if not isinstance(mutation,str): #it is the namedtuple `Variant` (gnomAD snp)!
                mutation = mutation.description
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
        for (single, triple, full) in self.names:
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
            else:
                raise ValueError('odd mutation of type deletion' + self.mutation)
        else:
            rex = re.match('(\w)(\d+)(\w)', self.mutation)
            if rex:
                assert rex.group(1) in self.aa_list, 'The from mutant is not a valid amino acid'
                assert rex.group(3) in self.aa_list, 'The to mutant is not a valid amino acid'
                self.from_residue = rex.group(1)
                self.residue_index = int(rex.group(2))
                self.to_residue = rex.group(3)
                self.file_friendly_mutation = self.from_residue + str(self.residue_index) + self.to_residue
                self.clean_mutation = self.from_residue + str(self.residue_index) + self.to_residue
                #self.blossum_score = self.blossum[self.from_residue][self.to_residue]
            else:
                raise ValueError(self.mutation + ' is an odd_mutation')

        ### classify apriori effect
        """
        {'S': 'smaller',
         'B': 'bigger',
         'E': 'equally sized',
         'F': 'more flexible',
         'R': 'more rigid',
         'C': 'differently charged',
         'H': 'more hydrophobic',
         'P': 'more polar',
         'I': 'identical',
         'D': 'differently shaped',
         'a': 'from an aromatic to a non-aromatic',
         'A': 'from a non-aromatic to an aromatic'}"""
        if '*' in self.to_residue:
            self.apriori_effect = 'The mutation results in a truncation.'
        elif '*' in self.from_residue:
            self.apriori_effect = 'The mutation results in a longer sequence.'
        elif len(self.to_residue) > 1:
            self.apriori_effect = 'The mutation results in a frameshift.'
        else:
            self.apriori_effect = 'The mutation changes one amino acid to another that is '+\
                                  ', '.join(self._apriori_data[self.from_residue+'>'+self.to_residue].split('|'))+\
                                  '.'
        return self

    @property
    def exposure_effect(self):
        lowconc = '(lowering protein concentrations)'
        lessbind = '(less binding means less activation or inhibition)'
        neglegible = 'The effect should be negligible.'
        if self.surface_expose == 'buried':
            if self.to_residue == 'P':
                return f'Proline is highly destabilising in protein cores (it cannot be part of &alpha;-helices for example) {lowconc}.'
            elif 'differently charged' in self.apriori_effect:
                return f'Charge changes in protein cores are highly destabilising {lowconc}.'
            elif 'bigger' in self.apriori_effect or 'shaped' in self.apriori_effect:
                return f'Larger residues in protein cores are highly destabilising {lowconc}.'
            elif 'polar' in self.apriori_effect:
                return f'Protein cores are generally hydrophobic, so a change in polarity is generally destabilising {lowconc}.'
            elif 'smaller' in self.apriori_effect:
                return f'Changes to smaller residues remove some interactions, thus weakly destabilising the protein (potentially lowering michelanglo_protein concentrations) but most likely have no effect.'
            else:
                return f'Mutations in protein cores are generally destabilising {lowconc}, but the mutation is very mild.'
        elif self.surface_expose == 'partially buried':
            if self.to_residue == 'P':
                return f'Proline is highly destabilising in protein cores (it cannot be part of &alpha;-helices for example) {lowconc}.'
            elif 'differently charged' in self.apriori_effect:
                return f'Charge changes are most often destabilising {lowconc}.'
            elif 'bigger' in self.apriori_effect or 'shaped' in self.apriori_effect:
                return f'Larger residues are most often destabilising {lowconc}.'
            elif 'polar' in self.apriori_effect:
                return f'A change in polarity is generally destabilising {lowconc} depending how buried it is.'
            elif 'smaller' in self.apriori_effect:
                return f'Changes to smaller residues remove some interactions, but most likely have no effect.'
            else:
                return neglegible
        elif self.surface_expose == 'surface':
            if 'differently charged' in self.apriori_effect:
                return f'A difference in charge on the surface may strongly affect membrane and partner binding {lessbind}.'
            elif 'more hydrophobic' in self.apriori_effect:
                return f'A more hydrophobic residue on the surface may result in aggregation {lowconc}.'
            elif 'bigger' in self.apriori_effect:
                return f'A larger residue on the surface may inhibit binding if it is part of a interface surface {lessbind}.'
            else:
                return neglegible
        else: #no surface_expose value. fill with the ProteinAnalyser
            return '*Could not be calculated*'

    @classmethod
    def long_name(cls, letter):
        """
        Single amino acid letter to a short string with the name spelt out three ways.
        :param letter: 1 AA letter
        :type letter: str
        :return: str
        """
        return ['{n} ({s}, {t})'.format(n=n,s=s,t=t) for s, t, n in cls.names if s == letter][0]


