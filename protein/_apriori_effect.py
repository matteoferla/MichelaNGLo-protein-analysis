from collections import defaultdict


class Changedex:
    inverter = {'S': 'B', 'B': 'S', 'E': 'E', 'D': 'D',
                'F': 'R', 'R': 'F',
                'C': 'C',
                'I': 'I',
                'A': 'a', 'a': 'A',
                'P': 'H', 'H': 'P'}
    full = {'S': 'smaller',
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
            'A': 'from a non-aromatic to an aromatic'}

    hydrophobic = tuple('AVILMFYW')
    polar = tuple('CGPSTNQ')
    positive = tuple('RHK')
    negative = tuple('DE')
    aa = tuple('AVILMFYWSTNQCGPRHKDE')
    aromatic = tuple('PYWH')

    def __init__(self):
        self._data = defaultdict(set)

    def __getitem__(self, items):
        return self._data['{0}>{1}'.format(*items)]

    def __setitem__(self, items, value):
        for i in items[0]:
            for j in items[1]:
                self._data['{0}>{1}'.format(i, j)].add(value)

    def to_dict(self):
        return {k: '|'.join([self.full[i] for i in self._data[k]]) for k in self._data}

    def __str__(self):
        return str(self.to_dict())

    def fill_inverse(self):
        for a in self.aa:
            for b in self.aa:
                if a != b:
                    for v in self[a, b]:
                        self[b, a] = self.inverter[v]
                else:
                    self._data[a + '>' + a] = set('I')

    def fill(self):
        for a in self.aa:
            self[a, 'G'] = 'S'  ## glycine is smallest
            if a != 'G':  ## glycine followed by alaine in smallness
                self[a, 'A'] = 'S'
                if a != 'A':  ## then serine
                    self[a, 'P'] = 'S'
                    if a != 'P':
                        self[a, 'SC'] = 'S'
                pass
            if a in positive:
                pass

        self['S', 'C'] = 'E'
        self['E', 'Q'] = 'E'
        self['D', 'N'] = 'E'
        self['L', 'DN'] = 'E'
        self['I', 'V'] = 'S'
        self['Y', 'F'] = 'S'
        self['EQ', 'DN'] = 'S'
        self['L', 'V'] = 'S'
        self['FYW', 'H'] = 'S'
        self['W', 'YF'] = 'D'
        self['EQ', 'M'] = 'S'
        self['FY', 'VTLDEM'] = 'S'
        self['W', 'VTLD'] = 'S'
        self['V', 'T'] = 'E'
        self['RK', 'MED'] = 'S'
        self['R', 'K'] = 'S'

        # inverse
        self.fill_inverse()
        for a in self.aa:
            # fill empties
            for b in self.aa:
                if not self[a, b]:
                    self[a, b] = 'D'
            # fill flex
            self[a, 'P'].add('F')  ## proline is most flexible.
            if a != 'P':  ##followed by glycine in flex
                self[a, 'G'].add('F')
                if a not in 'PG':
                    self[a, 'S'].add('F')
            # fill charge
            for charge in (self.positive, self.negative):
                if a in charge:
                    for b in self.aa:
                        if b not in charge:
                            self[a, b] = 'C'
            # fill polar
            if a in self.polar:
                for b in self.hydrophobic:
                    self[a, b] = 'H'
            # fill aroma
            if a in self.aromatic:
                for b in set(self.aa) - set(self.aromatic):
                    self[a, b] = 'a'
            elif a not in self.aromatic:
                for b in self.aromatic:
                    self[a, b] = 'A'

        self.fill_inverse()
        return self

if __name__ == '__main__':
    changedex = Changedex().fill()

    from pprint import PrettyPrinter
    pprint = PrettyPrinter().pprint

    pprint(changedex.to_dict())

    # {'A>G': 'IncreasedPolarity|Smaller|Flexibler', 'A>A': 'Identity', 'V>G': 'IncreasedPolarity|Smaller|Flexibler', 'V>A': 'Smaller', 'V>P': 'toAromatic|Smaller|Flexibler|IncreasedPolarity', 'V>S': 'IncreasedPolarity|Smaller|Flexibler', 'V>C': 'IncreasedPolarity|Smaller', 'I>G': 'IncreasedPolarity|Smaller|Flexibler', 'I>A': 'Smaller', 'I>P': 'toAromatic|Smaller|Flexibler|IncreasedPolarity', 'I>S': 'IncreasedPolarity|Smaller|Flexibler', 'I>C': 'IncreasedPolarity|Smaller', 'L>G': 'IncreasedPolarity|Smaller|Flexibler', 'L>A': 'Smaller', 'L>P': 'toAromatic|Smaller|Flexibler|IncreasedPolarity', 'L>S': 'IncreasedPolarity|Smaller|Flexibler', 'L>C': 'IncreasedPolarity|Smaller', 'M>G': 'IncreasedPolarity|Smaller|Flexibler', 'M>A': 'Smaller', 'M>P': 'toAromatic|Smaller|Flexibler|IncreasedPolarity', 'M>S': 'IncreasedPolarity|Smaller|Flexibler', 'M>C': 'IncreasedPolarity|Smaller', 'F>G': 'IncreasedPolarity|Smaller|Flexibler', 'F>A': 'Smaller', 'F>P': 'toAromatic|Smaller|Flexibler|IncreasedPolarity', 'F>S': 'IncreasedPolarity|Smaller|Flexibler', 'F>C': 'IncreasedPolarity|Smaller', 'Y>G': 'IncreasedPolarity|fromAromatic|Smaller|Flexibler', 'Y>A': 'fromAromatic|Smaller', 'Y>P': 'IncreasedPolarity|Smaller|Flexibler', 'Y>S': 'IncreasedPolarity|fromAromatic|Smaller|Flexibler', 'Y>C': 'IncreasedPolarity|fromAromatic|Smaller', 'W>G': 'IncreasedPolarity|fromAromatic|Smaller|Flexibler', 'W>A': 'fromAromatic|Smaller', 'W>P': 'IncreasedPolarity|Smaller|Flexibler', 'W>S': 'IncreasedPolarity|fromAromatic|Smaller|Flexibler', 'W>C': 'IncreasedPolarity|fromAromatic|Smaller', 'S>G': 'Smaller|Flexibler', 'S>A': 'Smaller|Rigider|IncreasedHydrophobicity', 'S>P': 'toAromatic|Smaller|Flexibler', 'S>S': 'Identity', 'S>C': 'Bigger|Smaller|Equal-sized|Rigider', 'T>G': 'Smaller|Flexibler', 'T>A': 'Smaller|IncreasedHydrophobicity', 'T>P': 'toAromatic|Smaller|Flexibler', 'T>S': 'Smaller|Flexibler', 'T>C': 'Smaller', 'N>G': 'Smaller|Flexibler', 'N>A': 'Smaller|IncreasedHydrophobicity', 'N>P': 'toAromatic|Smaller|Flexibler', 'N>S': 'Smaller|Flexibler', 'N>C': 'Smaller', 'Q>G': 'Smaller|Flexibler', 'Q>A': 'Smaller|IncreasedHydrophobicity', 'Q>P': 'toAromatic|Smaller|Flexibler', 'Q>S': 'Smaller|Flexibler', 'Q>C': 'Smaller', 'C>G': 'Smaller|Flexibler', 'C>A': 'Smaller|IncreasedHydrophobicity', 'C>P': 'toAromatic|Smaller|Flexibler', 'C>S': 'Bigger|Smaller|Equal-sized|Flexibler', 'C>C': 'Identity', 'G>G': 'Identity', 'P>G': 'fromAromatic|Smaller|Rigider', 'P>A': 'fromAromatic|Smaller|Rigider|IncreasedHydrophobicity', 'P>P': 'Identity', 'R>G': 'Charge-changed|Smaller|Flexibler', 'R>A': 'Charge-changed|Smaller', 'R>P': 'Charge-changed|toAromatic|Smaller|Flexibler', 'R>S': 'Charge-changed|Smaller|Flexibler', 'R>C': 'Charge-changed|Smaller', 'H>G': 'Charge-changed|fromAromatic|Smaller|Flexibler', 'H>A': 'Charge-changed|fromAromatic|Smaller', 'H>P': 'Charge-changed|Smaller|Flexibler', 'H>S': 'Charge-changed|fromAromatic|Smaller|Flexibler', 'H>C': 'Charge-changed|fromAromatic|Smaller', 'K>G': 'Charge-changed|Smaller|Flexibler', 'K>A': 'Charge-changed|Smaller', 'K>P': 'Charge-changed|toAromatic|Smaller|Flexibler', 'K>S': 'Charge-changed|Smaller|Flexibler', 'K>C': 'Charge-changed|Smaller', 'D>G': 'Charge-changed|Smaller|Flexibler', 'D>A': 'Charge-changed|Smaller', 'D>P': 'Charge-changed|toAromatic|Smaller|Flexibler', 'D>S': 'Charge-changed|Smaller|Flexibler', 'D>C': 'Charge-changed|Smaller', 'E>G': 'Charge-changed|Smaller|Flexibler', 'E>A': 'Charge-changed|Smaller', 'E>P': 'Charge-changed|toAromatic|Smaller|Flexibler', 'E>S': 'Charge-changed|Smaller|Flexibler', 'E>C': 'Charge-changed|Smaller', 'E>Q': 'Charge-changed|Equal-sized', 'D>N': 'Charge-changed|Equal-sized', 'L>D': 'Charge-changed|Equal-sized', 'L>N': 'IncreasedPolarity|Equal-sized', 'I>V': 'Smaller', 'Y>F': 'fromAromatic|Smaller', 'E>D': 'Smaller', 'E>N': 'Charge-changed|Smaller', 'Q>D': 'Charge-changed|Smaller', 'Q>N': 'Smaller', 'L>V': 'Smaller', 'F>H': 'Charge-changed|toAromatic|Smaller', 'Y>H': 'Charge-changed|Smaller', 'W>H': 'Charge-changed|Smaller', 'W>Y': 'Different-bulk', 'W>F': 'fromAromatic|Different-bulk', 'E>M': 'Charge-changed|Smaller', 'Q>M': 'Smaller|IncreasedHydrophobicity', 'F>V': 'Smaller', 'F>T': 'IncreasedPolarity|Smaller', 'F>L': 'Smaller', 'F>D': 'Charge-changed|Smaller', 'F>E': 'Charge-changed|Smaller', 'F>M': 'Smaller', 'Y>V': 'fromAromatic|Smaller', 'Y>T': 'IncreasedPolarity|fromAromatic|Smaller', 'Y>L': 'fromAromatic|Smaller', 'Y>D': 'Charge-changed|fromAromatic|Smaller', 'Y>E': 'Charge-changed|fromAromatic|Smaller', 'Y>M': 'fromAromatic|Smaller', 'W>V': 'fromAromatic|Smaller', 'W>T': 'IncreasedPolarity|fromAromatic|Smaller', 'W>L': 'fromAromatic|Smaller', 'W>D': 'Charge-changed|fromAromatic|Smaller', 'V>T': 'IncreasedPolarity|Equal-sized', 'R>M': 'Charge-changed|Smaller', 'R>E': 'Charge-changed|Smaller', 'R>D': 'Charge-changed|Smaller', 'K>M': 'Charge-changed|Smaller', 'K>E': 'Charge-changed|Smaller', 'K>D': 'Charge-changed|Smaller', 'R>K': 'Smaller', 'A>V': 'Bigger', 'A>I': 'Bigger', 'A>L': 'Bigger', 'A>M': 'Bigger', 'A>F': 'Bigger', 'A>Y': 'Bigger|toAromatic', 'A>W': 'Bigger|toAromatic', 'A>S': 'Bigger|Flexibler|IncreasedPolarity', 'A>T': 'Bigger|IncreasedPolarity', 'A>N': 'Bigger|IncreasedPolarity', 'A>Q': 'Bigger|IncreasedPolarity', 'A>C': 'Bigger|IncreasedPolarity', 'G>A': 'Bigger|Rigider|IncreasedHydrophobicity', 'A>P': 'Bigger|toAromatic|Flexibler|IncreasedPolarity', 'A>R': 'Charge-changed|Bigger', 'A>H': 'Charge-changed|Bigger|toAromatic', 'A>K': 'Charge-changed|Bigger', 'A>D': 'Charge-changed|Bigger', 'A>E': 'Charge-changed|Bigger', 'V>V': 'Identity', 'V>I': 'Bigger', 'V>L': 'Bigger', 'V>M': 'Different-bulk', 'V>F': 'Bigger', 'V>Y': 'Bigger|toAromatic', 'V>W': 'Bigger|toAromatic', 'S>V': 'Bigger|Rigider|IncreasedHydrophobicity', 'T>V': 'Equal-sized|IncreasedHydrophobicity', 'V>N': 'IncreasedPolarity|Different-bulk', 'V>Q': 'IncreasedPolarity|Different-bulk', 'C>V': 'Bigger|IncreasedHydrophobicity', 'G>V': 'Bigger|Rigider|IncreasedHydrophobicity', 'P>V': 'Bigger|fromAromatic|Rigider|IncreasedHydrophobicity', 'V>R': 'Charge-changed|Different-bulk', 'V>H': 'Charge-changed|toAromatic|Different-bulk', 'V>K': 'Charge-changed|Different-bulk', 'V>D': 'Charge-changed|Different-bulk', 'V>E': 'Charge-changed|Different-bulk', 'I>I': 'Identity', 'I>L': 'Different-bulk', 'I>M': 'Different-bulk', 'I>F': 'Different-bulk', 'I>Y': 'toAromatic|Different-bulk', 'I>W': 'toAromatic|Different-bulk', 'S>I': 'Bigger|Rigider|IncreasedHydrophobicity', 'I>T': 'IncreasedPolarity|Different-bulk', 'I>N': 'IncreasedPolarity|Different-bulk', 'I>Q': 'IncreasedPolarity|Different-bulk', 'C>I': 'Bigger|IncreasedHydrophobicity', 'G>I': 'Bigger|Rigider|IncreasedHydrophobicity', 'P>I': 'Bigger|fromAromatic|Rigider|IncreasedHydrophobicity', 'I>R': 'Charge-changed|Different-bulk', 'I>H': 'Charge-changed|toAromatic|Different-bulk', 'I>K': 'Charge-changed|Different-bulk', 'I>D': 'Charge-changed|Different-bulk', 'I>E': 'Charge-changed|Different-bulk', 'L>I': 'Different-bulk', 'L>L': 'Identity', 'L>M': 'Different-bulk', 'L>F': 'Bigger', 'L>Y': 'Bigger|toAromatic', 'L>W': 'Bigger|toAromatic', 'S>L': 'Bigger|Rigider|IncreasedHydrophobicity', 'L>T': 'IncreasedPolarity|Different-bulk', 'N>L': 'Equal-sized|IncreasedHydrophobicity', 'L>Q': 'IncreasedPolarity|Different-bulk', 'C>L': 'Bigger|IncreasedHydrophobicity', 'G>L': 'Bigger|Rigider|IncreasedHydrophobicity', 'P>L': 'Bigger|fromAromatic|Rigider|IncreasedHydrophobicity', 'L>R': 'Charge-changed|Different-bulk', 'L>H': 'Charge-changed|toAromatic|Different-bulk', 'L>K': 'Charge-changed|Different-bulk', 'D>L': 'Charge-changed|Equal-sized', 'L>E': 'Charge-changed|Different-bulk', 'M>V': 'Different-bulk', 'M>I': 'Different-bulk', 'M>L': 'Different-bulk', 'M>M': 'Identity', 'M>F': 'Bigger', 'M>Y': 'Bigger|toAromatic', 'M>W': 'toAromatic|Different-bulk', 'S>M': 'Bigger|Rigider|IncreasedHydrophobicity', 'M>T': 'IncreasedPolarity|Different-bulk', 'M>N': 'IncreasedPolarity|Different-bulk', 'M>Q': 'Bigger|IncreasedPolarity', 'C>M': 'Bigger|IncreasedHydrophobicity', 'G>M': 'Bigger|Rigider|IncreasedHydrophobicity', 'P>M': 'Bigger|fromAromatic|Rigider|IncreasedHydrophobicity', 'M>R': 'Charge-changed|Bigger', 'M>H': 'Charge-changed|toAromatic|Different-bulk', 'M>K': 'Charge-changed|Bigger', 'M>D': 'Charge-changed|Different-bulk', 'M>E': 'Charge-changed|Bigger', 'F>I': 'Different-bulk', 'F>F': 'Identity', 'F>Y': 'Bigger|toAromatic', 'F>W': 'toAromatic|Different-bulk', 'S>F': 'Bigger|Rigider|IncreasedHydrophobicity', 'T>F': 'Bigger|IncreasedHydrophobicity', 'F>N': 'IncreasedPolarity|Different-bulk', 'F>Q': 'IncreasedPolarity|Different-bulk', 'C>F': 'Bigger|IncreasedHydrophobicity', 'G>F': 'Bigger|Rigider|IncreasedHydrophobicity', 'P>F': 'Bigger|fromAromatic|Rigider|IncreasedHydrophobicity', 'F>R': 'Charge-changed|Different-bulk', 'H>F': 'Charge-changed|Bigger|fromAromatic', 'F>K': 'Charge-changed|Different-bulk', 'D>F': 'Charge-changed|Bigger', 'E>F': 'Charge-changed|Bigger', 'Y>I': 'fromAromatic|Different-bulk', 'Y>Y': 'Identity', 'Y>W': 'Different-bulk', 'S>Y': 'Bigger|toAromatic|Rigider|IncreasedHydrophobicity', 'T>Y': 'Bigger|toAromatic|IncreasedHydrophobicity', 'Y>N': 'IncreasedPolarity|fromAromatic|Different-bulk', 'Y>Q': 'IncreasedPolarity|fromAromatic|Different-bulk', 'C>Y': 'Bigger|toAromatic|IncreasedHydrophobicity', 'G>Y': 'Bigger|toAromatic|Rigider|IncreasedHydrophobicity', 'P>Y': 'Bigger|Rigider|IncreasedHydrophobicity', 'Y>R': 'Charge-changed|fromAromatic|Different-bulk', 'H>Y': 'Charge-changed|Bigger', 'Y>K': 'Charge-changed|fromAromatic|Different-bulk', 'D>Y': 'Charge-changed|Bigger|toAromatic', 'E>Y': 'Charge-changed|Bigger|toAromatic', 'W>I': 'fromAromatic|Different-bulk', 'W>M': 'fromAromatic|Different-bulk', 'W>W': 'Identity', 'S>W': 'Bigger|toAromatic|Rigider|IncreasedHydrophobicity', 'T>W': 'Bigger|toAromatic|IncreasedHydrophobicity', 'W>N': 'IncreasedPolarity|fromAromatic|Different-bulk', 'W>Q': 'IncreasedPolarity|fromAromatic|Different-bulk', 'C>W': 'Bigger|toAromatic|IncreasedHydrophobicity', 'G>W': 'Bigger|toAromatic|Rigider|IncreasedHydrophobicity', 'P>W': 'Bigger|Rigider|IncreasedHydrophobicity', 'W>R': 'Charge-changed|fromAromatic|Different-bulk', 'H>W': 'Charge-changed|Bigger', 'W>K': 'Charge-changed|fromAromatic|Different-bulk', 'D>W': 'Charge-changed|Bigger|toAromatic', 'W>E': 'Charge-changed|fromAromatic|Different-bulk', 'S>T': 'Bigger|Rigider', 'S>N': 'Bigger|Rigider', 'S>Q': 'Bigger|Rigider', 'G>S': 'Bigger|Rigider', 'P>S': 'Bigger|fromAromatic|Rigider', 'S>R': 'Charge-changed|Bigger|Rigider', 'S>H': 'Charge-changed|Bigger|toAromatic|Rigider', 'S>K': 'Charge-changed|Bigger|Rigider', 'S>D': 'Charge-changed|Bigger|Rigider', 'S>E': 'Charge-changed|Bigger|Rigider', 'T>I': 'Different-bulk|IncreasedHydrophobicity', 'T>L': 'Different-bulk|IncreasedHydrophobicity', 'T>M': 'Different-bulk|IncreasedHydrophobicity', 'T>T': 'Identity', 'T>N': 'Different-bulk', 'T>Q': 'Different-bulk', 'C>T': 'Bigger', 'G>T': 'Bigger|Rigider', 'P>T': 'Bigger|fromAromatic|Rigider', 'T>R': 'Charge-changed|Different-bulk', 'T>H': 'Charge-changed|toAromatic|Different-bulk', 'T>K': 'Charge-changed|Different-bulk', 'T>D': 'Charge-changed|Different-bulk', 'T>E': 'Charge-changed|Different-bulk', 'N>V': 'Different-bulk|IncreasedHydrophobicity', 'N>I': 'Different-bulk|IncreasedHydrophobicity', 'N>M': 'Different-bulk|IncreasedHydrophobicity', 'N>F': 'Different-bulk|IncreasedHydrophobicity', 'N>Y': 'toAromatic|Different-bulk|IncreasedHydrophobicity', 'N>W': 'toAromatic|Different-bulk|IncreasedHydrophobicity', 'N>T': 'Different-bulk', 'N>N': 'Identity', 'N>Q': 'Bigger', 'C>N': 'Bigger', 'G>N': 'Bigger|Rigider', 'P>N': 'Bigger|fromAromatic|Rigider', 'N>R': 'Charge-changed|Different-bulk', 'N>H': 'Charge-changed|toAromatic|Different-bulk', 'N>K': 'Charge-changed|Different-bulk', 'N>D': 'Charge-changed|Equal-sized', 'N>E': 'Charge-changed|Bigger', 'Q>V': 'Different-bulk|IncreasedHydrophobicity', 'Q>I': 'Different-bulk|IncreasedHydrophobicity', 'Q>L': 'Different-bulk|IncreasedHydrophobicity', 'Q>F': 'Different-bulk|IncreasedHydrophobicity', 'Q>Y': 'toAromatic|Different-bulk|IncreasedHydrophobicity', 'Q>W': 'toAromatic|Different-bulk|IncreasedHydrophobicity', 'Q>T': 'Different-bulk', 'Q>Q': 'Identity', 'C>Q': 'Bigger', 'G>Q': 'Bigger|Rigider', 'P>Q': 'Bigger|fromAromatic|Rigider', 'Q>R': 'Charge-changed|Different-bulk', 'Q>H': 'Charge-changed|toAromatic|Different-bulk', 'Q>K': 'Charge-changed|Different-bulk', 'D>Q': 'Charge-changed|Bigger', 'Q>E': 'Charge-changed|Equal-sized', 'G>C': 'Bigger|Rigider', 'P>C': 'Bigger|fromAromatic|Rigider', 'C>R': 'Charge-changed|Bigger', 'C>H': 'Charge-changed|Bigger|toAromatic', 'C>K': 'Charge-changed|Bigger', 'C>D': 'Charge-changed|Bigger', 'C>E': 'Charge-changed|Bigger', 'G>P': 'Bigger|toAromatic|Flexibler', 'G>R': 'Charge-changed|Bigger|Rigider', 'G>H': 'Charge-changed|Bigger|toAromatic|Rigider', 'G>K': 'Charge-changed|Bigger|Rigider', 'G>D': 'Charge-changed|Bigger|Rigider', 'G>E': 'Charge-changed|Bigger|Rigider', 'P>R': 'Charge-changed|Bigger|fromAromatic|Rigider', 'P>H': 'Charge-changed|Bigger|Rigider', 'P>K': 'Charge-changed|Bigger|fromAromatic|Rigider', 'P>D': 'Charge-changed|Bigger|fromAromatic|Rigider', 'P>E': 'Charge-changed|Bigger|fromAromatic|Rigider', 'R>V': 'Charge-changed|Different-bulk', 'R>I': 'Charge-changed|Different-bulk', 'R>L': 'Charge-changed|Different-bulk', 'R>F': 'Charge-changed|Different-bulk', 'R>Y': 'Charge-changed|toAromatic|Different-bulk', 'R>W': 'Charge-changed|toAromatic|Different-bulk', 'R>T': 'Charge-changed|Different-bulk', 'R>N': 'Charge-changed|Different-bulk', 'R>Q': 'Charge-changed|Different-bulk', 'R>R': 'Identity', 'R>H': 'toAromatic|Different-bulk', 'K>R': 'Bigger', 'D>R': 'Charge-changed|Bigger', 'E>R': 'Charge-changed|Bigger', 'H>V': 'Charge-changed|fromAromatic|Different-bulk', 'H>I': 'Charge-changed|fromAromatic|Different-bulk', 'H>L': 'Charge-changed|fromAromatic|Different-bulk', 'H>M': 'Charge-changed|fromAromatic|Different-bulk', 'H>T': 'Charge-changed|fromAromatic|Different-bulk', 'H>N': 'Charge-changed|fromAromatic|Different-bulk', 'H>Q': 'Charge-changed|fromAromatic|Different-bulk', 'H>R': 'fromAromatic|Different-bulk', 'H>H': 'Identity', 'H>K': 'fromAromatic|Different-bulk', 'H>D': 'Charge-changed|fromAromatic|Different-bulk', 'H>E': 'Charge-changed|fromAromatic|Different-bulk', 'K>V': 'Charge-changed|Different-bulk', 'K>I': 'Charge-changed|Different-bulk', 'K>L': 'Charge-changed|Different-bulk', 'K>F': 'Charge-changed|Different-bulk', 'K>Y': 'Charge-changed|toAromatic|Different-bulk', 'K>W': 'Charge-changed|toAromatic|Different-bulk', 'K>T': 'Charge-changed|Different-bulk', 'K>N': 'Charge-changed|Different-bulk', 'K>Q': 'Charge-changed|Different-bulk', 'K>H': 'toAromatic|Different-bulk', 'K>K': 'Identity', 'D>K': 'Charge-changed|Bigger', 'E>K': 'Charge-changed|Bigger', 'D>V': 'Charge-changed|Different-bulk', 'D>I': 'Charge-changed|Different-bulk', 'D>M': 'Charge-changed|Different-bulk', 'D>T': 'Charge-changed|Different-bulk', 'D>H': 'Charge-changed|toAromatic|Different-bulk', 'D>D': 'Identity', 'D>E': 'Bigger', 'E>V': 'Charge-changed|Different-bulk', 'E>I': 'Charge-changed|Different-bulk', 'E>L': 'Charge-changed|Different-bulk', 'E>W': 'Charge-changed|toAromatic|Different-bulk', 'E>T': 'Charge-changed|Different-bulk', 'E>H': 'Charge-changed|toAromatic|Different-bulk', 'E>E': 'Identity'}
