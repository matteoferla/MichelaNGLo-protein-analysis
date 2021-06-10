## This code is not part of the module, it is unlinked.

"""
Several methods return the list of entities. These look like these:

Protein

[{'entity_id': 1,
  'mutation_flag': None,
  'synonym': 'Guanine nucleotide-binding protein G(o) subunit alpha',
  'weight': 25286.922,
  'sequence': 'MTLSAEERAALERSKAIEKNLKEDGISAAKDVKLLLLGADNSGKSTIVKQMKIIHGGSGGSGGTTGIVETHFTFKNLHFRLFDVGGQRSERKKWIHCFEDVTAIIFCVDLSDYNRMHESLMLFDSICNNKFFIDTSIILFLNKKDLFGEKIKKSPLTICFPEYTGPNTYEDAAAYIQAQFESKNRSPNKEIYCHMTCATDTNNAQVIFDAVTDIIIANNLRGCGLY',
  'molecule_name': ['Guanine nucleotide-binding protein G(o) subunit alpha'],
  'pdb_sequence': 'MTLSAEERAALERSKAIEKNLKEDGISAAKDVKLLLLGADNSGKSTIVKQMKIIHGGSGGSGGTTGIVETHFTFKNLHFRLFDVGGQRSERKKWIHCFEDVTAIIFCVDLSDYNRMHESLMLFDSICNNKFFIDTSIILFLNKKDLFGEKIKKSPLTICFPEYTGPNTYEDAAAYIQAQFESKNRSPNKEIYCHMTCATDTNNAQVIFDAVTDIIIANNLRGCGLY',
  'ca_p_only': False,
  'source': [{'expression_host_scientific_name': 'Spodoptera frugiperda',
    'tax_id': 9606,
    'mappings': [{'start': {'residue_number': 1},
      'end': {'residue_number': 226}}],
    'organism_scientific_name': 'Homo sapiens',
    'expression_host_tax_id': 7108}],
  'length': 226,
  'in_chains': ['A'],
  'pdb_sequence_indices_with_multiple_residues': {},
  'molecule_type': 'polypeptide(L)',
  'in_struct_asyms': ['A'],
  'sample_preparation': 'Genetically manipulated',
  'gene_name': None,
  'number_of_copies': 1}]

Ligand

[{'entity_id': 5,
  'mutation_flag': None,
  'in_struct_asyms': ['E', 'G', 'H'],
  'weight': 256.424,
  'molecule_name': ['PALMITIC ACID'],
  'chem_comp_ids': ['PLM'],
  'ca_p_only': False,
  'in_chains': ['A', 'R'],
  'molecule_type': 'bound',
  'sample_preparation': 'Synthetically obtained',
  'number_of_copies': 3},]

"""

import requests
from typing import *


class PDBMeta:
    """
    Query the PDBe for what the code parts are.
    Herein by chain the chain letter is meant, while the data of the chain is called entity... terrible.
    """

    def __init__(self, entry, chain=None):
        if chain is not None:
            self.code = entry
            self.chain = chain
        elif entry.find('_') != -1:
            self.code, self.chain = entry.split('_')
        else:
            self.code = entry
            self.chain = '?'
        # cached property
        self._data = None
        self._summary = None
        self._publication = None
        self._experiment = None
        self.boring_ligands = ('WAT', 'HOH',  # `TP3` water is ambiguous and rare
                               'LI', 'NA', 'K', 'RB',  # group 1 cations
                               'BE', 'MG', 'CA', 'SR',  # earth metal cations
                               'F', 'CL', 'BR', 'I',  # halogens
                               'MN', 'FE', 'CO', 'NI', 'CU', 'ZN',  # top period transition metals
                               '3CO',  # cobalt (iii) ion
                               'BUQ',  # 4-hydroxy-2-butanone
                               # 'NAG',  # n-acetyl-d-glucosamine
                               # 'NAD',  # nicotinamide-adenine-dinucleotide
                               'CR',  # chromium ion
                               # 'SF4',  # iron/sulfur cluster
                               'EOH',  # ethanol
                               'ZNO',  # zinc ion, 2 waters coordinated
                               'NAO',  # sodium ion, 1 water coordinated
                               'EOM',  # ethyloxymethoxyl
                               'EHN',  # ethane
                               # 'NAP',  # nadp nicotinamide-adenine-dinucleotide phosphate
                               'CCN',  # acetonitrile
                               'NAW',  # sodium ion, 3 waters coordinated
                               'BR',  # bromide ion
                               'EGL',  # ethylene glycol
                               'NI2',  # nickel (ii) ion, 2 waters coordinated
                               # 'GSH',  # glutathione
                               'NI1',  # nickel ion, 1 water coordinated
                               # 'O2',  # oxygen molecule
                               'BA',  # barium ion
                               'RU',  # ruthenium ion
                               # 'SAH',  # s-adenosyl-l-homocysteine
                               'GU7',
                               # 2-amino-7-[2-(2-hydroxy-1-hydroxymethyl-ethylamino)-ethyl]-1,7-dihydro-purin-6-one
                               # 'SAM',  # s-adenosylmethionine
                               'TAS',  # trihydroxyarsenite(iii)
                               'DCE',  # 1,2-dichloroethane
                               '2BM',  # dibromomethane
                               # 'TP7',  # coenzyme b
                               'OF3',  # ferric ion, 1 water coordinated
                               'OF1',  # ferrous ion, 1 water coordinated
                               'RB',  # rubidium ion
                               'IOH',  # 2-propanol, isopropanol
                               'MW1',  # manganese ion, 1 water coordinated
                               'IOD',  # iodide ion
                               'C2O',  # cu-o-cu linkage
                               'BNZ',  # benzene
                               'TCN',  # tetracyanonickelate ion
                               'ARS',  # arsenic
                               'NH4',  # ammonium ion
                               'GD',  # gadolinium atom
                               # 'PER',  # peroxide ion
                               'GA',  # gallium (iii) ion
                               # 'TPP',  # thiamine diphosphate
                               'CHX',  # cyclohexane
                               # 'CME',  # s,s-(2-hydroxyethyl)thiocysteine
                               # 'THB',  # tetrahydrobiopterin
                               'IPA',  # isopropyl alcohol
                               'CD1',  # cadmium ion, 1 water coordinated
                               'OH',  # hydroxide ion
                               'SO4',  # sulfate ion
                               'DTT',  # 2,3-dihydroxy-1,4-dithiobutane
                               # 'PQN',  # phylloquinone
                               'CYN',  # cyanide ion
                               # 'PQQ',  # pyrroloquinoline quinone
                               'PYJ',  # phenylethane
                               # 'PEO',  # hydrogen peroxide
                               'NA6',  # sodium ion, 6 waters coordinated
                               'MBR',  # tribromomethane
                               'NA5',  # sodium ion, 5 waters coordinated
                               'OS',  # osmium ion
                               'MAN',  # alpha-d-mannose
                               'CMO',  # carbon monoxide
                               'OCL',  # cobalt ion, 1 water coordinated
                               'DMF',  # dimethylformamide
                               'OCN',  # cobalt ion, 2 waters coordinated
                               'MO3',  # magnesium ion, 3 waters coordinated
                               'NGN',  # nitrogen
                               'ACT',  # acetate ion
                               'U1',  # uranium atom
                               'HDZ',  # nitrogen molecule
                               'MO5',  # magnesium ion, 5 waters coordinated
                               'MO4',  # magnesium ion, 4 waters coordinated
                               'VO4',  # vanadate ion
                               'DMS',  # dimethyl sulfoxide
                               'FUC',  # alpha-l-fucose
                               'PCL',  # platinum(ii) di-chloride
                               'CB5',  # cobalt bis(1,2-dicarbollide)
                               'EEE',  # ethyl acetate
                               'HG',  # mercury (ii) ion
                               'NO2',  # nitrite ion
                               # 'CMP',  # adenosine-3',5'-cyclic-monophosphate
                               'PR',  # praseodymium ion
                               'BMA',  # beta-d-mannose
                               'IUM',  # uranyl (vi) ion
                               'PT',  # platinum (ii) ion
                               'ZN2',  # zinc ion on 3-fold crystal axis
                               # 'TTP',  # thymidine-5'-triphosphate
                               'NO3',  # nitrate ion
                               'YT3',  # yttrium (iii) ion
                               # 'TYS',  # o-sulfo-l-tyrosine
                               'PB',  # lead (ii) ion
                               'M2M',  # 1-methoxy-2-(2-methoxyethoxy)ethane
                               'ZO3',  # zinc ion, 3 waters coordinated
                               'PD',  # palladium ion
                               # 'AMP',  # adenosine monophosphate
                               'PI',  # hydrogenphosphate ion
                               'MH3',  # manganese ion, 1 hydroxyl coordinated
                               'AF3',  # aluminum fluoride
                               'ZN',  # zinc ion
                               'MN3',  # manganese (iii) ion
                               'OXY',  # oxygen molecule
                               'NI',  # nickel (ii) ion
                               # 'CSD',  # 3-sulfinoalanine
                               # 'OX',  # bound oxygen
                               'PS5',  # pentasulfide-sulfur
                               'MN5',  # manganese ion, 5 waters coordinated
                               'MN6',  # manganese ion, 6 waters coordinated
                               'S',  # sulfur atom
                               'HOH',  # water
                               'W',  # tungsten ion
                               'SB',  # antimony (iii) ion
                               # 'FOL',  # folic acid
                               'OXE',  # ortho-xylene
                               'PT4',  # platinum (iv) ion
                               'PBM',  # trimethyl lead ion
                               'O',  # oxygen atom
                               'MW2',  # manganese dihydrate ion
                               'MG',  # magnesium ion
                               '543',  # calcium ion, 6 waters plus ethanol coordinated
                               'MSM',  # (methylsulfanyl)methane
                               # 'C5P',  # cytidine-5'-monophosphate
                               'ANL',  # aniline
                               'MTO',  # bound water
                               'NO',  # nitric oxide
                               'TBU',  # tertiary-butyl alcohol
                               'OPY',  # (3s)-4-oxo-4-piperidin-1-ylbutane-1,3-diamine
                               'PC4',  # tetrachloroplatinate(ii)
                               # 'GU3',  # methyl 3-o-methyl-2,6-di-o-sulfo-alpha-d-glucopyranoside
                               # 'GU2',  # 2,3-di-o-methyl-alpha-l-idopyranuronic acid
                               # 'GU1',  # 2,3-di-o-methyl-beta-d-glucopyranuronic acid
                               'MOH',  # methanol
                               # 'ANP',  # phosphoaminophosphonic acid-adenylate ester
                               # 'GU6',  # 2,3,6-tri-o-sulfonato-alpha-d-glucopyranose
                               # 'GU5',  # 2,3-di-o-methyl-6-o-sulfonato-alpha-d-glucopyranose
                               # 'GU4',  # 2,3,4,6-tetra-o-sulfonato-alpha-d-glucopyranose
                               'AU',  # gold ion
                               'OC3',  # calcium ion, 3 waters coordinated
                               'BTN',  # biotin
                               'I42',  # hydroxy(dioxido)oxovanadium
                               'OC4',  # calcium ion, 4 waters coordinated
                               'OC7',  # calcium ion, 7 waters coordinated
                               'OC6',  # calcium ion, 6 waters coordinated
                               # 'TMP',  # thymidine-5'-phosphate
                               'RE',  # rhenium
                               'GD3',  # gadolinium ion
                               # 'CTP',  # cytidine-5'-triphosphate
                               'ACE',  # acetyl group
                               '3OF',  # hydrated fe (iii) ion, 2 waters coordinated
                               'ETZ',  # diethyl ether
                               'MM4',  # molybdenum (iv) oxide
                               'IN',  # indium (iii) ion
                               'ACN',  # acetone
                               'DOD',  # deuterated water
                               'AST',  # arsenite
                               # 'COA',  # coenzyme a
                               'EU',  # europium ion
                               'DOX',  # dioxane
                               # 'COB',  # co-methylcobalamin
                               # 'B12',  # cobalamin
                               'REO',  # perrhenate
                               # 'ATP',  # adenosine-5'-triphosphate
                               'CD3',  # cadmium ion, 3 waters coordinated
                               # 'U10',  # ubiquinone-10
                               'ACY',  # acetic acid
                               'PEG',  # di(hydroxyethyl)ether
                               'YB',  # ytterbium (iii) ion
                               # 'NDP',  # nadph dihydro-nicotinamide-adenine-dinucleotide phosphate
                               'NBZ',  # nitrobenzene
                               'ETI',  # iodoethane
                               'SER',  # serine
                               'C2C',  # cu-cl-cu linkage
                               'NA',  # sodium ion
                               'FMT',  # formic acid
                               'ASC',  # ascorbic acid
                               'AU3',  # gold 3+ ion
                               'FE2',  # fe (ii) ion
                               'LNK',  # pentane
                               'SEK',  # selenocyanate ion
                               'MO1',  # magnesium ion, 1 water coordinated
                               'EU3',  # europium (iii) ion
                               '1BO',  # 1-butanol
                               'AUC',  # gold (i) cyanide ion
                               'CLO',  # chloro group
                               'FE',  # fe (iii) ion
                               'DUM',  # dummy atoms
                               # 'ADP',  # adenosine-5'-diphosphate
                               'OF2',  # 2 ferric ion, 1 bridging oxygen
                               'BEF',  # beryllium trifluoride ion
                               'FEL',  # hydrated fe
                               'BF4',  # beryllium tetrafluoride ion
                               'HEX',  # hexane
                               'CUZ',  # (mu-4-sulfido)-tetra-nuclear copper ion
                               # 'NDG',  # 2-(acetylamino)-2-deoxy-a-d-glucopyranose
                               'XE',  # xenon
                               # 'FMN',  # flavin mononucleotide
                               'YAN',  # 1,2-dichlorobenzene
                               'CUA',  # dinuclear copper ion
                               'V',  # vanadium ion
                               'CUO',  # cu2-o2 cluster
                               # 'HEM',  # protoporphyrin ix containing fe
                               # 'GMP',  # guanosine
                               'CU',  # copper (ii) ion
                               'MGF',  # trifluoromagnesate
                               # 'GDP',  # guanosine-5'-diphosphate
                               'CFT',  # trifluoromethane
                               'SBT',  # 2-butanol
                               # 'PLP',  # pyridoxal-5'-phosphate
                               'SR',  # strontium ion
                               'FU1',  # tetrahydrofuran
                               'EDN',  # ethane-1,2-diamine
                               'EDO',  # 1,2-ethanediol
                               'H2S',  # hydrosulfuric acid
                               'ND4',  # ammonium cation with d
                               'BRO',  # bromo group
                               'KR',  # krypton
                               'CS',  # cesium ion
                               'NME',  # methylamine
                               # 'CDP',  # cytidine-5'-diphosphate
                               'HGI',  # mercury (ii) iodide
                               'SM',  # samarium (iii) ion
                               # 'ALY',  # n(6)-acetyllysine
                               # 'NMO',  # nitrogen monoxide
                               # 'TDP',  # thiamin diphosphate
                               'SE',  # selenium atom
                               'HO',  # holmium atom
                               '3CN',  # 3-aminopropane
                               'AZI',  # azide ion
                               # 'F42',  # coenzyme f420
                               'FLO',  # fluoro group
                               '6MO',  # molybdenum(vi) ion
                               'EMC',  # ethyl mercury ion
                               'Y1',  # yttrium ion
                               # 'MO7', # bis(mu4-oxo)-bis(mu3-oxo)-octakis(mu2-oxo)-dodecaoxo-heptamolybdenum (vi)
                               'SE4',  # selenate ion
                               'BF2',  # beryllium difluoride
                               'CO',  # cobalt (ii) ion
                               # 'NGD', # 3-(aminocarbonyl)-1-[(2r,3r,4s,5r)-5-({[(s)-{[(s)-{[(2r,3s,4r,5r)-5-(2-amino-6-oxo-1,6-dihydro-9h-purin-9-yl)-3,4-dihydroxytetrahydrofuran-2-yl]methoxy}(hydroxy)phosphoryl]oxy}(hydroxy)phosphoryl]oxy}methyl)-3,4-dihydroxytetrahydrofuran-2-yl]pyridinium
                               '2MO',  # molybdenum (iv)oxide
                               '202',  # bromic acid
                               'DIS',  # disordered solvent
                               'MBN',  # toluene
                               'LA',  # lanthanum (iii) ion
                               'PGO',  # s-1,2-propanediol
                               'CL',  # chloride ion
                               'HP6',  # heptane
                               'SO2',  # sulfur dioxide
                               'LI',  # lithium ion
                               # 'PPS',  # 3'-phosphate-adenosine-5'-phosphate sulfate
                               # 'TPO',  # phosphothreonine
                               'POL',  # n-propanol
                               # 'GU0',  # 2,3,6-tri-o-sulfonato-alpha-l-galactopyranose
                               'SGM',  # monothioglycerol
                               'DTU',  # (2r,3s)-1,4-dimercaptobutane-2,3-diol
                               'MOO',  # molybdate ion
                               'TE',  # tellurium
                               'TB',  # terbium(iii) ion
                               'CA',  # calcium ion
                               # 'FAD',  # flavin-adenine dinucleotide
                               'CNV',  # propanenitrile
                               'GOL',  # glycerol
                               'SCN',  # thiocyanate ion
                               'AG',  # silver ion
                               'PO4',  # phosphate ion
                               'IR',  # iridium ion
                               'DIO',  # 1,4-diethylene dioxide
                               'NH2',  # amino group
                               '8CL',  # chlorobenzene
                               '3NI',  # nickel (iii) ion
                               'IRI',  # iridium hexammine ion
                               # 'UTP',  # uridine 5'-triphosphate
                               'AR',  # argon
                               # 'N4M', # 5-formyltetrahydromethanopterin
                               'CE',  # cerium (iii) ion
                               'NH3',  # ammonia
                               'MN',  # manganese (ii) ion
                               'CNN',  # cyanamide
                               'HGC',  # methyl mercury ion
                               # 'GU8',  # 2,3,6-tri-o-methyl-beta-d-glucopyranose
                               # 'GTP',  # guanosine-5'-triphosphate
                               # 'UDP',  # uridine-5'-diphosphate
                               'OC2',  # calcium ion, 2 waters coordinated
                               'ART',  # arsenate
                               'TFH',  # nitrogen of trifluoro-ethylhydrazine
                               'MCH',  # trichloromethane
                               '2NO',  # nitrogen dioxide
                               '6WO',  # oxo-tungsten(vi)
                               'CD5',  # cadmium ion, 5 waters coordinated
                               # 'KCX',  # lysine nz-carboxylic acid
                               'E1H',  # ethanimine
                               'ARF',  # formamide
                               'TL',  # thallium (i) ion
                               'DXE',  # 1,2-dimethoxyethane
                               # 'GU9',  # 2,3,6-tri-o-methyl-alpha-d-glucopyranose
                               'IDO',  # iodo group
                               'KO4',  # potassium ion, 4 waters coordinated
                               'NRU',  # ruthenium (iii) hexaamine ion
                               '4MO'  # molybdenum(iv) ion
                               )

    # ---- cached props (cannot use 3.8 method)

    def _cached(self, attr_name, url, first=True):
        if getattr(self, attr_name) is not None:
            return getattr(self, attr_name)
        else:
            reply = requests.get(url + self.code).json()
            reply_inner = reply[self.code.lower()]
            if not first:
                setattr(self, attr_name, reply_inner)
            elif len(reply_inner) == 0:
                setattr(self, attr_name, None)
            else:
                setattr(self, attr_name, reply[self.code.lower()][0])
            return getattr(self, attr_name)

    @property
    def data(self):
        return self._cached('_data', 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/', False)

    @property
    def summary(self):
        # keys= ['related_structures', 'experimental_method', 'assemblies', 'title', 'release_date', 'split_entry',
        # 'experimental_method_class', 'revision_date', 'entry_authors', 'deposition_site', 'number_of_entities',
        # 'deposition_date', 'processing_site']
        return self._cached('_summary', 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/', True)

    @property
    def publication(self):
        # keys= ['related_structures', 'experimental_method', 'assemblies', 'title', 'release_date', 'split_entry',
        # 'experimental_method_class', 'revision_date', 'entry_authors', 'deposition_site', 'number_of_entities',
        # 'deposition_date', 'processing_site']
        return self._cached('_publication', 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/publications/', True)

    @property
    def experiment(self):
        return self._cached('_experiment', 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/experiment/', True)

    # ---------------------

    def get_data_by_chain(self, chain=None):
        if not chain:
            chain = self.chain
        return [entity for entity in self.data if chain in entity['in_chains']][0]

    def get_range_by_chain(self, chain=None):
        if not chain:
            chain = self.chain
        entity = self.get_data_by_chain(chain)
        return self.get_range_by_entity(entity)

    def get_range_by_entity(self, entity):
        if 'source' in entity:
            if len(entity['source']):
                mappings = entity['source'][0]['mappings']
                if len(entity['source'][0]['mappings']) > 1:
                    raise ValueError('MULTIPLE MAPPINGS?!')
                s = mappings[0]['start']['residue_number']
                e = mappings[0]['end']['residue_number']
                return (s, e)
            else:  # synthetic.
                return (1, len(entity['sequence']))
        else:
            raise ValueError('This is not a peptide')

    def get_proteins(self):
        return [entity for entity in self.data if self.is_peptide(entity)]

    def get_polymers(self):
        return [entity for entity in self.data if self.is_polymer(entity)]

    def get_nonproteins(self):
        return [entity for entity in self.data if not self.is_peptide(entity)]

    def is_peptide(self, entity):
        return entity['molecule_type'] == 'polypeptide(L)'

    def is_DNA(self, entity):
        return entity['molecule_type'] == 'polydeoxyribonucleotide'

    def is_RNA(self, entity):
        return entity['molecule_type'] == 'polyribonucleotide'

    def is_polymer(self, entity):
        return self.is_peptide(entity) or self.is_DNA(entity) or self.is_RNA(entity)

    def wordy_describe_entity(self, entity):
        if self.is_peptide(entity):
            return f'{"/".join(entity["molecule_name"])} as chain {"/".join(entity["in_chains"])} [{"-".join([str(n) for n in self.get_range_by_entity(entity)])}]'
        else:
            return f'{"/".join(entity["molecule_name"])} in chain {"/".join(entity["in_chains"])}'

    def is_boring_ligand(self, entity):
        # DNA and RNA do not have chem_comp_ids key so will be classed as "boring"
        if 'chem_comp_ids' not in entity:
            return True  # this entity isnt even a thing
        elif len(entity['chem_comp_ids']) == 0:
            return True  # this entity isnt even a thing
        else:
            return entity['chem_comp_ids'][0] in self.boring_ligands

    def wordy_describe(self, delimiter=' + '):
        descr = delimiter.join([self.wordy_describe_entity(entity) for entity in self.get_proteins()])
        descr += ' &mdash; ' + delimiter.join(
            [self.wordy_describe_entity(entity) for entity in self.get_nonproteins() if
             not self.is_boring_ligand(entity)])
        return f'<span class="prolink" name="pdb" data-code="{self.code}" data-chain="{self.chain}">{self.code}</span> ({descr})'

    def describe(self):
        peptide = [(f'{"-".join([str(n) for n in self.get_range_by_entity(entity)])}:{chain}',
                    "/".join(entity["molecule_name"])) for entity in self.get_proteins() for chain in
                   entity["in_chains"]]
        hetero = [(f'{entity["chem_comp_ids"][0]} and :{chain}',
                   "/".join(entity["molecule_name"])) for entity in self.get_nonproteins() for chain in
                  entity["in_chains"] if not self.is_boring_ligand(entity)]
        return {'peptide': peptide, 'hetero': hetero}

    def get_other_proteins(self, chain: str) -> List:
        """
        Get the peptide chains that are not this.
        """
        return [entity for entity in self.get_proteins() if chain not in entity['in_chains']]

    def get_other_polymers(self, chain: str) -> List:
        """
        Get the polymer chains that are not this.
        """
        return [entity for entity in self.get_polymers() if chain not in entity['in_chains']]

    def is_antibody(self, entity: dict) -> bool:
        if 'molecule_name' not in entity:
            return False  # it is not an antibody... but is as un-useful as one ?!
        elif any([term in entity['molecule_name'] for term in ['antibody', 'anti-body', 'nanobody']]):
            return True
        else:
            return False

    def get_polymer_chains(self, first_only: bool = False) -> Set:
        if first_only:
            return {entity['in_chains'][0] for entity in self.get_polymers() if not self.is_antibody(entity)}
        else:
            return {e_chain for entity in self.get_polymers() for e_chain in entity['in_chains']
                    if not self.is_antibody(entity)}

    def get_other_chains(self, chain: str, first_only: bool = False) -> Set:
        if first_only:
            return {entity['in_chains'][0] for entity in self.get_other_polymers(chain) if not self.is_antibody(entity)}
        else:
            return {e_chain for entity in self.get_other_polymers(chain) for e_chain in entity['in_chains']
                    if not self.is_antibody(entity)}

    def get_chains(self, first_only: bool = False) -> Set:
        if first_only:
            return {entity['in_chains'][0] for entity in self.get_polymers() if not self.is_antibody(entity)}
        else:
            return {e_chain for entity in self.get_polymers() for e_chain in entity['in_chains']
                    if not self.is_antibody(entity)}

    def get_interesting_ligands(self) -> List:
        return [entity for entity in self.get_nonproteins() if not self.is_boring_ligand(entity)]

    def get_interesting_ligand_names(self) -> Set:
        return {entity['chem_comp_ids'][0] for entity in self.get_nonproteins() if not self.is_boring_ligand(entity)}

    def remove_chain(self, chain):
        for entity in list(self.data):
            if 'in_chains' not in entity:
                continue  # What is this though??
            elif chain in entity['in_chains'] and len(entity['in_chains']) == 1:
                self.data.remove(entity)
            elif chain in entity['in_chains']:
                entity['in_chains'].remove(chain)
            else:
                pass  # no match

    def move_chain(self, old_chain, new_chain):
        for entity in list(self.data):
            if 'in_chains' not in entity:
                continue  # What is this though??
            elif old_chain in entity['in_chains'][0]:
                entity['in_chains'].remove(old_chain)
                entity['in_chains'].append(new_chain)
            else:
                pass  # no match

    @property
    def experimental_method(self):
        return self.experiment['experimental_method']

    @property
    def resolution(self):
        if 'resolution' in self.experiment:
            return self.experiment['resolution']
        else:
            return None

    @classmethod
    def bulk_resolution(cls, codes: list):
        reply = requests.post('https://www.ebi.ac.uk/pdbe/api/pdb/entry/experiment/', data=','.join(codes)).json()
        return {code: info[0]['resolution'] if 'resolution' in info[0] else 0. for code, info in reply.items()}
