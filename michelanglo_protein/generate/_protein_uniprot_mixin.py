############### UNIPROT PARSING METHODS ###############
from michelanglo_protein.generate._protein_base_mixin import _BaseMixin
# this class requires mixin classes from `_protein_*_mixin.py`. The XML parser for Uniprot requires a special ET from `.ET_monkeypatched`.
from ..core import Structure

_failsafe = _BaseMixin._failsafe

class _UniprotMixin:

    @_failsafe
    def _parse_protein_element(self, elem):
        for name_el in elem:
            if name_el.is_tag('recommendedName') and name_el.has_text():
                self.recommended_name = name_el.text.rstrip()
            elif name_el.is_tag('recommendedName'):
                for subname_el in name_el:
                    if subname_el.is_tag('fullName') and subname_el.text:
                        self.recommended_name = subname_el.text.rstrip()
                    elif subname_el.text:
                        self.alternative_shortname_list.append(subname_el.text.rstrip())
            else:
                for subname_el in name_el:
                    if subname_el.is_tag('fullName') and subname_el.text:
                        self.alternative_fullname_list.append(subname_el.text.rstrip())
                    elif subname_el.text:
                        self.alternative_shortname_list.append(subname_el.text.rstrip())

    @_failsafe
    def _parse_protein_dbReference(self, elem):
        if elem.has_attr('type','Pfam'):
            pass
            # no need as pfam is parsed separately as it has better info.
            #self.Uniprot_pfam = [(xref['id'], xref['property'][0]['value']) for xref in clean_dict[''] if xref['type'] == 'Pfam']
        elif elem.has_attr('type','GO'):
            pass
        elif elem.has_attr('type','PDB'):
            chain = elem.get_sub_by_type('chains')
            if chain is not None:  # this is so unpredictable. It needs to be done by blast.
                loca=chain.attrib['value'].split('=')[1].split('-')
                if loca[0] != '':
                    chainid=chain.attrib['value'].split('=')[0].split('/')[0]
                    model = Structure(description=elem.attrib['id'],
                                       id=elem.attrib['id']+'_'+chainid,
                                       chain=chainid,
                                       code=elem.attrib['id'],
                                       x=int(loca[0].split(',')[0]),
                                       y=int(loca[1].split(',')[0]),
                                       type='rcsb')
                    self.pdbs.append(model)
        elif elem.has_attr('type', 'Ensembl'):
            self.ENST = elem.attrib['id']
            for subelem in elem:
                if subelem.is_tag('molecule'):
                    if subelem.attrib['id'][-2:] != '-1':
                        return None # this is not the isoform 1 !!!
                elif subelem.has_attr('type', 'protein sequence ID'):
                    self.ENSP = subelem.attrib['value']
                elif subelem.has_attr('type', 'gene ID'):
                    self.ENSG = subelem.attrib['value']
        else:
            pass

    @_failsafe
    def _parse_protein_gene(self, elem):
        for name_el in elem:
            if name_el.is_tag('name'):
                if name_el.has_attr('type', 'primary'):
                    self.gene_name = name_el.text
                elif name_el.has_attr('type', 'synonym'):
                    self.alt_gene_name_list.append(name_el.text)
        return self

    @_failsafe
    def _parse_protein_comment(self,elem):
        if elem.has_attr('type','interaction') or elem.has_attr('type','subunit'):
            for subelem in elem:
                if subelem.is_tag('interactant'):
                    partner = subelem.get_subtag('label')
                    if partner is not None:
                        self.partners['interactant'].append(partner.text)
                elif subelem.is_tag('text'):  # some entries are badly annotated and have only a text line..
                    self.partners['interactant'].append(subelem.text)
        elif elem.has_attr('type','disease'):
            for subelem in elem:
                if subelem.is_tag('disease'): #yes. the comment type=disease has a tag disease. wtf
                    disease = {'id': subelem.attrib['id']}
                    for key in ('description', 'name'):
                        subsub = subelem.get_subtag(key)
                        if subsub is not None:
                            disease[key] = subsub.text
                    mim = subelem.get_sub_by_type('MIM')
                    if min:
                        disease['MIM'] = mim.attrib['id']
                    self.diseases.append(disease)
        elif elem.has_attr('type','alternative products'):
            pass
            # The XML comment type does not contain the data, feature does kind of.
            #if 'isoform' not in self.features:
            #    self.features['isoform'] = []
            #for subelem in elem:
            #    if subelem.is_tag('isoform'):
            #        self.features['isoform'].append(self._get_location(subelem))

    @_failsafe
    def _parse_protein_feature(self,elem):
        """ These are the possible feature types:
        * active site
        * binding site
        * calcium-binding region
        * chain
        * coiled-coil region
        * compositionally biased region
        * cross-link
        * disulfide bond
        * DNA-binding region
        * domain
        * glycosylation site
        * helix
        * initiator methionine
        * lipid moiety-binding region
        * metal ion-binding site
        * modified residue
        * mutagenesis site
        * non-consecutive residues
        * non-terminal residue
        * nucleotide phosphate-binding region
        * peptide
        * propeptide
        * region of interest
        * repeat
        * non-standard amino acid
        * sequence conflict
        * sequence variant
        * short sequence motif
        * signal peptide
        * site
        * splice variant
        * strand
        * topological domain
        * transit peptide
        * transmembrane region
        * turn
        * unsure residue
        * zinc finger region
        * intramembrane region"""
        if elem.attrib['type'] not in self.features:  #avoiding defaultdictionary to avoid JSON issue.
            self.features[elem.attrib['type']]=[]
        locadex=self._get_location(elem)
        if locadex:
            self.features[elem.attrib['type']].append(self._get_location(elem))
        elif elem.has_attr('type', 'chain'):
            pass
        else:
            #print('no location?') TODO fix this!
            #print(elem.attrib)
            pass
        return self

    @_failsafe
    def _parse_organism(self, elem):
        for organism_el in list(elem):
            for mode in ("scientific", "common"):
                if 'type' in organism_el.attrib and organism_el.attrib['type'] == mode and organism_el.text:
                    self.organism[mode] = organism_el.text.rstrip().lstrip()
                    break
            else:
                if 'type' in organism_el.attrib and organism_el.attrib['type'] == "NCBI Taxonomy":
                    self.organism["NCBI Taxonomy"] = organism_el.attrib['id']
                elif organism_el.text and organism_el.text.rstrip().lstrip():
                    self.organism["other"] = organism_el.text.rstrip().lstrip()
                else:
                    pass #lineage.

    def _get_location(self, elem):
        location = elem.get_subtag('location')
        if location is None:
            return None
        position = location.get_subtag('position')
        start = location.get_subtag('start')
        if start is None:
            start = location.get_subtag('begin')
        end =  location.get_subtag('end')
        if elem.has_attr('description'):
            description = elem.attrib['description']
        else:
            description = '-'
        if position is not None:  # single residue
            x = position.attrib['position']
            return {'x':int(x), 'y': int(x), 'description': description, 'id': '{t}_{s}'.format(s=x, t=elem.attrib['type'].replace(' ','').replace('-',''))}
        elif start is not None and end is not None and not start.has_attr('status') and not end.has_attr('status'):
            x = start.attrib['position']
            y = end.attrib['position']
            return {'x': int(x), 'y': int(y), 'description': description, 'id': '{t}_{x}_{y}'.format(x=x, y=y, t=elem.attrib['type'].replace(' ', '').replace('-',''))}
        elif end is not None and not end.has_attr('status'):
            x = 1
            y = end.attrib['position']
            return {'x': int(x), 'y': int(y), 'description': description, 'id': '{t}_{x}_{y}'.format(x=x, y=y, t=elem.attrib['type'].replace(' ', '').replace('-',''))}
        elif elem.has_attr('type', 'chain'):
            return None #pointless chain.
        else:
            #print('Unexpected location entry') TODO fix this!
            #print(elem.attrib, position, start, end)
            return None

    def _parse_uniprot_xml(self, entry):
        """
        loads the Protein instance with the data from the uniprot XML entry element.
        Unlike the previous version the elemtn tree is parsed directly as opposed to converitng into a dictionary that seemed at first a wiser strategy but wasn't.
        Do note htat the ET.Element has to be monkeypatched. See `help(ElementalExpansion)`
        :param entry: the entry element of the XML parsed by Element Tree.
        :return:
        """
        self.features = {}
        self.uniprot_dataset = entry.get_attr('dataset')
        for elem in entry:
            if elem.is_tag('accession') and elem.text:
                self.accession_list.append(elem.text.rstrip())
            elif elem.is_tag('organism'):
                self._parse_organism(elem)
            elif elem.is_tag('name') and elem.text:
                self.uniprot_name = elem.text.rstrip()
            elif elem.is_tag('sequence') and elem.text:
                self.sequence = elem.text.rstrip().replace(' ','')
            elif elem.is_tag('protein'):
                self._parse_protein_element(elem)
            elif elem.is_tag('dbReference'):
                self._parse_protein_dbReference(elem)
            elif elem.is_tag('gene'):
                self._parse_protein_gene(elem)
            elif elem.is_tag('feature'):
                self._parse_protein_feature(elem)
            elif elem.is_tag('comment'):
                self._parse_protein_comment(elem)
            elif elem.is_tag('keyword'):
                pass
        # Ater all the parsing some weirdnesses can arise
        for group in self.features:
            for i in range(len(self.features[group])):
                self.features[group][i]['type'] = group  # in case of flattening.
                if 'description' in self.features[group][i] and self.features[group][i]['description'] == '-':
                    self.features[group][i]['description'] = group
        if self.uniprot == '':
            if len(self.accession_list):
                self.uniprot = self.accession_list[0]
            else:
                self.uniprot = self.uniprot_name
        return self
