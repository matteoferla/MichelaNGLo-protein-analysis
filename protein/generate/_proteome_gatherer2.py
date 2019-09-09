"""
This file parses the uniprot FTP file and can do various things. such as making a small one that is only human.
But mainly the `UniprotReader.convert('uniprot_sprot.xml')` method whcih generates the JSON files required. In future these will be databases...
Be warned that ET.Element is a monkeypatched version.
"""

import os, json, re
from .ET_monkeypatch import ET
from ._protein_gatherer import ProteinGatherer as Protein

from collections import defaultdict

##### Uniprot reader
class UniprotReader:
    """
    see generator iter_human
    NB. The ET.Element has been expanded. See `help(ElementalExpansion)`
    """

    def iter_human(self):
        """
        Interates across a LARGE Uniprot XML file and returns *only* the humans.
        :return: ET.Element()
        """
        for event, elem in ET.iterparse(self.file, events=('end',)):
            if elem is not None and isinstance(elem, ET.Element):
                if elem.ns_strip() == 'entry':
                    if elem.is_human():
                        yield elem
                    elem.clear()

    def iter_all(self, dataset=None):
        """
        dataset = Swiss-Prot is better than TrEMBL
        Interates across a LARGE Uniprot XML file and returns entries regardless of humanity.
        :return: ET.Element()
        """
        for event, elem in ET.iterparse(self.file, events=('end',)):
            if elem is not None and isinstance(elem, ET.Element):
                if elem.ns_strip() == 'entry':
                    if dataset and dataset != elem.get_attr('dataset'):
                        continue
                    yield elem
                    elem.clear()

    def shrink(self, outfile='human_proteome.xml'):
        """
        Make a smaller XML file, but with only the human proteome.
        :param outfile:
        :return:
        """
        with open(outfile, 'wb') as proteo:
            proteo.write(
                b'<?xml version="1.0" encoding="UTF-8"?><uniprot xmlns="http://uniprot.org/uniprot" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" ' +
                b'xsi:schemaLocation="http://uniprot.org/uniprot http://www.uniprot.org/docs/uniprot.xsd">')
            for entry in self.iter_human():
                proteo.write(ET.tostring(entry))
            proteo.write(b'</uniprot>')

    def parse_human(self):
        return [Protein(entry) for entry in self.iter_human()]

    @classmethod
    def convert(cls, uniprot_master_file=None, first_n_protein=0):
        """
        DO NOT USE!!!
        :param uniprot_master_file:
        :param first_n_protein: set to zero for all, to interger to get the first n.
        :return:
        """
        raise DeprecationWarning('DO NOT USE!')
        if uniprot_master_file is None:
            uniprot_master_file = os.path.join(Protein.settings.reference_folder, 'uniprot_sprot.xml')
            Protein.settings.retrieve_references(ask = True)
        count=0
        greater_namedex = {}
        lesser_namedex = {}
        seqdex = {}
        genedex = {}
        for entry in cls(uniprot_master_file).iter_human():
            count+=1
            if count == first_n_protein:
                break
            prot = Protein.from_uniprot(entry)
            chosen_name = prot.accession_list[0] #ideally prot.uniprot_name or the first acc id. But for code usage going for gene name.
            ## fill namedex
            greater_namedex[prot.uniprot_name] = chosen_name
            greater_namedex[prot.recommended_name] = chosen_name
            greater_namedex[prot.gene_name] = chosen_name
            for group in [prot.alt_gene_name_list, prot.alternative_fullname_list, prot.alternative_shortname_list]:
                for name in group:
                    if re.match('[\d\-]+\.[\d\-]+\.[\d\-]+\.[\d\-]+',name):
                        continue # no EC numbers!
                    lesser_namedex[name] = chosen_name
            ## fill seqdex
            seqdex[chosen_name] = prot.sequence
            genedex[chosen_name] = prot.gene_name
            ## save
            prot.write_uniprot(os.path.join(Protein.settings.uniprot_folder, chosen_name+'_uniprot.xml'))
            # prot.
        namedex = {**lesser_namedex, **greater_namedex}
        ## cleanup
        for k in ('\n      ', '\n     ', '\n    ', '\n   ', '', '\n', ' '):
            if k in namedex:
                del namedex[k]
        json.dump(namedex, open(os.path.join(Protein.settings.data_folder, 'human_prot_namedex.json'), 'w'))  ## name to uniprot
        #json.dump(genedex, open(os.path.join(Protein.settings.temp_folder, 'human_prot_genedex.json'), 'w'))  ## uniprot to name. not needed.
        #json.dump(seqdex, open(os.path.join(Protein.settings.temp_folder, 'human_prot_seqdex.json'), 'w'))    ## not needed.
        open(os.path.join(Protein.settings.temp_folder, 'human.fa'), 'w').writelines(['>{i}\n{s}\n\n'.format(s=seqdex[k], i=k) for k in seqdex])

    def __init__(self, uniprot_master_file=None, first_n_protein=0, chosen_attribute='uniprot'):
        """
        THIS IS FOR MICHELANGLO
        :param uniprot_master_file:
        :param first_n_protein: set to zero for all, to interger to get the first n.
        :return:
        """
        if uniprot_master_file is None:
            uniprot_master_file = os.path.join(Protein.settings.reference_folder, 'uniprot_sprot.xml')
            Protein.settings.retrieve_references(ask=False)
        self.file = uniprot_master_file
        count=0
        uniprot_pdbdex = defaultdict(list)
        uniprot_datasetdex = defaultdict(dict)
        organism_greater_namedex = defaultdict(dict)
        organism_lesser_namedex = defaultdict(dict)
        uniprot_namedex = {}
        uniprot_speciesdex = {}
        organismdex = {}
        seqdex = {}
        genedex = {}
        for entry in self.iter_all():
            # debugging overrider....
            count+=1
            if count == first_n_protein:
                break
            ### parser...
            prot = Protein.from_uniprot(entry)
            if prot.organism['common'] == 'Human':
                prot.parse_swissmodel().get_offsets().get_resolutions().parse_gNOMAD()
            prot.dump()  # gdump??
            ### dict
            chosen_name = getattr(prot, chosen_attribute)
            # update the organism dex
            org = prot.organism['NCBI Taxonomy']
            if not org in organismdex:
                for k in prot.organism:
                    organismdex[prot.organism[k]] = org
            ## fill namedex
            if prot.uniprot in uniprot_datasetdex:
                print('Repeated uniprot!', prot.uniprot)
                continue
            ## make dictionaries...
            uniprot_datasetdex[prot.uniprot] = prot.uniprot_dataset
            uniprot_pdbdex[prot.uniprot].extend([p.id for p in prot.pdbs])
            if prot.gene_name and prot.gene_name in organism_greater_namedex[org]:
                if prot.organism['scientific'] == 'Homo sapiens':
                    print('#'*20)
                    print('CLASH!!!!!', prot.gene_name, prot.uniprot, prot.organism['scientific'], prot.recommended_name, uniprot_datasetdex[prot.uniprot], uniprot_datasetdex[organism_greater_namedex[org][prot.gene_name]])
                if prot.uniprot_dataset == 'TrEMBL' and uniprot_datasetdex[organism_greater_namedex[org][prot.gene_name]] == 'Swiss-Prot':
                    continue
            if prot.recommended_name and prot.recommended_name in organism_greater_namedex[org]:
                if prot.organism['scientific'] == 'Homo sapiens':
                    print('#' * 20)
                    print('CLASH!!!!!', prot.gene_name, prot.uniprot, prot.recommended_name,
                      uniprot_datasetdex[prot.uniprot],
                      uniprot_datasetdex[organism_greater_namedex[org][prot.recommended_name]])
                if prot.uniprot_dataset == 'TrEMBL' and uniprot_datasetdex[organism_greater_namedex[org][prot.recommended_name]] == 'Swiss-Prot':
                    continue
            organism_greater_namedex[org][prot.uniprot_name] = chosen_name
            organism_greater_namedex[org][prot.recommended_name] = chosen_name
            organism_greater_namedex[org][prot.gene_name] = chosen_name
            uniprot_namedex[prot.uniprot] = prot.recommended_name
            uniprot_speciesdex[prot.uniprot] = org
            for group in [prot.alt_gene_name_list, prot.alternative_fullname_list, prot.alternative_shortname_list]:
                for name in group:
                    if re.match('[\d\-]+\.[\d\-]+\.[\d\-]+\.[\d\-]+',name):
                        continue # no EC numbers!
                    organism_lesser_namedex[org][name] = chosen_name
        ## final touches to the whole sets...
        for org in organism_greater_namedex:
            namedex = {**organism_lesser_namedex[org], **organism_greater_namedex[org]}
            ## cleanup
            for k in ('\n      ', '\n     ', '\n    ', '\n   ', '', '\n', ' '):
                if k in namedex:
                    del namedex[k]
            #namedex = {k.lower(): namedex[k] for k in namedex}
            fn = os.path.join(Protein.settings.dictionary_folder, f'taxid{org}-names2{chosen_attribute}.json')
            json.dump(namedex, open(fn, 'w'))  ## name to pdbs
        fn = os.path.join(Protein.settings.dictionary_folder, f'organism.json')
        json.dump(organismdex, open(fn, 'w'))  ## organism to taxid
        ## lighten
        for dex, fn in ((uniprot_pdbdex, 'uniprot2pdb.json'),
                        (uniprot_namedex, 'uniprot2name.json'),
                        (uniprot_speciesdex, 'uniprot2species.json')):
            fp = os.path.join(Protein.settings.dictionary_folder, fn)
            json.dump({k: dex[k] for k in dex if dex[k]}, open(fp, 'w'))
