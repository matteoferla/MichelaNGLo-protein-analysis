"""
This file parses the uniprot FTP file and can do various things. such as making a small one that is only human.
But mainly the `UniprotReader.convert('uniprot_sprot.xml')` method whcih generates the JSON files required. In future these will be databases...
Be warned that ET.Element is a monkeypatched version.
"""

import os, json, re
from protein.ET_monkeypatch import ET
from protein import Protein

##### Uniprot reader
class UniprotReader:
    """
    see generator iter_human
    NB. The ET.Element has been expanded. See `help(ElementalExpansion)`
    """

    def __init__(self, file):
        self.file = file

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
    def convert(cls, uniprot_master_file, first_n_protein=0):
        """

        :param uniprot_master_file:
        :param first_n_protein: set to zero for all, to interger to get the first n.
        :return:
        """
        if not os.path.isfile(uniprot_master_file):
            print('ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz')
            Protein.settings.retrieve_references(ask = True)
        count=0
        namedex = {}
        seqdex = {}
        genedex = {}
        for entry in cls(uniprot_master_file).iter_human():
            count+=1
            if count == first_n_protein:
                break
            prot = Protein(entry)
            chosen_name = prot.accession_list[0] #ideally prot.uniprot_name or the first acc id. But for code usage going for gene name.
            ## fill namedex
            namedex[prot.uniprot_name] = chosen_name
            namedex[prot.gene_name] = chosen_name
            namedex[prot.recommended_name] = chosen_name
            for group in [prot.alt_gene_name_list, prot.alternative_fullname_list, prot.alternative_shortname_list]:
                for name in group:
                    if re.match('[\d\-]+\.[\d\-]+\.[\d\-]+\.[\d\-]+',name):
                        continue # no EC numbers!
                    namedex[name] = chosen_name
            ## fill seqdex
            seqdex[chosen_name] = prot.sequence
            genedex[chosen_name] = prot.gene_name
            ## save
            prot.write_uniprot(os.path.join(Protein.settings.uniprot_folder, chosen_name+'_uniprot.xml'))
            # prot.
        ## cleanup
        for k in ('\n      ', '\n     ', '\n    ', '\n   ', '', '\n', ' '):
            if k in namedex:
                del namedex[k]
        json.dump(namedex, open(os.path.join(Protein.settings.data_folder, 'human_prot_namedex.json'), 'w'))  ## name to uniprot
        #json.dump(genedex, open(os.path.join(Protein.settings.temp_folder, 'human_prot_genedex.json'), 'w'))  ## uniprot to name. not needed.
        #json.dump(seqdex, open(os.path.join(Protein.settings.temp_folder, 'human_prot_seqdex.json'), 'w'))    ## not needed.
        open(os.path.join(Protein.settings.temp_folder, 'human.fa'), 'w').writelines(['>{i}\n{s}\n\n'.format(s=seqdex[k], i=k) for k in seqdex])
