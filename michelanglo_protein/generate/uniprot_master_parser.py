"""
This file parses the uniprot FTP file and can do various things. such as making a small one that is only human.
But mainly the `UniprotMasterReader.convert('uniprot_sprot.xml')` method whcih generates the JSON files required. In future these will be databases...
Be warned that ET.Element is a monkeypatched version.
"""

import os, json, re, time
from threading import Thread, Semaphore, Lock, active_count
import itertools
from .ET_monkeypatch import ET
from ._protein_gatherer import ProteinGatherer as Protein
from warnings import warn

from michelanglo_protein.generate.split_gnomAD import gnomAD

from collections import defaultdict

### Uniprot reader
class UniprotMasterReader:
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
        count = 0
        for event, elem in ET.iterparse(self.file, events=('end',)):
            if elem is not None and isinstance(elem, ET.Element):
                if elem.ns_strip() == 'entry':
                    if dataset and dataset != elem.get_attr('dataset'):
                        continue
                    count += 1
                    if count == self.first_n_protein:
                        break
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
            # fill namedex
            greater_namedex[prot.uniprot_name] = chosen_name
            greater_namedex[prot.recommended_name] = chosen_name
            greater_namedex[prot.gene_name] = chosen_name
            for group in [prot.alt_gene_name_list, prot.alternative_fullname_list, prot.alternative_shortname_list]:
                for name in group:
                    if re.match('[\d\-]+\.[\d\-]+\.[\d\-]+\.[\d\-]+',name):
                        continue # no EC numbers!
                    lesser_namedex[name] = chosen_name
            # fill seqdex
            seqdex[chosen_name] = prot.sequence
            genedex[chosen_name] = prot.gene_name
            # save
            prot.write_uniprot(os.path.join(Protein.settings.uniprot_folder, chosen_name+'_uniprot.xml'))
            # prot.
        namedex = {**lesser_namedex, **greater_namedex}
        # cleanup
        for k in ('\n      ', '\n     ', '\n    ', '\n   ', '', '\n', ' '):
            if k in namedex:
                del namedex[k]
        json.dump(namedex, open(os.path.join(Protein.settings.data_folder, 'human_prot_namedex.json'), 'w'))  # name to uniprot
        #json.dump(genedex, open(os.path.join(Protein.settings.temp_folder, 'human_prot_genedex.json'), 'w'))  # uniprot to name. not needed.
        #json.dump(seqdex, open(os.path.join(Protein.settings.temp_folder, 'human_prot_seqdex.json'), 'w'))    # not needed.
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
            #Protein.settings.retrieve_references(ask=False)
            if not os.path.exists(os.path.join(Protein.settings.data_folder, 'gnomAD')):
                #gnomAD().split().write('gnomAD')
                warn('I turned off this step!!!')
        self.file = uniprot_master_file
        self.first_n_protein = first_n_protein
        self.chosen_attribute = chosen_attribute
        self._uniprot_pdbdex = defaultdict(list)  #: dict of key=uniprot and value=list of PDB codes
        self._uniprot_datasetdex = defaultdict(str)  #: dict of key=uniprot and value=dataset type Swiss-Prot | TrEMBL
        self._organism_greater_namedex = defaultdict(dict)  #: dict of taxid as keys, with dict value with keys names and values prefered name (uniprot id) for the good names
        self._organism_lesser_namedex = defaultdict(dict)  #: dict of taxid as keys, with dict value with keys names and values prefered name (uniprot id) for the potentially stinky names
        self._uniprot_namedex = {}  #: dict of key=uniprot and value=human name for title
        self._uniprot_speciesdex = {}  #: dict of key=uniprot and value=taxid
        self._organismdex = {}  #: dict of key=orgnism names and value=taxid
        self._semaphore = Semaphore(50)
        self._lock = Lock()
        idleness = active_count()
        # preload some steps to speed up
        with Protein.settings.open('resolution') as fh:
            resolutions = {entry['IDCODE']: float(entry['RESOLUTION']) for entry in json.load(fh) if
                           entry['RESOLUTION'].strip()}
        self.resolutions = resolutions
        # run
        for entry in self.iter_all():
            Thread(target=self.parse, args=[entry]).start()
            while active_count() > 50:
                time.sleep(1)
        while active_count() > idleness: #no idea why it is two not one...
            print('DEBUG: waiting...', active_count(), idleness)
            print(self._lock)
            print(self._semaphore)
            time.sleep(1)
        # final touches to the whole sets...
        for org in list(self._organism_greater_namedex.keys()):
            namedex = {**self._organism_lesser_namedex[org], **self._organism_greater_namedex[org]}
            # cleanup
            for k in ('\n      ', '\n     ', '\n    ', '\n   ', '', '\n', ' '):
                if k in namedex:
                    del namedex[k]
            #namedex = {k.lower(): namedex[k] for k in namedex}
            fn = os.path.join(Protein.settings.dictionary_folder, f'taxid{org}-names2{chosen_attribute}.json')
            json.dump(namedex, open(fn, 'w'))  # name to pdbs
        fn = os.path.join(Protein.settings.dictionary_folder, f'organism.json')
        json.dump(self._organismdex, open(fn, 'w'))  # organism to taxid
        # lighten
        for dex, fn in ((self._uniprot_pdbdex, 'uniprot2pdb.json'),
                        (self._uniprot_namedex, 'uniprot2name.json'),
                        (self._uniprot_speciesdex, 'uniprot2species.json')):
            fp = os.path.join(Protein.settings.dictionary_folder, fn)
            json.dump({k: dex[k] for k in dex if dex[k]}, open(fp, 'w'))

    def parse(self, entry):
        ## parser...
        #print('waiting for semaphore')
        self._semaphore.acquire()
        #print('waited for semaphore')
        prot = Protein.from_uniprot(entry)
        if prot.uniprot in self._uniprot_datasetdex:
            #print('Repeated uniprot!', prot.uniprot)
            self._semaphore.release()
            return None
        #print('getting offsets')
        prot.get_offsets()
        #.get_resolutions() too slow.
        self.get_resolutions_for_prot(prot)
        if prot.organism['common'] == 'Human':
            prot.parse_swissmodel()
            pass
        prot.compute_params()
        ## dict
        chosen_name = getattr(prot, self.chosen_attribute)
        # update the organism dex
        org = prot.organism['NCBI Taxonomy']
        #print('waiting for lock')
        self._lock.acquire()
        #print('waited for lock')
        try:
            prot.dump()  # gdump??
            #print('dumping')
            if org not in self._organismdex:
                for k in prot.organism:
                    self._organismdex[prot.organism[k]] = org
            # make dictionaries...
            self._uniprot_datasetdex[prot.uniprot] = prot.uniprot_dataset
            self._uniprot_pdbdex[prot.uniprot].extend([p.id for p in prot.pdbs])
            if prot.gene_name and prot.gene_name in self._organism_greater_namedex[org]:
                if prot.organism['scientific'] == 'Homo sapiens':
                    print('#' * 20)
                    print('CLASH!!!!!', prot.gene_name, prot.uniprot, prot.organism['scientific'],
                          prot.recommended_name,
                          self._uniprot_datasetdex[prot.uniprot],
                          self._uniprot_datasetdex[self._organism_greater_namedex[org][prot.gene_name]])
                if prot.uniprot_dataset == 'TrEMBL' and self._uniprot_datasetdex[
                    self._organism_greater_namedex[org][prot.gene_name]] == 'Swiss-Prot':
                    self._lock.release()
                    return None
            if prot.recommended_name and prot.recommended_name in self._organism_greater_namedex[org]:
                if prot.organism['scientific'] == 'Homo sapiens':
                    print('#' * 20)
                    print('CLASH!!!!!', prot.gene_name, prot.uniprot, prot.recommended_name,
                          self._uniprot_datasetdex[prot.uniprot],
                          self._uniprot_datasetdex[self._organism_greater_namedex[org][prot.recommended_name]])
                if prot.uniprot_dataset == 'TrEMBL' and self._uniprot_datasetdex[
                    self._organism_greater_namedex[org][prot.recommended_name]] == 'Swiss-Prot':
                    self._lock.release()
                    return None
            self._organism_greater_namedex[org][prot.uniprot] = chosen_name
            self._organism_greater_namedex[org][prot.uniprot_name] = chosen_name
            self._organism_greater_namedex[org][prot.recommended_name] = chosen_name
            self._organism_greater_namedex[org][prot.gene_name] = chosen_name
            self._uniprot_namedex[prot.uniprot] = prot.recommended_name
            self._uniprot_speciesdex[prot.uniprot] = org
            for group in [prot.alt_gene_name_list, prot.alternative_fullname_list, prot.alternative_shortname_list]:
                for name in group:
                    if re.match('[\d\-]+\.[\d\-]+\.[\d\-]+\.[\d\-]+', name):
                        continue  # no EC numbers!
                    self._organism_lesser_namedex[org][name] = chosen_name
            self._lock.release()
            self._semaphore.release()
            #print(prot.uniprot)
            return None
        except Exception as error:
            print(error.__class__.__name__, str(error))
            self._lock.release()
            self._semaphore.release()



    def get_resolutions_for_prot(self, prot):
        # modified from lookup_resolution
        for pdb in prot.pdbs:
            if pdb.type != 'rcsb':
                pass
            elif pdb.code in self.resolutions:
                pdb.resolution = self.resolutions[pdb.code]
            else:
                warn(f'No resolution info for {pdb.code}')
