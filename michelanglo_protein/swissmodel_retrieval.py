from typing import *
import requests
from .structure import Structure
from .metadata_from_PDBe import PDBMeta
import logging

log = logging.getLogger()

"""
The data return from a swissmodel query is full of information but is as complicated as a uniprot one...

.. code-block:: python
    fs = FromSwissmodel('NDKB_HUMAN')
    # ProteinCore(uniprot='NDKB_HUMAN')
    data = fs.retrieve_structures_from_swissmodel()
    data.keys() # dict_keys(['api_version', 'query', 'query_date', 'result'])
    data['result'].keys() # dict_keys(['crc64', 'md5', 'sequence', 'sequence_length', 'structures', 'uniprot_entries'])
    data['result']['uniprot_entries'] # [{'ac': 'P22392', 'id': 'NDKB_HUMAN', 'isoid': 1}]
    data['result']['structures'] # list one per structure
    data['result']['structures'][0].keys() # dict_keys(['chains', 'coordinates', 'coverage', 'created_date', 'from', 'in_complex_with', 'ligand_chains', 'md5', 'method', 'oligo-state', 'provider', 'template', 'to'])
    get_key = lambda entry: {key: entry[key] for key in ('method', 'oligo-state','provider','template')}
    for entry in data['result']['structures']:
        print(get_key(entry))
    #{'method': 'X-RAY DIFFRACTION', 'oligo-state': 'homo-6-mer', 'provider': 'PDB', 'template': '1nsk'}
    #{'method': 'HOMOLOGY MODELLING', 'oligo-state': 'homo-6-mer', 'provider': 'SWISSMODEL', 'template': '1jxv.1.A'}

The parsing is done in `Structure.from_swissmodel_query()
"""

# regen pdb and swissmodel data
class FromSwissmodel:
    def __init__(self, uniprot:str):
        # do not inherit. here solely for mental clarity
        self.uniprot = uniprot
        # {'description': elem.attrib['id'], 'id': elem.attrib['id'], 'x': loca[0], 'y': loca[1]}
        self.pdbs = []
        self.swissmodel = []

    def retrieve_structures_from_swissmodel(self, blank_previous: bool = True) -> dict:
        pdbs = []
        swissmodel = []
        url = f'https://swissmodel.expasy.org/repository/uniprot/{self.uniprot}.json'
        reply = requests.get(url,
                             # params=dict(provider='swissmodel')  # do both.
                             )
        if reply.status_code != 200:
            msg = f'Swissmodel retrieval failed (code {reply.status_code}): {url}'
            log.info(msg)
            raise ConnectionError(msg)
        data = reply.json()
        for structural_data in data['result']['structures']:
            structure = Structure.from_swissmodel_query(structural_data, self.uniprot)
            if structure.type == 'rcsb':
                pdbs.append(structure)
            else:
                swissmodel.append(structure)
        log.debug(f'{len(pdbs)} PDBs and {len(swissmodel)} SWISSMODELs available')
        # This does not used the cached data.
        res = PDBMeta.bulk_resolution([structure.code for structure in pdbs if len(structure.code) == 4])
        for structure in pdbs:
            if structure.code in res:
                structure.resolution = res[structure.code.lower()]
            else:
                log.debug(f'PDBe resolution not found for {structure.code} (options: {list(res.keys())}')
        # add...
        if blank_previous:
            self.pdbs = pdbs
            self.swissmodel = swissmodel
        else:
            self.pdbs.extend(pdbs)
            self.swissmodel.extend(swissmodel)
        return data
