from .structure import Structure

alphamodels = {'Arabidopsis thaliana': 3702,
               'Caenorhabditis elegans': 6239,
               'Candida albicans': 5476,
               'Danio rerio': 7955,
               'Dictyostelium discoideum': 44689,
               'Drosophila melanogaster': 7227,
               'Escherichia coli': 562,
               'Glycine max': 3847,
               'Homo sapiens': 9606,
               'Leishmania infantum': 5671,
               'Methanocaldococcus jannaschii': 243232,
               'Mus musculus': 10090,
               'Mycobacterium tuberculosis': 1773,
               'Oryza sativa': 4530,
               'Plasmodium falciparum': 5833,
               'Rattus norvegicus': 10116,
               'Saccharomyces cerevisiae': 4932,
               'Schizosaccharomyces pombe': 4896,
               'Staphylococcus aureus': 1280,
               'Trypanosoma cruzi': 5693}

# this has to be external as it is used by app.
def is_alphafold_taxon(taxid:int) -> bool:
    try:
        if int(taxid) in alphamodels.values():
            return True
        else:
            return False
    except Exception:
        return False

# --------------------------------------------------------------------------------

class FromAlphaFold2:  # to be inherited by ProteinCore
    def add_alphafold2(self):
        if 'NCBI Taxonomy' not in self.organism:
            return
        if not is_alphafold_taxon(self.organism['NCBI Taxonomy']):
            return
        suffix = 'F1-model_v1' # will this change in future?
        url = f'https://alphafold.ebi.ac.uk/files/AF-{self.uniprot}-{suffix}.pdb'
        structure = Structure(
            # these two do ziltch:
            id='alphafold',
            description=suffix,
            # this is shown by venus
            code=f'AF-{self.uniprot}',
            # range in uniprot
            x=1,
            y=len(self),
            offset=0,
            type='alphafold2',
            url=url,
            chain='A'
        )
        structure.chain_definitions = [{'chain': 'A',
                                           'uniprot': self.uniprot,
                                           'x': 1,
                                           'y': len(self),
                                           'offset': 0,
                                           'range': f'1-{len(self)}',
                                           'name': self.uniprot,
                                           'description': 'AlphaFold2'}]
        self.alphafold2.append(structure)


