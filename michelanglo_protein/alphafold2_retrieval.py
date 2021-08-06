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

def is_alphafold_taxon(taxid:int) -> bool:
    try:
        if int(taxid) in alphamodels.values():
            return True
        else:
            return False
    except Exception:
        return False
