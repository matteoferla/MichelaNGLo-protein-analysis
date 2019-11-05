from protein import ProteinAnalyser, ProteinCore, Mutation, Structure
from protein.settings_handler import global_settings
from protein.generate import ProteinGatherer, ProteomeGatherer, split_gNOMAD
from protein.protein_analysis import StructureAnalyser

import pickle


def test_ProteinAnalyser():
    p = ProteinAnalyser(uniprot = ' Q86V25').load()
    print(p)
    p.mutation = Mutation('p.N127W')
    p.analyse_structure()
    print(p.get_features_near_position())
    print(p.get_gNOMAD_near_position())
    print(p.model.get_structure_neighbours())
    print(p2.get_superficiality())


# p=ProteinGatherer(uniprot='Q6ZN55').parse_uniprot().parse_pdb_blast()

# from protein.apriori_effect import WikiTable
# print(WikiTable(WikiTable.grantham).ndata)

def main():
    ## make everything!

    global_settings.error_tolerant = True

    ProteomeGatherer(skip=True, remake_pickles=True)

from protein.generate._proteome_gatherer2 import UniprotReader
import os, json
def mini_gene_data():
    genes = '''DOCK180
    DOCK2
    DOCK3
    DOCK4
    DOCK5
    DOCK6
    DOCK7
    DOCK8
    DOCK9
    DOCK10
    DOCK11
    '''.split()


    data = {}
    from pprint import PrettyPrinter
    pprint = PrettyPrinter().pprint
    namedex = json.load(open('data/human_prot_namedex.json'))
    for uni in set(namedex.values()):
        g = ProteinGatherer(uniprot=uni).parse_uniprot()
        data[g.gene_name] = {'name': g.gene_name, 'uniprot': g.uniprot, 'len': len(g), 'domains': {k: g.features[k] for k in ('active site','modified residue','topological domain','domain','region of interest','transmembrane region') if k in g.features}, 'disease': g.diseases}
        #print(g.gene_name,g.uniprot,len(g))
    json.dump(data,open('map.json','w'))

def make_pdb_dex():
    #I need to make a uniprot to pdb dex.
    from protein.generate._proteome_gatherer2 import UniprotReader
    master_file = os.path.join(ProteinGatherer.settings.temp_folder, 'uniprot_sprot.xml')
    UniprotReader.make_dictionary(uniprot_master_file=master_file, first_n_protein=0, chosen_attribute='uniprot')

def iterate_taxon(taxid):
    path = os.path.join(global_settings.pickle_folder,f'taxid{taxid}')
    for pf in os.listdir(path):
        try:
            protein = ProteinGatherer().load(file=os.path.join(path, pf))
            protein.gNOMAD = []
            protein.parse_gNOMAD()
            protein.dump()
            #protein.get_offsets().parse_gNOMAD().compute_params()
            #protein.dump()
        except:
            pass


if __name__ == '__main__':
    global_settings.verbose = False
    global_settings.init(data_folder='../protein-data')
        #.retrieve_references(ask=False, refresh=False)
    #UniprotReader()

    #global_settings.init()

    #make_pdb_dex()
    split_gNOMAD.gNOMAD()
    iterate_taxon('9606')

    p = ProteinAnalyser(taxid='9606', uniprot='Q9BZ29').load()
    p.mutation = 'P23W'
    print(p.check_mutation())
    print(p.mutation_discrepancy())
    print(p.predict_effect())
    print(p.elmdata)
    print(p._elmdata)

    # fetch_binders is too slow. Pre-split the data like for gnomad.

if __name__ == '__main__':
    global_settings.verbose = False
    global_settings.init(data_folder='../protein-data')

if 1 == 0:
    s = Structure(id='2WM9', description='', x=-1, y=-1, code='2WM9').lookup_sifts()


if __name__ == '__main__':
    #test_ProteinAnalyser()
    p = ProteinCore(uniprot='Q96N67').load()

    def ncbize(n):
        if n < 1419:
            return n
        elif n < 1832:
            return n+9
        else:
            return n+11

    s = [str(ncbize(v.x)) for v in p.gNOMAD]
    for i in range(len(s)//100 + 1):
        print('color yellow, chain A and resi '+'+'.join(s[100*i:100*(i+1)]))
    print('color red, chain A and resi ' + '+'.join([str(ncbize(v.x)) for v in p.gNOMAD if v.impact == 'HIGH']))
