from protein import ProteinAnalyser, ProteinCore, Mutation, Structure
from protein.settings_handler import global_settings
from protein.generate import ProteinGatherer, ProteomeGatherer, split_gnomAD
from protein.protein_analysis import StructureAnalyser

import pickle


def test_ProteinAnalyser():
    p = ProteinAnalyser(uniprot = ' Q86V25').load()
    print(p)
    p.mutation = Mutation('p.N127W')
    p.analyse_structure()
    print(p.get_features_near_position())
    print(p.get_gnomAD_near_position())
    print(p.model.get_structure_neighbours())
    print(p2.get_superficiality())


# p=ProteinGatherer(uniprot='Q6ZN55').parse_uniprot().parse_pdb_blast()

# from protein.apriori_effect import WikiTable
# print(WikiTable(WikiTable.grantham).ndata)

def main():
    ## make everything!

    global_settings.error_tolerant = True

    ProteomeGatherer(skip=True, remake_pickles=True)

from protein.generate.uniprot_master_parser import UniprotReader
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
    from protein.generate.uniprot_master_parser import UniprotReader
    master_file = os.path.join(ProteinGatherer.settings.temp_folder, 'uniprot_sprot.xml')
    UniprotReader.make_dictionary(uniprot_master_file=master_file, first_n_protein=0, chosen_attribute='uniprot')

def iterate_taxon(taxid=9606):
    """
    This is an ad hoc fix to fix humans or similar. For full deployment use ProteomeParser.
    :param taxid:
    :return:
    """
    path = os.path.join(global_settings.pickle_folder,f'taxid{taxid}')
    for pf in os.listdir(path):
        try:
            protein = ProteinGatherer().load(file=os.path.join(path, pf))
            protein.gnomAD = []
            protein.parse_gnomAD()
            protein.get_PTM()
            protein.compute_params()
            protein.dump()
            #protein.get_offsets().parse_gnomAD().compute_params()
            #protein.dump()
        except:
            pass


if __name__ == '__main__':
    global_settings.verbose = True #False
    global_settings.startup(data_folder='../protein-data')
if 1==1:
    from protein.generate.split_phosphosite import Phoshosite
    #ph = Phoshosite().split().write('phosphosite')
    p = ProteinGatherer(taxid='9606', uniprot='P62879').load()
    print(':B and ('+' or '.join([str(m.x) for m in p.gnomAD if m.homozygous])+')')

    print([m for m in p.gnomAD if m.homozygous])
    print(' '.join([str(m.description.split()[0]) for m in p.gnomAD])+')')
    print([m for m in p.gnomAD if m.homozygous])
elif 1==0:
    iterate_taxon('9606')
        #.retrieve_references(ask=False, refresh=False)
    #UniprotReader()

    #global_settings.startup()

    #make_pdb_dex()
    #split_gnomAD.gnomAD().write()
    #iterate_taxon('9606')

    p = ProteinAnalyser(taxid='9606', uniprot='Q9BZ29').load()
    p.mutation = 'P23W'
    print('check_mutation', p.check_mutation())
    print('mutation_discrepancy',p.mutation_discrepancy())
    print('predict_effect', p.predict_effect())
    print('elmdata', p.elmdata)
    print(p.mutation)
    # fetch_binders is too slow. Pre-split the data like for gnomAD.
elif 1==0:
    print('retrieving...')
    global_settings.retrieve_references(ask=False, refresh=False)
else:
    #test_ProteinAnalyser()
    p = ProteinGatherer(uniprot='Q96N67').load()
    p.parse_gnomAD()
    def ncbize(n):
        if n < 1419:
            return n
        elif n < 1832:
            return n+9
        else:
            return n+11
    s = [str(ncbize(v.x)) for v in p.gnomAD]
    for i in range(len(s)//100 + 1):
        print('color yellow, chain A and resi '+'+'.join(s[100*i:100*(i+1)]))
    print('color orange, chain A and resi ' + '+'.join([str(ncbize(v.x)) for v in p.gnomAD if v.impact == 'HIGH']))
    print('color red, chain A and resi '+'+'.join([str(ncbize(v.x)) for v in p.gnomAD if v.homozygous]))
