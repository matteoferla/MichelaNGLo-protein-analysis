from michelanglo_protein import ProteinAnalyser, ProteinCore, Mutation, structure
from michelanglo_protein.settings_handler import global_settings
from michelanglo_protein.generate import ProteinGatherer, ProteomeGatherer
from michelanglo_protein.generate.split_gnomAD import gnomAD
from michelanglo_protein.protein_analysis import StructureAnalyser
# Settings = namedtuple('settings', 'dictionary_folder', 'reference_folder', 'temp_folder')
import pickle
import sys, traceback

from pprint import PrettyPrinter
pprint = PrettyPrinter().pprint


def test_ProteinAnalyser():
    p = ProteinAnalyser(uniprot = 'Q86V25').load()
    print(p)
    p.mutation = Mutation('p.N127W')
    p.analyse_structure()
    print(p.get_features_near_position())
    print(p.get_gnomAD_near_position())
    print(p.model.get_structure_neighbours())
    print(p.get_superficiality())


# p=ProteinGatherer(uniprot='Q6ZN55').parse_uniprot().parse_pdb_blast()

# from michelanglo_protein.apriori_effect import WikiTable
# print(WikiTable(WikiTable.grantham).ndata)


from michelanglo_protein.generate.uniprot_master_parser import UniprotMasterReader
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
    from michelanglo_protein.generate.uniprot_master_parser import UniprotMasterReader
    master_file = os.path.join(ProteinGatherer.settings.temp_folder, 'uniprot_sprot.xml')
    UniprotMasterReader.make_dictionary(uniprot_master_file=master_file, chosen_attribute='uniprot')

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
            #michelanglo_protein.get_offsets().parse_gnomAD().compute_params()
            #michelanglo_protein.dump()
        except:
            pass

def how_many_empty(taxid=9606):
    from collections import Counter
    global_settings.verbose = False
    empty = 0
    full = 0
    path = os.path.join(global_settings.pickle_folder, f'taxid{taxid}')
    for pf in os.listdir(path):
        p = ProteinGatherer().load(file=os.path.join(path, pf))
        if len(p.sequence) == 0:
            print(p)
            empty += 1
        else:
            full +=1
    print(full, empty)

def compress(taxid=9606, target='../gpickle'):
    if not os.path.exists(target):
        os.mkdir(target)
    source = os.path.join(global_settings.pickle_folder, f'taxid{taxid}')
    for pf in os.listdir(source):
        p = ProteinCore().load(file=os.path.join(source, pf))
        p.gdump(file=os.path.join(target, os.path.splitext(pf)[0]+'.gzp'))


def fix_empty(taxid=9606):
    from collections import Counter
    global_settings.verbose = False
    glitchy = 0
    fine = 0
    fixed = 0
    path = os.path.join(global_settings.pickle_folder, f'taxid{taxid}')
    for pf in os.listdir(path):
        p = ProteinGatherer().load(file=os.path.join(path, pf))
        if len(p.sequence) == 0:
            print('****************************************')
            print(f'Attempting to fix {p.gene_name}')
            try:
                global_settings.verbose = True
                p.parse_uniprot()
                p.parse_swissmodel()
                p.compute_params()
                p.parse_gnomAD()
                p.get_PTM()
                assert len(p.sequence) > 0, 'Darn. Sequence is zero AA long'
                p.dump()
                fixed += 1
                global_settings.verbose = False
            except Exception:
                traceback.print_exc(file=sys.stdout)
                glitchy += 1
        else:
            fine +=1
    print('****************************************')
    print(f'Fine: {fine:,}, Fixed {fixed:,}, Glitchy: {glitchy:,}')


def describe(uniprot):
    print('***************** DESCRIPTION *******************************')
    p = ProteinCore(taxid='9606', uniprot=uniprot).load()  # gnb1 P62873 gnb2 P62879
    pprint(p.asdict())

def jsonable(self):
    def deobjectify(x):
        if isinstance(x, dict):
            return {k: deobjectify(x[k]) for k in x}
        elif isinstance(x, list) or isinstance(x, set):
            return [deobjectify(v) for v in x]
        elif isinstance(x, int) or isinstance(x, float):
            return x
        else:
            return str(x) # really ought to deal with falseys.
    return {a: deobjectify(getattr(self, a, '')) for a in self.__dict__}

def analyse(uniprot):
    print('***************** ANALYSIS *******************************')
    p = ProteinAnalyser(taxid='9606', uniprot=uniprot).load()
    p.mutation = f'{p.sequence[65]}66W'
    p.predict_effect()
    p.analyse_structure()
    print(p)
    import json
    json.dump(p.asdict(), open('test.json','w'))

    #print('Best one: ',p.get_best_model())
    # print('analysed', {**jsonable(p.structural),
    #     'superficiality': p.structural.get_superficiality(),
    #     'structural_neighbours': list(p.structural.get_structure_neighbours())})
    # http://0.0.0.0:8088/venus_analyse?uniprot=P62879&species=9606&mutation=A73T

def reparse_gene(name):
    human = json.load(open(os.path.join(global_settings.dictionary_folder, 'taxid9606-names2uniprot.json')))
    target = human[name]
    p = ProteinGatherer(uniprot=target)
    p.parse_uniprot()
    print(p.sequence)

if __name__ == '__main__':
    global_settings.verbose = True #False
    #global_settings.startup(data_folder='../MichelaNGLo-protein-data')
    global_settings.startup(data_folder='../protein-data')
#### workspace!
if 1==1:
    #describe('P01112')
    #analyse('P62873')
    #how_many_empty()
    #fix_empty()
    compress()

    #parse_uniprot(
elif 1==9:
    p = ProteinGatherer(taxid='9606', uniprot='P62873').load()
    print(p.gnomAD)
    print(p.parse_gnomAD())
    print(p.gnomAD)
    print(p.features['PSP_modified_residues'])
    from michelanglo_protein.generate.split_phosphosite import Phoshosite
    #ph = Phoshosite().split().write('phosphosite')
    p = ProteinGatherer(taxid='9606', uniprot='P62879').load()
    print(':B and ('+' or '.join([str(m.x) for m in p.gnomAD if m.homozygous])+')')

    print([m for m in p.gnomAD if m.homozygous])
    print(' '.join([str(m.description.split()[0]) for m in p.gnomAD])+')')
    print([m for m in p.gnomAD if m.homozygous])
elif 1==0:
    iterate_taxon('9606')
        #.retrieve_references(ask=False, refresh=False)
    #UniprotMasterReader()

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
    ## dock 9 ops.
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
