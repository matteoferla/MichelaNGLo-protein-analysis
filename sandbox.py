from protein import ProteinAnalyser, ProteinCore, Mutation, settings_handler
from protein.generate import ProteinGatherer, ProteomeGatherer

p = ProteinAnalyser(uniprot = 'O75015').load()
print(p.get_structures_with_position(124))
