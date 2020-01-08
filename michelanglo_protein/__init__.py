__doc__ = """
This module has several classes.
* ProteinCore is used by the rest and provides the backbone. It can read and write itself (even in compressed form, see .gdump and .gload) but not generate protein data. for that there is
* generate.ProteinGatherer, which parses data from various sources. generate.ProteomeGatherer starts everything up and parses the whole proteome.
* there are some classes in generate too, but they don't see the light of day. For those see generate._protein_gatherer
* Mutation handles the mutation
* the `global_settings` variable, declaired in `.settings_handler` has handles the config stuff.
* ProteinAnalysis analyses a mutation. 
* Variant and Structure are just namedtuples

The submodule generate has a method `generate`, which generates all the datasets for the human proteome. Although the files required are in settings_handler

The Mutation class is mostly used by ProteinAnalysis but itself holds wordy *_effects attributes and uses a variable that was originally generated in _apriori_effect.py.

The script unitests.py does... unitests.

The script michelanglo_protein.aprior_effect generates the dictionary that is used to say what the apriori effect are. Namely, what amino acid is smaller etc.
This is already added to michelanglo_protein, so does not need to be redone.
>>> from michelanglo_protein.apriori_effect import Changedex
>>> pprint(Changedex().fill().to_dict())

To start everything going...
>>> from michelanglo_protein.generate.ProteomeGatherer
>>> ProteomeGatherer()
and wait for the missing file. There surely is one. generate.__init__ has a list of references.
"""

from .settings_handler import global_settings
from .core import ProteinCore, Variant, Structure
from .protein_analysis import ProteinAnalyser
from .protein_manual import ProteinManual
from .mutation import Mutation

__version__ = '0.4'
