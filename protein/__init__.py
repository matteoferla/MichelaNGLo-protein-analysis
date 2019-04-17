__doc__ = """
This module has several classes.
* ProteinCore is used by the rest and provides the backbone. It can read and write but not generate protein date. for that there is
* generate.ProteinGatherer, which parses data from various sources. generate.ProteomeGatherer starts everything up and parses the whole proteome.
* Mutation handles the mutation
* the `global_settings` variable, declaired in `.settings_handler` has handles the config stuff.

The submodule generate has a method `generate`, which generates all the datasets for the human proteome. Although the files required are in settings_handler

The Mutation class uses a variable that was originally generated in _apriori_effect.py.

The script unitests.py does... unitests.


The script protein.aprior_effect generates the dictionary that is used to say what the apriori effect are. Namely, what amino acid is smaller etc.
This is already added to protein, so does not need to be redone.
>>> from protein.apriori_effect import Changedex
>>> pprint(Changedex().fill().to_dict())

To start everything going...
>>> from protein.generate.ProteomeGatherer
>>> ProteomeGatherer()
and wait for the missing file. There surely is one. generate.__init__ has a list of references.

"""

from .core import ProteinCore
from .protein_analysis import ProteinAnalyser
from .mutation import Mutation
from .settings_handler import global_settings
