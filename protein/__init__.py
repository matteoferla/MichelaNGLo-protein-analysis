__docs__ = """
Protein was formerly called Variant.
This module has three classes.
* ProteinLite, which handles Protein that have a pickled dataset
* Protein, which inherits the former and fetches data to make the datasets.
* Mutation handles the mutation
* the `global_settings` variable, declaired in `.settings_handler` has handles the config stuff.

Additionally, the Protein class requires mixin classes from `_protein_*_mixin.py`. The XML parser for Uniprot requires a special ET from `.ET_monkeypatched`.

The submodule generate has a method `generate`, which generates all the datasets for the human proteome. Although the files required are in settings_handler

The Mutation class uses a variable that was originally generated in _apriori_effect.py.
"""

from ._protein_lite import ProteinLite
from ._protein_full import Protein
from ._mutation import Mutation

from .ET_monkeypatch import ET #monkeypatched version
from .settings_handler import global_settings
