## Protein module of VENUS (SNV analsyer)

The `protein` module collects all the data needed for analysing variants.

    >>> from protein import Protein
    >>> Protein(uniprot = 'Q9NWZ3').parse_all('parallel')

### files within protein module

The `protein` module's `__init__.py` file contains the `Protein` and `Mutation` classes. But `Protein` is so big it is split across three files, it's mixin base classes are in `protein._protein_uniprot_mixin` (handles uniprot) and `protein._protein_base_mixin.py` (magic methods).

The `_UniprotMixin` requires Element-tree to be monkeypatched, which is done in `ET_monkeypatch`.

    >>> from ET_monkeypatch import ET

The where-is-what logistics and other settings is stored in `Protein.settings`, which is an instance of `protein.settings_handler.GlobalSettings()`.

    >>> from protein import Protein
    >>> Protein.settings.verbose = True
    >>> Protein.settings.missing_attribute_tolerant = True
    >>> Protein(uniprot = 'Q9NWZ3').foo
    warning... but does not die.

For where the error tolerance gets applied see the mixin `protein._protein_base_mixin._BaseMixin`

### Slimming

The folder protein is the module. It has two classes Protein and Mutation.

Unfortunately it has gotten bloated. So it needs tidying.

The folder `ref` is the starting data.

The folder `temp` is where temp files to make the final data

The folder `data` contains the final data.

### Future ideas
For what is needed for VENUS (SNV analyser see that)

The uniprot parser part of this module could be integrated with Biopython. Currently there is no Uniprot parsing capability for it.
