## Protein module of VENUS (SNV analsyer)

The `protein` module collects all the data needed for analysing variants.

    >>> from protein import Protein
    >>> Protein(uniprot = 'Q9NWZ3').parse_all('parallel')

### files within protein module

The `protein` module's `__init__.py` file contains the `Protein` and `Mutation` classes. But `Protein` is so big it is split across three files, it's mixin base classes are in `protein._protein_uniprot_mixin` (handles uniprot) and `protein._protein_base_mixin.py` (magic methods except `__init__`).

The `_UniprotMixin` requires Element-tree to be monkeypatched, which is done in `ET_monkeypatch`.

    >>> from ET_monkeypatch import ET

The where-is-what logistics and other settings is stored in `Protein.settings`, which is an instance of `protein.settings_handler.GlobalSettings()`.

    >>> from protein import Protein
    >>> Protein.settings.verbose = True
    >>> Protein.settings.missing_attribute_tolerant = True
    >>> Protein(uniprot = 'Q9NWZ3').foo
    warning... but does not die.
    >>> Protein.settings.error_tolerant = True ##see the decorator failsafe for more.
    >>> Protein(uniprot = 'foo').parse_uniprot() #will fail (extreme case)
    warning... but does not die
    >>> Protein.settings.data_folder = 'new_location' ## all subfolders will be created too.

For where the error tolerance gets applied see the mixin `protein._protein_base_mixin._BaseMixin`

This script used to be called `Variant` and generated static html pages (_cf._ HIFC2-Tracker_variant_classifier repo). It still can (maybe).
The pages folder is `.page_folder` and `.wipe_html()` clears them.

To change species see `ET_monkeypatch.ET.is_human()`.

### Data

The data is a big. It needs slimming as unfortunately it has gotten bloated. So it needs tidying.

The module `protein.generate` creates the needed files.

At some point I started writing `settings.retrieve_references()` but stopped as it would have been for the once and some reference require passwords and some change faster than stairways in Hogwarts (_ie._ NCBI).


NB. The script without the data will fetch off the web the Uniprot (only `.parse_uniprot()` will work without screaming).

#### Current requirements that are implemented

| Resource | Location | Reason | Last checked |
| --- | --- | --- | --- |
| Uniprot | [ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz) | Most of the data comes from Uniprot | 15/03/2019 |
| NCBI PDB BLAST DB | [ftp://ftp.ncbi.nlm.nih.gov/blast/db/pdbaa.tar.gz](ftp://ftp.ncbi.nlm.nih.gov/blast/db/pdbaa.tar.gz)  | Needed to find homogues with crystal structures | 15/03/2019 |


If any changes please change in `protein/settings_handler` the class attribute`GlobalSettings.addresses` (or `Protein.settings.addresses`).

### Future ideas
For what is needed for VENUS (SNV analyser see that)

The uniprot parser part of this module could be integrated with Biopython. Currently there is no Uniprot parsing capability for it.
