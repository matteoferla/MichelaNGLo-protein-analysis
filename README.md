## Protein module for Michelanglo and VENUS (formerly SNV analyser)

The `michelanglo_protein` module collects all the data needed for analysing variants.
    
    p = ProteinAnalyser(uniprot = 'Q86V25').load()
    print(p)
    p.mutation = Mutation('p.N127W')
    p.analyse_structure()
    print(p.get_features_near_position())
    print(p.get_gnomAD_near_position())
    print(p.model.get_structure_neighbours())
    print(p.get_superficiality())
    
The data loaded is either gatherered from various databases, some of which need splitting (_vide infra_ or `create.py`).
If it is just the one gene, you can use the following, which will retrieve the Uniprot data of the one gene (failing to get external data if unavailable):

    p = ProteinGatherer(uniprot='Q86V25')
    p.parse_uniprot()
    

The site [michelanglo.sgc.ox.ac.uk](https://michelanglo.sgc.ox.ac.uk) depends on three repos:
* [MichelaNGLo-app](https://github.com/matteoferla/MichelaNGLo)
* [MichelaNGLo-transpiler](https://github.com/matteoferla/MichelaNGLo-transpiler)
* **MichelaNGLo-protein-module**

This module can be used independently of Michelanglo app module, but requires the transpiler module.

If you are interested in the gene/protein name synonyms to uniprot mapping files see 
[Name synomyms to Uniprot](https://github.com/matteoferla/Name-synomyms-to-Uniprot).

If you want human protein and don't want to generate the data (i.e. by running `create.py`) see [human protein](https://github.com/matteoferla/MichelaNGLo-human-protein-data).

For installation instructions see relevant sections in [Michelanglo deployment](https://github.com/matteoferla/MichelaNGLo-app/blob/master/git_docs/deploy.md).

### files within protein module

The `michelanglo_protein` module contains a `.generate` subfolder.
For speed at time of a request the module preparses a lot of things thanks to the `.generate` submodule.
Initially the Uniprot data was fetched online and parsed on the fly —technically this should still be possible with a few tweaks.

There are a few "protein" classes.

* These have different roles, but are all based off `ProteinCore` class.
* The parser protein class from `.generate` is `ProteinGatherer`, this has many fetching parts that are fault tolerant (changeable in settings).
* For analysis of a mutation `ProteinAnalyser` is used.

Additionally, there are classes whose instances are bound to these:

* `.settings`, instance of `GlobalSettings`, a singleton class, which controls the settings and folders and actually fetches the initial data (see below).
* `Variant` instances appear in lists. gnomAD etc.
* whereas `ProteinAnalyser.mutation` is an instance of `Mutation`
* `Structure` instances appear in lists. Stores crystal data.
* `StructureAnalyser` is a special case of the above for `ProteinAnalyser.structure`

## GlobalSettings
GlobalSettings appears as `.settings` of various classes. It is a singleton class, so the changes are global.
It is instatiated in `settings_handler.py`, but it does not create folders etc. until the `.startup` method is called (directly or indirectly by attempting anything).

    global_settings.verbose = True False
    global_settings.startup(data_folder='My-data')
    global_settings.retrieve_references(ask=False, refresh=False)

## Mixin classes making ProteinGatherer
`ProteinGatherer` is so big it is split across three files, it's mixin base classes are in `protein._protein_uniprot_mixin` (handles uniprot) and `protein._protein_base_mixin.py` (magic methods except `__init__`).

The `_UniprotMixin` requires Element-tree to be monkeypatched, which is done in `ET_monkeypatch`.

    >>> from ET_monkeypatch import ET

The where-is-what logistics and other settings is stored in `Protein.settings`, which is an instance of `protein.settings_handler.GlobalSettings()`.

    >>> from protein.generate import ProteinGatherer as Protein
    >>> Protein.settings.startup('temp_just_for_today')  # create the folders.
    >>> Protein.settings.verbose = True
    >>> Protein.settings.missing_attribute_tolerant = True
    >>> Protein(uniprot = 'Q9NWZ3').foo
    warning... but does not die.
    >>> Protein.settings.error_tolerant = True ##see the decorator failsafe for more.
    >>> Protein(uniprot = 'foo').parse_uniprot() #will fail (extreme case)
    warning... but does not die
    >>> Protein.settings.data_folder = 'new_location' ## all subfolders will be created too.

For where the error tolerance gets applied see the mixin `protein._protein_base_mixin._BaseMixin`

The settings does not create on import the folders, but when the `.init(data_folder='data')` method is called or any operations are performed.
To change the folder of `data` midway use `.data_folder = xxx`, don't reinitialise.

This script used to be called `Variant` and generated static html pages (_cf._ HIFC2-Tracker_variant_classifier repo). It still can (maybe).
The pages folder is `.page_folder` and `.wipe_html()` clears them.

To change species see `ET_monkeypatch.ET.is_human()`.

### Data

The data is a big. It needs slimming as unfortunately it has gotten bloated. So it needs tidying.

The module `protein.generate` creates the needed files.

At some point I started writing `settings.retrieve_references()` but stopped as it would have been for the once and some reference require passwords and some change faster than stairways in Hogwarts (_ie._ NCBI).
It mostly works. But downloading large datasets such as gnomAD will run out of RAM as it currently is not streamed.
Also as of December 2019 there is no exome gnomAD v3. So it has to be made from the gargantuan dataset.


NB. The script without the data will fetch off the web the Uniprot (only `.parse_uniprot()` will work without screaming).


#### Current requirements that are implemented

See `protein/settings_handler` the class attribute`GlobalSettings.addresses` (or `Protein.settings.addresses`).

### Offset

The offset of a PDB chain in the coordinate resi relative to Uniprot is based on `SIFT` data. Mostly.

In some cases there is both a missing N and C terminus (see [offset_issue.png](/offset_issue)), which means that there is no PDB_BEG or PDB_END.
This is where the bound method of Structure `.get_offset_from_PDB(detail, sequence)` comes in. This is a fix to an already generated file.
It requires the sequence of the chain in question.

See also: [blog post](http://blog.matteoferla.com/2019/09/pdb-numbering-rollercoaster.html).


### Future ideas

See [Notes](notes.md)

For what is needed for VENUS (SNV analyser see that)

The uniprot parser part of this module could be integrated with Biopython. Currently there is no Uniprot parsing capability for it.


## More!

Sphinx documentation:
* [index.md](./index.md)
* [protein module](./protein.md)
* [protein.generate module](./protein.generate.md)
* [modules](./modules.md)