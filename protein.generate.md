# protein.generate package

## Submodules

## protein.generate.ET_monkeypatch module


#### class protein.generate.ET_monkeypatch.NewElement(tag, attrib={}, \*\*extra)
Bases: `xml.etree.ElementTree.Element`

This is a collection of methods that helps handle better the ET.Element instnaces. They are monkeypatched to the class object itself.


### describe()

### get_attr(key)

### get_sub_by_type(type_attr)
Gets only the first subtag.


### get_subtag(tag)
Gets only the first subtag.


### has_attr(key, value=None)

### has_text()

### is_human()

### is_tag(tag)

### ns_strip(ns='{http://uniprot.org/uniprot}')
## protein.generate.PDB_blast module


#### class protein.generate.PDB_blast.Blaster()
Bases: `object`

This is just a container.


### static extract_db()

### classmethod full_blaster(outfolder_name, db)
Given the list of genes in in seqdex.json. do a blast against pdbaa from NCBI ftp.
:return:


### classmethod make_fastas()

### static parse(infolder, outfolder)

### static part_blaster(todo)

### classmethod pdb_blaster()

### classmethod self_blaster()
## protein.generate.protParam_mod module


#### protein.generate.protParam_mod.mod(sequence)
This is a not implemented function. It is a fix for ProtParam.ProteinAnalysis().protein_scale and the DIWV scale.
As the latter requires knowldge of the preceeding amino acid it will fail.
>>> p = ProtParam.ProteinAnalysis(sequence)
>>> p.protein_scale(ProtParamData.DIWV, window=9, edge=.4)
hashtag epicfail.
So this is the repalacement.
:param sequence: sequence to score
:type sequence: str
:return: DIWV score.
:rtype: list[int]

## protein.generate.split_gnomAD module

Entry point:

```python
>>> gnomAD().split().write()
```


#### class protein.generate.split_gnomAD.gnomAD()
Bases: `object`


### \__init__()
Instantiation starts the settings.
but the settings can be changed.
split splits the file into the self.data dict containing gene acc id as key and list of gnomADVariant.
But the bound method write writes and the gnomADVariant as regular dictionary.


### split()

### write(folder='gnomAD')

#### class protein.generate.split_gnomAD.gnomADVariant(symbol, identifier, from_residue, residue_index, to_residue, impact, count, homozygous)
Bases: `object`

This is the same as the namedtuple but with more stuff. It does not get written. to_dict does.


### \__init__(symbol, identifier, from_residue, residue_index, to_residue, impact, count, homozygous)
Initialize self.  See help(type(self)) for accurate signature.


### classmethod from_line(line)

### static parse_line(line)

### to_dict()
## protein.generate.split_phosphosite module

The old code was slow:

> def get_PTM(self):

>     assert self.uniprot, ‘Uniprot Acc. required. Kind of.’
>     modified_residues = []
>     for f in os.listdir(self.settings.reference_folder):

>     > if ‘_site_dataset’ in f and ‘.gz’ not in f:  # it’s a Phosphosite plus entry.

>     >     with open(os.path.join(self.settings.reference_folder, f)) as fh:

>     >         next(fh)  # date
>     >         next(fh)  # licence
>     >         next(fh)  # blankline
>     >         for row in csv.DictReader(fh, delimiter=’   ‘):

>     >         > if row[‘ACC_ID’] == self.uniprot: ## this will not pick up mice!

>     >         >     modified_residues.append(row[“MOD_RSD”])

>     self.features[‘PSP_modified_residues’] = modified_residues ## list of str (e.g. ‘K30-m2’)

Unfortunately, you have to agree to the CC-by licence at [https://www.phosphosite.org/staticDownloads](https://www.phosphosite.org/staticDownloads) at phosphosite.
Then you have to manually download all the files with _site_dataset.gz.
Afterwards run the settings method retrieve_references.


#### class protein.generate.split_phosphosite.Phoshosite()
Bases: `object`


### \__init__()
Initialize self.  See help(type(self)) for accurate signature.


### split()

### write(folder='phosphosite')
## protein.generate.uniprot_master_parser module

This file parses the uniprot FTP file and can do various things. such as making a small one that is only human.
But mainly the UniprotMasterReader.convert(‘uniprot_sprot.xml’) method whcih generates the JSON files required. In future these will be databases…
Be warned that ET.Element is a monkeypatched version.


#### class protein.generate.uniprot_master_parser.UniprotMasterReader(uniprot_master_file=None, first_n_protein=0, chosen_attribute='uniprot')
Bases: `object`

see generator iter_human
NB. The ET.Element has been expanded. See help(ElementalExpansion)


### \__init__(uniprot_master_file=None, first_n_protein=0, chosen_attribute='uniprot')
THIS IS FOR MICHELANGLO
:param uniprot_master_file:
:param first_n_protein: set to zero for all, to interger to get the first n.
:return:


### classmethod convert(uniprot_master_file=None, first_n_protein=0)
DO NOT USE!!!
:param uniprot_master_file:
:param first_n_protein: set to zero for all, to interger to get the first n.
:return:


### iter_all(dataset=None)
dataset = Swiss-Prot is better than TrEMBL
Interates across a LARGE Uniprot XML file and returns entries regardless of humanity.
:return: ET.Element()


### iter_human()
Interates across a LARGE Uniprot XML file and returns *only* the humans.
:return: ET.Element()


### parse_human()

### shrink(outfile='human_proteome.xml')
Make a smaller XML file, but with only the human proteome.
:param outfile:
:return:

## Module contents
