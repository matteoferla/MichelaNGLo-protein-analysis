# protein package

## Subpackages


* protein.generate package


    * Submodules


    * protein.generate.ET_monkeypatch module


    * protein.generate.PDB_blast module


    * protein.generate.protParam_mod module


    * protein.generate.split_gnomAD module


    * protein.generate.split_phosphosite module


    * protein.generate.uniprot_master_parser module


    * Module contents


## Submodules

## protein.ET_monkeypatch module

## protein.apriori_effect module

The script protein.aprior_effect generates the dictionary that is used to say what the apriori effect are. Namely, what amino acid is smaller etc.
>>> from protein.apriori_effect import Changedex
>>> pprint(Changedex().fill().to_dict())

Also the scores from wikipedia.


#### class protein.apriori_effect.Changedex()
Bases: `object`


### \__init__()
Initialize self.  See help(type(self)) for accurate signature.


### aa( = ('A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W', 'S', 'T', 'N', 'Q', 'C', 'G', 'P', 'R', 'H', 'K', 'D', 'E'))

### aromatic( = ('P', 'Y', 'W', 'H'))

### fill()

### fill_inverse()

### full( = {'A': 'from a non-aromatic to an aromatic', 'B': 'bigger', 'C': 'differently charged', 'D': 'differently shaped', 'E': 'equally sized', 'F': 'more flexible', 'H': 'more hydrophobic', 'I': 'identical', 'P': 'more polar', 'R': 'more rigid', 'S': 'smaller', 'a': 'from an aromatic to a non-aromatic'})

### hydrophobic( = ('A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W'))

### inverter( = {'A': 'a', 'B': 'S', 'C': 'C', 'D': 'D', 'E': 'E', 'F': 'R', 'H': 'P', 'I': 'I', 'P': 'H', 'R': 'F', 'S': 'B', 'a': 'A'})

### negative( = ('D', 'E'))

### polar( = ('C', 'G', 'P', 'S', 'T', 'N', 'Q'))

### positive( = ('R', 'H', 'K'))

### to_dict()

#### class protein.apriori_effect.WikiTable(table)
Bases: `object`


### \__init__(table)
Initialize self.  See help(type(self)) for accurate signature.


### epstein( = '{| class="wikitable"\\n|\\n|Phe\\n|Met\\n|Leu\\n|Ile\\n|Val\\n|Pro\\n|Tyr\\n|Trp\\n|Cys\\n|Ala\\n|Gly\\n|Ser\\n|Thr\\n|His\\n|Glu\\n|Gln\\n|Asp\\n|Asn\\n|Lys\\n|Arg\\n|-\\n|Phe\\n|\\n|0.05\\n|0.08\\n|0.08\\n|0.1\\n|0.1\\n|0.21\\n|0.25\\n|0.22\\n|0.43\\n|0.53\\n|0.81\\n|0.81\\n|0.8\\n|1\\n|1\\n|1\\n|1\\n|1\\n|1\\n|-\\n|Met\\n|0.1\\n|\\n|0.03\\n|0.03\\n|0.1\\n|0.1\\n|0.25\\n|0.32\\n|0.21\\n|0.41\\n|0.42\\n|0.8\\n|0.8\\n|0.8\\n|1\\n|1\\n|1\\n|1\\n|1\\n|1\\n|-\\n|Leu\\n|0.15\\n|0.05\\n|\\n|0\\n|0.03\\n|0.03\\n|0.28\\n|0.36\\n|0.2\\n|0.43\\n|0.51\\n|0.8\\n|0.8\\n|0.81\\n|1\\n|1\\n|1\\n|1\\n|1\\n|1.01\\n|-\\n|Ile\\n|0.15\\n|0.05\\n|0\\n|\\n|0.03\\n|0.03\\n|0.28\\n|0.36\\n|0.2\\n|0.43\\n|0.51\\n|0.8\\n|0.8\\n|0.81\\n|1\\n|1\\n|1\\n|1\\n|1\\n|1.01\\n|-\\n|Val\\n|0.2\\n|0.1\\n|0.05\\n|0.05\\n|\\n|0\\n|0.32\\n|0.4\\n|0.2\\n|0.4\\n|0.5\\n|0.8\\n|0.8\\n|0.81\\n|1\\n|1\\n|1\\n|1\\n|1\\n|1.02\\n|-\\n|Pro\\n|0.2\\n|0.1\\n|0.05\\n|0.05\\n|0\\n|\\n|0.32\\n|0.4\\n|0.2\\n|0.4\\n|0.5\\n|0.8\\n|0.8\\n|0.81\\n|1\\n|1\\n|1\\n|1\\n|1\\n|1.02\\n|-\\n|Tyr\\n|0.2\\n|0.22\\n|0.22\\n|0.22\\n|0.24\\n|0.24\\n|\\n|0.1\\n|0.13\\n|0.27\\n|0.36\\n|0.62\\n|0.61\\n|0.6\\n|0.8\\n|0.8\\n|0.81\\n|0.81\\n|0.8\\n|0.8\\n|-\\n|Trp\\n|0.21\\n|0.24\\n|0.25\\n|0.25\\n|0.27\\n|0.27\\n|0.05\\n|\\n|0.18\\n|0.3\\n|0.39\\n|0.63\\n|0.63\\n|0.61\\n|0.81\\n|0.81\\n|0.81\\n|0.81\\n|0.81\\n|0.8\\n|-\\n|Cys\\n|0.28\\n|0.22\\n|0.21\\n|0.21\\n|0.2\\n|0.2\\n|0.25\\n|0.35\\n|\\n|0.25\\n|0.31\\n|0.6\\n|0.6\\n|0.62\\n|0.81\\n|0.81\\n|0.8\\n|0.8\\n|0.81\\n|0.82\\n|-\\n|Ala\\n|0.5\\n|0.45\\n|0.43\\n|0.43\\n|0.41\\n|0.41\\n|0.4\\n|0.49\\n|0.22\\n|\\n|0.1\\n|0.4\\n|0.41\\n|0.47\\n|0.63\\n|0.63\\n|0.62\\n|0.62\\n|0.63\\n|0.67\\n|-\\n|Gly\\n|0.61\\n|0.56\\n|0.54\\n|0.54\\n|0.52\\n|0.52\\n|0.5\\n|0.58\\n|0.34\\n|0.1\\n|\\n|0.32\\n|0.34\\n|0.42\\n|0.56\\n|0.56\\n|0.54\\n|0.54\\n|0.56\\n|0.61\\n|-\\n|Ser\\n|0.81\\n|0.8\\n|0.8\\n|0.8\\n|0.8\\n|0.8\\n|0.62\\n|0.63\\n|0.6\\n|0.4\\n|0.3\\n|\\n|0.03\\n|0.1\\n|0.21\\n|0.21\\n|0.2\\n|0.2\\n|0.21\\n|0.24\\n|-\\n|Thr\\n|0.81\\n|0.8\\n|0.8\\n|0.8\\n|0.8\\n|0.8\\n|0.61\\n|0.63\\n|0.6\\n|0.4\\n|0.31\\n|0.03\\n|\\n|0.08\\n|0.21\\n|0.21\\n|0.2\\n|0.2\\n|0.21\\n|0.22\\n|-\\n|His\\n|0.8\\n|0.8\\n|1\\n|1\\n|0.8\\n|0.8\\n|0.6\\n|0.61\\n|0.61\\n|0.42\\n|0.34\\n|0.1\\n|0.08\\n|\\n|0.2\\n|0.2\\n|0.21\\n|0.21\\n|0.2\\n|0.2\\n|-\\n|Glu\\n|1\\n|1\\n|1\\n|1\\n|1\\n|1\\n|0.8\\n|0.81\\n|0.8\\n|0.61\\n|0.52\\n|0.22\\n|0.21\\n|0.2\\n|\\n|0\\n|0.03\\n|0.03\\n|0\\n|0.05\\n|-\\n|Gln\\n|1\\n|1\\n|1\\n|1\\n|1\\n|1\\n|0.8\\n|0.81\\n|0.8\\n|0.61\\n|0.52\\n|0.22\\n|0.21\\n|0.2\\n|0\\n|\\n|0.03\\n|0.03\\n|0\\n|0.05\\n|-\\n|Asp\\n|1\\n|1\\n|1\\n|1\\n|1\\n|1\\n|0.81\\n|0.81\\n|0.8\\n|0.61\\n|0.51\\n|0.21\\n|0.2\\n|0.21\\n|0.03\\n|0.03\\n|\\n|0\\n|0.03\\n|0.08\\n|-\\n|Asn\\n|1\\n|1\\n|1\\n|1\\n|1\\n|1\\n|0.81\\n|0.81\\n|0.8\\n|0.61\\n|0.51\\n|0.21\\n|0.2\\n|0.21\\n|0.03\\n|0.03\\n|0\\n|\\n|0.03\\n|0.08\\n|-\\n|Lys\\n|1\\n|1\\n|1\\n|1\\n|1\\n|1\\n|0.8\\n|0.81\\n|0.8\\n|0.61\\n|0.52\\n|0.22\\n|0.21\\n|0.2\\n|0\\n|0\\n|0.03\\n|0.03\\n|\\n|0.05\\n|-\\n|Arg\\n|1\\n|1\\n|1\\n|1\\n|1.01\\n|1.01\\n|0.8\\n|0.8\\n|0.81\\n|0.62\\n|0.53\\n|0.24\\n|0.22\\n|0.2\\n|0.05\\n|0.05\\n|0.08\\n|0.08\\n|0.05\\n|\\n|}\\n|}')

### grantham( = '{| class="wikitable"\\n|Arg\\n|Leu\\n|Pro\\n|Thr\\n|Ala\\n|Val\\n|Gly\\n|Ile\\n|Phe\\n|Tyr\\n|Cys\\n|His\\n|Gln\\n|Asn\\n|Lys\\n|Asp\\n|Glu\\n|Met\\n|Trp\\n|\\n|-\\n|110\\n|145\\n|74\\n|58\\n|99\\n|124\\n|56\\n|142\\n|155\\n|144\\n|112\\n|89\\n|68\\n|46\\n|121\\n|65\\n|80\\n|135\\n|177\\n|Ser\\n|-\\n|\\n|102\\n|103\\n|71\\n|112\\n|96\\n|125\\n|97\\n|97\\n|77\\n|180\\n|29\\n|43\\n|86\\n|26\\n|96\\n|54\\n|91\\n|101\\n|Arg\\n|-\\n|\\n|\\n|98\\n|92\\n|96\\n|32\\n|138\\n|5\\n|22\\n|36\\n|198\\n|99\\n|113\\n|153\\n|107\\n|172\\n|138\\n|15\\n|61\\n|Leu\\n|-\\n|\\n|\\n|\\n|38\\n|27\\n|68\\n|42\\n|95\\n|114\\n|110\\n|169\\n|77\\n|76\\n|91\\n|103\\n|108\\n|93\\n|87\\n|147\\n|Pro\\n|-\\n|\\n|\\n|\\n|\\n|58\\n|69\\n|59\\n|89\\n|103\\n|92\\n|149\\n|47\\n|42\\n|65\\n|78\\n|85\\n|65\\n|81\\n|128\\n|Thr\\n|-\\n|\\n|\\n|\\n|\\n|\\n|64\\n|60\\n|94\\n|113\\n|112\\n|195\\n|86\\n|91\\n|111\\n|106\\n|126\\n|107\\n|84\\n|148\\n|Ala\\n|-\\n|\\n|\\n|\\n|\\n|\\n|\\n|109\\n|29\\n|50\\n|55\\n|192\\n|84\\n|96\\n|133\\n|97\\n|152\\n|121\\n|21\\n|88\\n|Val\\n|-\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|135\\n|153\\n|147\\n|159\\n|98\\n|87\\n|80\\n|127\\n|94\\n|98\\n|127\\n|184\\n|Gly\\n|-\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|21\\n|33\\n|198\\n|94\\n|109\\n|149\\n|102\\n|168\\n|134\\n|10\\n|61\\n|Ile\\n|-\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|22\\n|205\\n|100\\n|116\\n|158\\n|102\\n|177\\n|140\\n|28\\n|40\\n|Phe\\n|-\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|194\\n|83\\n|99\\n|143\\n|85\\n|160\\n|122\\n|36\\n|37\\n|Tyr\\n|-\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|174\\n|154\\n|139\\n|202\\n|154\\n|170\\n|196\\n|215\\n|Cys\\n|-\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|24\\n|68\\n|32\\n|81\\n|40\\n|87\\n|115\\n|His\\n|-\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|46\\n|53\\n|61\\n|29\\n|101\\n|130\\n|Gln\\n|-\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|94\\n|23\\n|42\\n|142\\n|174\\n|Asn\\n|-\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|101\\n|56\\n|95\\n|110\\n|Lys\\n|-\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|45\\n|160\\n|181\\n|Asp\\n|-\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|126\\n|152\\n|Glu\\n|-\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|67\\n|Met\\n|} \\n|}')

### miyata( = '{| class="wikitable"\\n|Cys\\n|Pro\\n|Ala\\n|Gly\\n|Ser\\n|Thr\\n|Gln\\n|Glu\\n|Asn\\n|Asp\\n|His\\n|Lys\\n|Arg\\n|Val\\n|Leu\\n|Ile\\n|Met\\n|Phe\\n|Tyr\\n|Trp\\n|\\n|-\\n|\\n|1.33\\n|1.39\\n|2.22\\n|2.84\\n|1.45\\n|2.48\\n|3.26\\n|2.83\\n|3.48\\n|2.56\\n|3.27\\n|3.06\\n|0.86\\n|1.65\\n|1.63\\n|1.46\\n|2.24\\n|2.38\\n|3.34\\n|Cys\\n|-\\n|\\n|\\n|0.06\\n|0.97\\n|0.56\\n|0.87\\n|1.92\\n|2.48\\n|1.8\\n|2.4\\n|2.15\\n|2.94\\n|2.9\\n|1.79\\n|2.7\\n|2.62\\n|2.36\\n|3.17\\n|3.12\\n|4.17\\n|Pro\\n|-\\n|\\n|\\n|\\n|0.91\\n|0.51\\n|0.9\\n|1.92\\n|2.46\\n|1.78\\n|2.37\\n|2.17\\n|2.96\\n|2.92\\n|1.85\\n|2.76\\n|2.69\\n|2.42\\n|3.23\\n|3.18\\n|4.23\\n|Ala\\n|-\\n|\\n|\\n|\\n|\\n|0.85\\n|1.7\\n|2.48\\n|2.78\\n|1.96\\n|2.37\\n|2.78\\n|3.54\\n|3.58\\n|2.76\\n|3.67\\n|3.6\\n|3.34\\n|4.14\\n|4.08\\n|5.13\\n|Gly\\n|-\\n|\\n|\\n|\\n|\\n|\\n|0.89\\n|1.65\\n|2.06\\n|1.31\\n|1.87\\n|1.94\\n|2.71\\n|2.74\\n|2.15\\n|3.04\\n|2.95\\n|2.67\\n|3.45\\n|3.33\\n|4.38\\n|Ser\\n|-\\n|\\n|\\n|\\n|\\n|\\n|\\n|1.12\\n|1.83\\n|1.4\\n|2.05\\n|1.32\\n|2.1\\n|2.03\\n|1.42\\n|2.25\\n|2.14\\n|1.86\\n|2.6\\n|2.45\\n|3.5\\n|Thr\\n|-\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|0.84\\n|0.99\\n|1.47\\n|0.32\\n|1.06\\n|1.13\\n|2.13\\n|2.7\\n|2.57\\n|2.3\\n|2.81\\n|2.48\\n|3.42\\n|Gln\\n|-\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|0.85\\n|0.9\\n|0.96\\n|1.14\\n|1.45\\n|2.97\\n|3.53\\n|3.39\\n|3.13\\n|3.59\\n|3.22\\n|4.08\\n|Glu\\n|-\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|0.65\\n|1.29\\n|1.84\\n|2.04\\n|2.76\\n|3.49\\n|3.37\\n|3.08\\n|3.7\\n|3.42\\n|4.39\\n|Asn\\n|-\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|1.72\\n|2.05\\n|2.34\\n|3.4\\n|4.1\\n|3.98\\n|3.69\\n|4.27\\n|3.95\\n|4.88\\n|Asp\\n|-\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|0.79\\n|0.82\\n|2.11\\n|2.59\\n|2.45\\n|2.19\\n|2.63\\n|2.27\\n|3.16\\n|His\\n|-\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|0.4\\n|2.7\\n|2.98\\n|2.84\\n|2.63\\n|2.85\\n|2.42\\n|3.11\\n|Lys\\n|-\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|2.43\\n|2.62\\n|2.49\\n|2.29\\n|2.47\\n|2.02\\n|2.72\\n|Arg\\n|-\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|0.91\\n|0.85\\n|0.62\\n|1.43\\n|1.52\\n|2.51\\n|Val\\n|-\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|0.14\\n|0.41\\n|0.63\\n|0.94\\n|1.73\\n|Leu\\n|-\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|0.29\\n|0.61\\n|0.86\\n|1.72\\n|Ile\\n|-\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|0.82\\n|0.93\\n|1.89\\n|Met\\n|-\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|0.48\\n|1.11\\n|Phe\\n|-\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|1.06\\n|Tyr\\n|-\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|Trp\\n|} \\n|-\\n|\\n|}')

### sneath( = '{| class="wikitable"\\n|\\n|Leu\\n|Ile\\n|Val\\n|Gly\\n|Ala\\n|Pro\\n|Gln\\n|Asn\\n|Met\\n|Thr\\n|Ser\\n|Cys\\n|Glu\\n|Asp\\n|Lys\\n|Arg\\n|Tyr\\n|Phe\\n|Trp\\n|-\\n|Ile\\n|5\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|-\\n|Val\\n|9\\n|7\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|-\\n|Gly\\n|24\\n|25\\n|19\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|-\\n|Ala\\n|15\\n|17\\n|12\\n|9\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|-\\n|Pro\\n|23\\n|24\\n|20\\n|17\\n|16\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|-\\n|Glu\\n|22\\n|24\\n|25\\n|32\\n|26\\n|33\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|-\\n|Asn\\n|20\\n|23\\n|23\\n|26\\n|25\\n|31\\n|10\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|-\\n|Met\\n|20\\n|22\\n|23\\n|34\\n|25\\n|31\\n|13\\n|21\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|-\\n|Thr\\n|23\\n|21\\n|17\\n|20\\n|20\\n|25\\n|24\\n|19\\n|25\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|-\\n|Ser\\n|23\\n|25\\n|20\\n|19\\n|16\\n|24\\n|21\\n|15\\n|22\\n|12\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|-\\n|Cys\\n|24\\n|26\\n|21\\n|21\\n|13\\n|25\\n|22\\n|19\\n|17\\n|19\\n|13\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|-\\n|Glu\\n|30\\n|31\\n|31\\n|37\\n|34\\n|43\\n|14\\n|19\\n|26\\n|34\\n|29\\n|33\\n|\\n|\\n|\\n|\\n|\\n|\\n|\\n|-\\n|Asp\\n|25\\n|28\\n|28\\n|33\\n|30\\n|40\\n|22\\n|14\\n|31\\n|29\\n|25\\n|28\\n|7\\n|\\n|\\n|\\n|\\n|\\n|\\n|-\\n|Lys\\n|23\\n|24\\n|26\\n|31\\n|26\\n|31\\n|21\\n|27\\n|24\\n|34\\n|31\\n|32\\n|26\\n|34\\n|\\n|\\n|\\n|\\n|\\n|-\\n|Arg\\n|33\\n|34\\n|36\\n|43\\n|37\\n|43\\n|23\\n|31\\n|28\\n|38\\n|37\\n|36\\n|31\\n|39\\n|14\\n|\\n|\\n|\\n|\\n|-\\n|Tyr\\n|30\\n|34\\n|36\\n|36\\n|34\\n|37\\n|29\\n|28\\n|32\\n|32\\n|29\\n|34\\n|34\\n|34\\n|34\\n|36\\n|\\n|\\n|\\n|-\\n|Phe\\n|19\\n|22\\n|26\\n|29\\n|26\\n|27\\n|24\\n|24\\n|24\\n|28\\n|25\\n|29\\n|35\\n|35\\n|28\\n|34\\n|13\\n|\\n|\\n|-\\n|Trp\\n|30\\n|34\\n|37\\n|39\\n|36\\n|37\\n|31\\n|32\\n|31\\n|38\\n|35\\n|37\\n|43\\n|45\\n|34\\n|36\\n|21\\n|13\\n|\\n|-\\n|His\\n|25\\n|28\\n|31\\n|34\\n|29\\n|36\\n|27\\n|24\\n|30\\n|34\\n|28\\n|31\\n|27\\n|35\\n|27\\n|31\\n|23\\n|18\\n|25\\n|} \\n|}')

### classmethod wiki_to_csv(table)

### yampolsky( = '\\n{| class="wikitable"\\n|\\n|Cys\\n|Ser\\n|Thr\\n|Pro\\n|Ala\\n|Gly\\n|Asn\\n|Asp\\n|Glu\\n|Gln\\n|His\\n|Arg\\n|Lys\\n|Met\\n|Ile\\n|Leu\\n|Val\\n|Phe\\n|Tyr\\n|Trp\\n|Ex<sub>src</sub>\\n|-\\n|Cys\\n|.\\n|258\\n|121\\n|201\\n|334\\n|288\\n|109\\n|109\\n|270\\n|383\\n|258\\n|306\\n|252\\n|169\\n|109\\n|347\\n|89\\n|349\\n|349\\n|139\\n|280\\n|-\\n|Ser\\n|373\\n|.\\n|481\\n|249\\n|490\\n|418\\n|390\\n|314\\n|343\\n|352\\n|353\\n|363\\n|275\\n|321\\n|270\\n|295\\n|358\\n|334\\n|294\\n|160\\n|351\\n|-\\n|Thr\\n|325\\n|408\\n|.\\n|164\\n|402\\n|332\\n|240\\n|190\\n|212\\n|308\\n|246\\n|299\\n|256\\n|152\\n|198\\n|271\\n|362\\n|273\\n|260\\n|66\\n|287\\n|-\\n|Pro\\n|345\\n|392\\n|286\\n|.\\n|454\\n|404\\n|352\\n|254\\n|346\\n|384\\n|369\\n|254\\n|231\\n|257\\n|204\\n|258\\n|421\\n|339\\n|298\\n|305\\n|335\\n|-\\n|Ala\\n|393\\n|384\\n|312\\n|243\\n|.\\n|387\\n|430\\n|193\\n|275\\n|320\\n|301\\n|295\\n|225\\n|549\\n|245\\n|313\\n|319\\n|305\\n|286\\n|165\\n|312\\n|-\\n|Gly\\n|267\\n|304\\n|187\\n|140\\n|369\\n|.\\n|210\\n|188\\n|206\\n|272\\n|235\\n|178\\n|219\\n|197\\n|110\\n|193\\n|208\\n|168\\n|188\\n|173\\n|228\\n|-\\n|Asn\\n|234\\n|355\\n|329\\n|275\\n|400\\n|391\\n|.\\n|208\\n|257\\n|298\\n|248\\n|252\\n|183\\n|236\\n|184\\n|233\\n|233\\n|210\\n|251\\n|120\\n|272\\n|-\\n|Asp\\n|285\\n|275\\n|245\\n|220\\n|293\\n|264\\n|201\\n|.\\n|344\\n|263\\n|298\\n|252\\n|208\\n|245\\n|299\\n|236\\n|175\\n|233\\n|227\\n|103\\n|258\\n|-\\n|Glu\\n|332\\n|355\\n|292\\n|216\\n|520\\n|407\\n|258\\n|533\\n|.\\n|341\\n|380\\n|279\\n|323\\n|219\\n|450\\n|321\\n|351\\n|342\\n|348\\n|145\\n|363\\n|-\\n|Gln\\n|383\\n|443\\n|361\\n|212\\n|499\\n|406\\n|338\\n|68\\n|439\\n|.\\n|396\\n|366\\n|354\\n|504\\n|467\\n|391\\n|603\\n|383\\n|361\\n|159\\n|386\\n|-\\n|His\\n|331\\n|365\\n|205\\n|220\\n|462\\n|370\\n|225\\n|141\\n|319\\n|301\\n|.\\n|275\\n|332\\n|315\\n|205\\n|364\\n|255\\n|328\\n|260\\n|72\\n|303\\n|-\\n|Arg\\n|225\\n|270\\n|199\\n|145\\n|459\\n|251\\n|67\\n|124\\n|250\\n|288\\n|263\\n|.\\n|306\\n|68\\n|139\\n|242\\n|189\\n|213\\n|272\\n|63\\n|259\\n|-\\n|Lys\\n|331\\n|376\\n|476\\n|252\\n|600\\n|492\\n|457\\n|465\\n|272\\n|441\\n|362\\n|440\\n|.\\n|414\\n|491\\n|301\\n|487\\n|360\\n|343\\n|218\\n|409\\n|-\\n|Met\\n|347\\n|353\\n|261\\n|85\\n|357\\n|218\\n|544\\n|392\\n|287\\n|394\\n|278\\n|112\\n|135\\n|.\\n|612\\n|513\\n|354\\n|330\\n|308\\n|633\\n|307\\n|-\\n|Ile\\n|362\\n|196\\n|193\\n|145\\n|326\\n|160\\n|172\\n|27\\n|197\\n|191\\n|221\\n|124\\n|121\\n|279\\n|.\\n|417\\n|494\\n|331\\n|323\\n|73\\n|252\\n|-\\n|Leu\\n|366\\n|212\\n|165\\n|146\\n|343\\n|201\\n|162\\n|112\\n|199\\n|250\\n|288\\n|185\\n|171\\n|367\\n|301\\n|.\\n|275\\n|336\\n|295\\n|152\\n|248\\n|-\\n|Val\\n|382\\n|326\\n|398\\n|201\\n|389\\n|269\\n|108\\n|228\\n|192\\n|280\\n|253\\n|190\\n|197\\n|562\\n|537\\n|333\\n|.\\n|207\\n|209\\n|286\\n|277\\n|-\\n|Phe\\n|176\\n|152\\n|257\\n|112\\n|236\\n|94\\n|136\\n|90\\n|62\\n|216\\n|237\\n|122\\n|85\\n|255\\n|181\\n|296\\n|291\\n|.\\n|332\\n|232\\n|193\\n|-\\n|Tyr\\n|142\\n|173\\n|.\\n|194\\n|402\\n|357\\n|129\\n|87\\n|176\\n|369\\n|197\\n|340\\n|171\\n|392\\n|.\\n|362\\n|.\\n|360\\n|.\\n|303\\n|258\\n|-\\n|Trp\\n|137\\n|92\\n|17\\n|66\\n|63\\n|162\\n|.\\n|.\\n|65\\n|61\\n|239\\n|103\\n|54\\n|110\\n|.\\n|177\\n|110\\n|364\\n|281\\n|.\\n|142\\n|-\\n|Ex<sub>dest</sub>\\n|315\\n|311\\n|293\\n|192\\n|411\\n|321\\n|258\\n|225\\n|262\\n|305\\n|290\\n|255\\n|225\\n|314\\n|293\\n|307\\n|305\\n|294\\n|279\\n|172\\n|291\\n|} \\n|}')
## protein.core module


#### class protein.core.ProteinCore(gene_name='', uniprot='', uniprot_name='', sequence='', organism=None, taxid=None, \*\*other)
Bases: `object`

This is a lightweight version of Protein that is intended to run off pre parsed pickles.
It forms the base of Protein. This does no protein analyses.
It has IO powers though .dump/.gdump saves an instance .load/.gload loads and can work as a class method if the filename is provided as an argument.
The gzipped forms (.gdump and .gload) are about 1/3 the size. 50 KB.


### property ExAC_type()

### \__init__(gene_name='', uniprot='', uniprot_name='', sequence='', organism=None, taxid=None, \*\*other)
Initialize self.  See help(type(self)) for accurate signature.


### complete()
Make sure that all subthreads are complete. Not used for Core!


### dump(file=None)

### exists(file=None)
Method to check if file exists already.
Actually loads it sneakily!
:return:


### gdump(file=None)

### get_species_for_uniprot()

### gload(file=None)

### load(file=None)

### log(text)
Logging is primarily for protein_full
:param text:
:return:


### settings( = <protein.settings_handler.GlobalSettings object>)

### version( = 1.0)

#### class protein.core.Structure(id, description, x: int, y: int, code, type='rcsb', chain='A', offset: int = 0, coordinates=None, extra=None, url='')
Bases: `object`

No longer a namedtuple.
Stores the structural data for easy use by FeatureViewer and co. Can be converted to StructureAnalyser
type = rcsb | swissmodel | homologue


### \__init__(id, description, x: int, y: int, code, type='rcsb', chain='A', offset: int = 0, coordinates=None, extra=None, url='')
Stores the structural data for easy use by FeatureViewer and co. Can be converted to StructureAnalyser
type = rcsb | swissmodel | homologue | www | local


### get_coordinates()

### includes(position, offset=0)
Generally there should not be an offset as x and y are from Uniprot data so they are already fixed!
:param position:
:param offset:
:return:


### lookup_ligand()

### lookup_resolution()

### lookup_sifts()
SIFTS data. for PDBe query see elsewhere.
There are four start/stop pairs that need to be compared to get a good idea of a protein.
For a lengthy discussion see [https://blog.matteoferla.com/2019/09/pdb-numbering-rollercoaster.html](https://blog.matteoferla.com/2019/09/pdb-numbering-rollercoaster.html)
Also for a good list of corner case models see [https://proteopedia.org/wiki/index.php/Unusual_sequence_numbering](https://proteopedia.org/wiki/index.php/Unusual_sequence_numbering)
:return: self


### settings( = <protein.settings_handler.GlobalSettings object>)

### to_dict()

#### class protein.core.Variant()
Bases: `tuple`

Stores the gnomAD data for easy use by FeatureViewer and co. Can be converted to Mutation.


### property description()
Alias for field number 4


### property homozygous()
Alias for field number 5


### property id()
Alias for field number 0


### property impact()
Alias for field number 3


### property x()
Alias for field number 1


### property y()
Alias for field number 2

## protein.generate_AA_MCS module

## protein.metadata_from_PDBe module


#### class protein.metadata_from_PDBe.PDBMeta(entry)
Bases: `object`

Query the PDBe for what the code parts are.
Herein by chain the chain letter is meant, while the data of the chain is called entity… terrible.


### \__init__(entry)
Initialize self.  See help(type(self)) for accurate signature.


### describe()

### get_data_by_chain(chain=None)

### get_nonproteins()

### get_proteins()

### get_range_by_chain(chain=None)

### get_range_by_entity(entity)

### is_boring_ligand(entity)

### is_peptide(entity)

### wordy_describe(delimiter=' + ')

### wordy_describe_entity(entity)
## protein.mutation module


#### class protein.mutation.Mutation(mutation=None)
Bases: `object`

Stores the mutation. Not to be confused with the namedtuple Variant, which stores gnomAD mutations.
>>> Mutation(‘p.M12D’)
>>> Mutation(‘M12D’)
>>> Mutation(gnomAD_variant_instance)
This class does not do analyses with Protein, but ProteinAnalysis do. Here however, wordy conversions happen.


### \__init__(mutation=None)
Initialize self.  See help(type(self)) for accurate signature.


### aa_list( = ('Q', 'W', 'E', 'R', 'T', 'Y', 'I', 'P', 'A', 'S', 'D', 'F', 'G', 'H', 'K', 'L', 'C', 'V', 'N', 'M', '\*'))

### property exposure_effect()

### classmethod long_name(letter)
Single amino acid letter to a short string with the name spelt out three ways.
:param letter: 1 AA letter
:type letter: str
:return: str


### names( = (('A', 'Ala', 'Alanine'), ('B', 'Asx', 'Aspartate/asparagine'), ('C', 'Cys', 'Cysteine'), ('D', 'Asp', 'Aspartate'), ('E', 'Glu', 'Glutamate'), ('F', 'Phe', 'Phenylanine'), ('G', 'Gly', 'Glycine'), ('H', 'His', 'Histidine'), ('I', 'Ile', 'Isoleucine'), ('K', 'Lys', 'Lysine'), ('L', 'Leu', 'Leucine'), ('M', 'Met', 'Methionine'), ('N', 'Asn', 'Asparagine'), ('P', 'Pro', 'Proline'), ('Q', 'Gln', 'Glutamine'), ('R', 'Arg', 'Arginine'), ('S', 'Ser', 'Serine'), ('T', 'Thr', 'Threonine'), ('U', 'Sel', 'Selenocysteine'), ('V', 'Val', 'Valine'), ('W', 'Trp', 'Tryptophan'), ('X', 'Xaa', 'Any'), ('Y', 'Tyr', 'Tyrosine'), ('Z', 'Glx', 'Glutamate/glutamine'), ('\*', 'Stop', 'Stop')))

### parse_mutation(mutation)
## protein.protein_analysis module


#### class protein.protein_analysis.ProteinAnalyser(\*args, \*\*kwargs)
Bases: `protein.core.ProteinCore`


### \__init__(\*args, \*\*kwargs)
Initialize self.  See help(type(self)) for accurate signature.


### analyse_structure()

### check_elm()

### check_mutation()

### property elmdata()

### get_best_model()
This currently just gets the first PDB. It ought to check what is the best.
:return:


### get_features_at_position(position=None)

* **Parameters**

    **position** – mutation, str or position



* **Returns**

    list of gnomAD mutations, which are dictionary e.g. {‘id’: ‘gnomAD_19_19_rs562294556’, ‘description’: ‘R19Q (rs562294556)’, ‘x’: 19, ‘y’: 19, ‘impact’: ‘MODERATE’}



### get_features_near_position(position=None, wobble=10)

### get_gnomAD_near_position(position=None, wobble=5)

* **Parameters**

    
    * **position** – mutation, str or position


    * **wobble** – int, number of residues before and after.



* **Returns**

    list of gnomAD mutations, which are dictionary e.g. {‘id’: ‘gnomAD_19_19_rs562294556’, ‘description’: ‘R19Q (rs562294556)’, ‘x’: 19, ‘y’: 19, ‘impact’: ‘MODERATE’}



### property mutation()

### mutation_discrepancy()

### predict_effect()
main entry point for analyses.
Do note that there is another class called StructureAnalyser which deals with the model specific details.
:return:


#### class protein.protein_analysis.StructureAnalyser(structure, mutation)
Bases: `object`


### \__init__(structure, mutation)

* **Parameters**

    
    * **structure** – a instance of Structure, a former namedtuple and is in core.py


    * **mutation** – a instance of mutation.



### get_structure_neighbours(threshhold: float = 3)

* **Parameters**

    **threshhold** – &Aring;ngstrom distance.



* **Returns**

    


### get_superficiality()

### is_core()
Uses half solvent exposure within PDB module.


### is_full_surface()
Uses half solvent exposure within PDB module.


### is_partial_surface()
Uses half solvent exposure within PDB module.


### property target_hse()
## protein.protein_manual module


#### class protein.protein_manual.ProteinManual(gene_name='', uniprot='', uniprot_name='', sequence='', organism=None, taxid=None, \*\*other)
Bases: `protein.core.ProteinCore`

ProteinCore but with manual data! For a very specific use.
It may have rusted into non functionality.


### add_manual_data()
## protein.settings_handler module


#### class protein.settings_handler.GlobalSettings(home_url='/')
Bases: `object`

This class is container for the paths, which are used by both Variant and Tracker classes.
Hence why in these two is the attribute .settings


### \__init__(home_url='/')
Initialize self.  See help(type(self)) for accurate signature.


### addresses( = ('ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz', 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/pdbaa.tar.gz', 'ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/functional_gene_constraint/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt', 'ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz', 'ftp://ftp.wwpdb.org/pub/pdb/derived_data/index/resolu.idx', 'https://swissmodel.expasy.org/repository/download/core_species/9606_meta.tar.gz', 'https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.exome_calling_intervals.sites.vcf.bgz', 'https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz', 'http://www.sbg.bio.ic.ac.uk/~missense3d/download/1052_hom_str.zip'))

### create_json_from_idx(infile, outfile)

### property data_folder()

### degunk()
Removes the zero sized files that may ahve arised from error or keyboard interuption.
:return:


### error_tolerant( = False)

### fetch( = True)

### get_folder_of(name)

### manual_task_note( = "'## Manual TASKS\\nRemember that user has manually downloaded _site_dataset.gz files from https://www.phosphosite.org/staticDownloads at phosphosite.'")

### missing_attribute_tolerant( = True)

### open(kind)

### retrieve_references(ask=True, refresh=False, issue='')

### startup(data_folder='data')

### subdirectory_names( = ('reference', 'temp', 'uniprot', 'pdbblast', 'pickle', 'binders', 'dictionary'))

### verbose( = False)

### wipe_html()
No longer needed.
:return:


#### class protein.settings_handler.Singleton()
Bases: `type`

There can only be one setting.

## protein.to_reimplement module

This chunk of code is for looking at the neighnbours
It will be implemented in to the code but not now.


#### protein.to_reimplement.junk()
## protein.unitests module


#### class protein.unitests.TestProteinCore(methodName='runTest')
Bases: `unittest.case.TestCase`


### test_error_tollerant()

### test_full_parse()

### test_parse()

### test_warn()
## Module contents

This module has several classes.
\* ProteinCore is used by the rest and provides the backbone. It can read and write itself (even in compressed form, see .gdump and .gload) but not generate protein data. for that there is
\* generate.ProteinGatherer, which parses data from various sources. generate.ProteomeGatherer starts everything up and parses the whole proteome.
\* there are some classes in generate too, but they don’t see the light of day. For those see generate._protein_gatherer
\* Mutation handles the mutation
\* the global_settings variable, declaired in .settings_handler has handles the config stuff.
\* ProteinAnalysis analyses a mutation. 
\* Variant and Structure are just namedtuples

The submodule generate has a method generate, which generates all the datasets for the human proteome. Although the files required are in settings_handler

The Mutation class is mostly used by ProteinAnalysis but itself holds wordy 

```
*
```

_effects attributes and uses a variable that was originally generated in _apriori_effect.py.

The script unitests.py does… unitests.

The script protein.aprior_effect generates the dictionary that is used to say what the apriori effect are. Namely, what amino acid is smaller etc.
This is already added to protein, so does not need to be redone.
>>> from protein.apriori_effect import Changedex
>>> pprint(Changedex().fill().to_dict())

To start everything going…
>>> from protein.generate.ProteomeGatherer
>>> ProteomeGatherer()
and wait for the missing file. There surely is one. generate.__init__ has a list of references.
