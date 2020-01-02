from .core import ProteinCore
import os, csv, markdown
from warnings import warn
from shutil import copyfile

class ProteinManual(ProteinCore):
    """
    ProteinCore but with manual data! For a very specific use.
    It may have rusted into non functionality.
    """

    def add_manual_data(self):
        man = {x['Gene']: x for x in csv.DictReader(open(os.path.join(self.settings.manual_folder, 'manual.csv')))}
        ### PDB
        pdb_candidate = os.path.join(self.settings.manual_folder, self.uniprot + '.pdb')
        if os.path.isfile(pdb_candidate):
            pdb_man_file = pdb_candidate
        elif self.uniprot in man and man[self.uniprot]['PDB']:
            pdb_candidate = os.path.join(self.settings.manual_folder, man[self.uniprot]['PDB'].lstrip().rstrip())
            if os.path.isfile(pdb_candidate):
                pdb_man_file = pdb_candidate
            else:  # csv is wrong
                pdb_man_file = None
                warn('Cannot find structure file, claimed by manual.csv: '+pdb_candidate)
        else:
            pdb_man_file = None
        if pdb_man_file:
            self.pdb_file = os.path.join(self.settings.page_folder, self.uniprot + '.pdb')
            self.pdb_resi = 0  # self.resi - int(man[self.gene]['PDB_offset'])
            if not os.path.isfile(self.pdb_file):  # make the html copy
                copyfile(pdb_man_file, self.pdb_file)
        ### MD
        txt_candidate = os.path.join(self.settings.manual_folder, self.uniprot + '.md')
        if os.path.isfile(txt_candidate):
            txt = txt_candidate
        elif self.gene_name in man and man[self.uniprot]['Text']:
            txt_candidate = os.path.join(self.settings.manual_folder, man[self.uniprot]['Text'])
            if os.path.isfile(txt_candidate):
                txt = txt_candidate
            else:
                warn('Cannot find text file, ' + txt_candidate + ' claimed by manual.csv.')
                txt = None
        else:
            txt = None
        if txt:
            self.user_text = markdown.markdown(open(txt).read())
        # log
        if self.user_text and self.pdb_file:
            self.log('Manual data added.')
        elif self.user_text and not self.pdb_file:
            self.log('Manual data. No structure.')
        else:  # this is unnedded except for legacy pickled files.
            self.user_text = ''
            self.pdb_file = ''
        # done
        return self
