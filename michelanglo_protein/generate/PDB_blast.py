__docs__ = """
This part requires NCBI blast tools!
"""
import json
import os
import tarfile
from ..settings_handler import global_settings
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline #the worlds most pointless wrapper...

from warnings import warn

class Blaster:
    """
    This is just a container.
    """

    @staticmethod
    def extract_db():
        file = os.path.join(global_settings.reference_folder, 'pdbaa.tar.gz')
        tar = tarfile.open(file)
        tar.extractall(global_settings.temp_folder)
        tar.close()

    @classmethod
    def make_fastas(cls):
        # load
        genes = {}
        name = 'error'
        for line in open(os.path.join(global_settings.temp_folder, 'human.fa')):
            if '>' in line:
                name = line[1:].rstrip()
                genes[name] = ''
            else:
                genes[name] += line.rstrip()
        # dump
        folder = os.path.join(global_settings.temp_folder, 'fasta')
        os.mkdir(folder)
        for acc in genes:
            with open(os.path.join(folder, acc + '.fa'), 'w') as w:
                w.write('>' + acc + '\n' + genes[acc] + '\n')
        return cls

    @classmethod
    def pdb_blaster(cls):
        return cls.full_blaster('blastpdb', os.path.join(global_settings.temp_folder,'pdbaa'))

    @classmethod
    def self_blaster(cls):
        return cls.full_blaster('blastself', os.path.join(global_settings.temp_folder,'human'))

    @classmethod
    def full_blaster(cls, outfolder_name, db):
        """
        Given the list of genes in in seqdex.json. do a blast against pdbaa from NCBI ftp.
        :return:
        """
        outfolder = os.path.join(global_settings.temp_folder, outfolder_name)
        if not os.path.exists(outfolder):
            os.mkdir(outfolder)
        infolder = os.path.join(global_settings.temp_folder, 'fasta')
        if not os.path.exists(infolder):
            cls.make_fastas()
        for infile in os.listdir(infolder):
            outfile = infile.replace('.fa', '.xml')
            if global_settings.verbose:
                print(infile)
            NcbiblastpCommandline(cmd='C:\\Program Files\\NCBI\\blast-2.8.1+\\bin\\blastp.exe',
                                  query=os.path.join(infolder, infile),
                                  db=db,
                                  evalue=0.001,
                                  outfmt = 5,
                                  out = os.path.join(outfolder, outfile))()
        return cls

    @staticmethod
    def part_blaster(todo):
        raise DeprecationWarning
        """
        Like full blaster but for the fails.
        :param todo: todo is a set of ids that may have failed.
        :return:
        """
        dex = json.load(open('human_prot_seqdex.json', 'r'))
        for k in todo:
            file = 'blastpdb/' + k + '.fa'
            open(file, 'w').write('>{i}\n{s}\n\n'.format(s=dex[k], i=k))
            os.system('blastp -query {infile} -db pdbaa -outfmt 5 -num_threads 6 > {outfile}'.format(infile=file, outfile=file.replace('.fa', '_blastPDB.xml')))

    @staticmethod
    def parse(infolder,outfolder):
        if not os.path.exists(os.path.join(global_settings.temp_folder, outfolder)):
            os.mkdir(os.path.join(global_settings.temp_folder, outfolder))
        for file in os.listdir(os.path.join(global_settings.temp_folder, infolder)):
            if '.xml' in file:
                try:
                    blast_record = NCBIXML.read(open(os.path.join(global_settings.temp_folder, infolder, file)))
                    matches = []
                    for align in blast_record.alignments:
                        for hsp in align.hsps:
                            if hsp.score > 100:
                                pdb = align.title.split('|')[3]
                                chain = align.title.split('|')[4][0]
                                d = {'x': int(hsp.query_start),
                                     'y': int(hsp.align_length + hsp.query_start),
                                     'description': align.hit_def.split('&gt;')[0],
                                     'id': 'blastpdb_{p}_{x}_{y}_{c}'.format(p=pdb, c=chain, x=hsp.query_start, y=hsp.align_length + hsp.query_start),
                                     'chain': chain,
                                     'url': pdb,
                                     'offset': hsp.sbjct_start - hsp.query_start,
                                     'extra': {'match': align.title[0:50],
                                                  'match_score': hsp.score,
                                                  'match_start': hsp.query_start,
                                                  'match_length': hsp.align_length,
                                                  'match_identity': hsp.identities / hsp.align_length}}
                                matches.append(d)
                    with open(os.path.join(global_settings.temp_folder, outfolder, file.replace('.xml','.json')),'w') as w:
                        json.dump(matches,w)
                except ValueError as err:
                    warn(f'Value error for {file}: {err}')  ##why art thou so empty?

    @staticmethod
    def _test_describe():
        blast_record = NCBIXML.read(open(os.path.join(global_settings.temp_folder, 'blastpdb', 'S438966_blast.xml')))
        print(blast_record.alignments)
        print(blast_record.alignments[0].title)
