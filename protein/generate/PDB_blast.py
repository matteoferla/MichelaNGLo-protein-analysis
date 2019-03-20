__docs__="""
This part requires NCBI blast tools!
"""
import json
import os
import tarfile
from ..settings_handler import global_settings
from Bio.Blast import NCBIXML

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
        #load
        genes={}
        name='error'
        for line in open(os.path.join(global_settings.temp_folder, 'human.fa')):
            if '>' in line:
                name = line[1:].rstrip()
                genes[name]=''
            else:
                genes[name]+=line.rstrip()
        #dump
        folder = os.path.join(global_settings.temp_folder, 'fasta')
        os.mkdir(folder)
        for acc in genes:
            with open(os.path.join(folder,acc+'.fa'),'w') as w:
                w.write('>'+acc+'\n'+genes[acc]+'\n')
        return cls

    @classmethod
    def pdb_blaster(cls):
        return cls.full_blaster('blastpdb','pdbaa')

    @classmethod
    def self_blaster(cls):
        return cls.full_blaster('blastself', 'human')


    @classmethod
    def full_blaster(cls, outfolder_name, db):
        """
        Given the list of genes in in seqdex.json. do a blast against pdbaa from NCBI ftp.
        :return:
        """
        outfolder = os.path.join(global_settings.temp_folder, outfolder_name)
        os.mkdir(outfolder)
        infolder = os.path.join(global_settings.temp_folder, 'fasta')
        if not os.path.exists(infolder):
            cls.make_fastas()
        for infile in os.listdir(infolder):
            outfile = infile.replace('.fa','.xml')
            os.system('blastp -query {infile} -db {db} -outfmt 5 -num_threads 6 > {outfile}'.format(db=db,
                                                                                                    infile=os.path.join(infolder,infile),
                                                                                                     outfile=os.path.join(outfolder,outfile)))
        return cls

    @staticmethod
    def part_blaster(todo):
        """
        Like full blaster but for the fails.
        :param todo: todo is a set of ids that may have failed.
        :return:
        """
        dex = json.load(open('human_prot_seqdex.json','r'))
        for k in todo:
            file='blastpdb/'+k+'.fa'
            open(file,'w').write('>{i}\n{s}\n\n'.format(s=dex[k],i=k))
            os.system('blastp -query {infile} -db pdbaa -outfmt 5 -num_threads 6 > {outfile}'.format(infile=file, outfile=file.replace('.fa','_blastPDB.xml')))

    @staticmethod
    def parse(folder):
        for file in os.listdir(os.path.join(global_settings.temp_folder,folder)):
            if '_blast.xml' in file:
                print(file)
                try:
                    blast_record = NCBIXML.read(open(os.path.join(global_settings.temp_folder,folder,file)))
                    matches = []
                    for align in blast_record.alignments:
                        print(align.title)
                        for hsp in align.hsps:
                            if hsp.score > 100:
                                d = {'match': align.title[0:50],
                                     'match_score': hsp.score,
                                     'match_start': hsp.query_start,
                                     'match_length': hsp.align_length,
                                     'match_identity': hsp.identities / hsp.align_length}
                                d['formatted'] = {'x': hsp.query_start,
                                                  'y': hsp.align_length + hsp.query_start,
                                                  'description': align.title[0:20],
                                                  'id': 'blastpdb_'}
                                matches.append(d)
                except ValueError as err:
                    print('Value error: '+str(err)) ##why art thou so empty?

    @staticmethod
    def _test_describe():
        blast_record = NCBIXML.read(open(os.path.join(global_settings.temp_folder, 'blastpdb', 'S438966_blast.xml')))
        print(blast_record.alignments)
        print(blast_record.alignments[0].title)



