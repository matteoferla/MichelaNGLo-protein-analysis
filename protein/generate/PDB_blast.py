__docs__="""
This part requires NCBI blast tools
"""
import json
import os

def full_blaster():
    """
    Given the list of genes in in seqdex.json. do a blast against pdbaa from NCBI ftp.
    :return:
    """
    os.mkdir('genes')
    dex = json.load(open('human_prot_seqdex.json','r'))
    done = json.load(open('done.json','r'))
    for k in dex:
        if k not in done:
            file='genes/'+k+'.fa'
            open(file,'w').write('>{i}\n{s}\n\n'.format(s=dex[k],i=k))
            os.system('blastp -query {infile} -db pdbaa -outfmt 5 -num_threads 6 > {outfile}'.format(infile=file, outfile=file.replace('.fa','_blastPDB.xml')))

def part_blaster(todo):
    """
    Like full blaster but for the fails.
    :param todo: todo is a set of ids that may have failed.
    :return:
    """
    dex = json.load(open('human_prot_seqdex.json','r'))
    for k in todo:
        file='genes/'+k+'.fa'
        open(file,'w').write('>{i}\n{s}\n\n'.format(s=dex[k],i=k))
        os.system('blastp -query {infile} -db pdbaa -outfmt 5 -num_threads 6 > {outfile}'.format(infile=file, outfile=file.replace('.fa','_blastPDB.xml')))


def extract():
    for file in os.listdir('fasta'):
        if '_blast.xml' in file:
            raise Exception
