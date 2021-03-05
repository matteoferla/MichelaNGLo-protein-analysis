"""
The splitting of gnomAD data is the slowest step. On my 2015 core i5 desktop it takes 30 secs per gene (3s on a single thread on the cluster).
Therefore, this script is a hack to run it on n cores.
NB. The Uniprot splitter runs on threads but within one core.
If you set an env variable 'SLACK_WEBHOOK', you get messages.
Annoyingly, ``taxid9606-names2uniprot.json`` is required.
"""

cores = 25

# Get the following. Generally servers throttle https requests. So the general strategy is to do this locally and rsync.
'''wget https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.exome_calling_intervals.sites.vcf.bgz;
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.1.vcf.bgz &
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.2.vcf.bgz;
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.3.vcf.bgz &
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.4.vcf.bgz;
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.5.vcf.bgz &
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.6.vcf.bgz;
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.7.vcf.bgz &
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.8.vcf.bgz;
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.9.vcf.bgz &
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.10.vcf.bgz;
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.11.vcf.bgz &
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.12.vcf.bgz;
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.13.vcf.bgz &
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.14.vcf.bgz;
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.15.vcf.bgz &
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.16.vcf.bgz;
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.17.vcf.bgz &
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.18.vcf.bgz;
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.19.vcf.bgz &
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.20.vcf.bgz;
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.21.vcf.bgz;
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.22.vcf.bgz;
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.23.vcf.bgz;
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.X.vcf.bgz;
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.Y.vcf.bgz;'''

################

def gwrap(file):
    from michelanglo_protein.generate.split_gnomAD import gnomAD
    import os, requests, re, unicodedata
    def message(msg):
        if 'SLACK_WEBHOOK' in os.environ:
            msg = unicodedata.normalize('NFKD', msg).encode('ascii', 'ignore').decode('ascii')
            msg = re.sub('[^\w\s\-.,;?!@#()\[\]]', '', msg)
            r = requests.post(url=os.environ['SLACK_WEBHOOK'],
                              headers={'Content-type': 'application/json'},
                              data=f"{{'text': '{msg}'}}")
        else:
            print(msg)
    gnomAD(masterfiles=[file],
           namedexfile='taxid9606-names2uniprot.json',
           folder='gnomAD').split()
    message(f'{file} parsed!')


###################
import os
from multiprocessing import Pool

p = Pool(cores)
p.map(gwrap, [file for file in os.listdir() if 'gnomad.' in file])
