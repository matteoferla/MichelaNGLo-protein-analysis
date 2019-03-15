def old_retrieve_references(self):
    print('*' * 20)
    print('CORE REFERENCE DATA IS MISSING')
    print('There are two options, you have never ever run this script before or the folder {0} is not corrent'.format(self.reference_folder))
    print('this is super experimental (i.e. I\'ve never bother')
    i = input('Continue y/[n] _')
    if not i or i in ('N', 'n'):
        print('Exiting...')
        exit()
    addresses = ('ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz')
    for url in addresses:
        file = os.path.split(url)[1]
        if os.path.isfile(os.path.join(self.raw_folder, file)):
            continue
        else:
            req = requests.get(url)  # headers={"Accept": "application/xml"}
            if req.status_code != 200:
                raise ConnectionError('Could not retrieve data: ' + req.text)
            data = req.text
            open(file, 'w').write(data)
            if os.path.splitext(url) == 'gz':
                raise Exception('Okay. this will not work on windows... Can you unzip this file?? cd {0}; tar -x {1}; cd ../..'.format(self.reference_folder, file))
        raise NotImplementedError('Due to crappy windows 8 computer... this part is manual in VM: cat *.psi > cat.psi where psi files are from http://interactome.baderlab.org/data/')
