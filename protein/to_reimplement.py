@_failsafe
def query_ELM_variant(self):
    # TODO complete
    mfile = os.path.join(self.settings.ELM_variant_folder, self.uniprot + '_' + self.file_friendly_mutation + '_ELM_variant.tsv')
    if os.path.isfile(mfile):
        self.log('Reading ELM variant data from file')
    else:
        self._assert_fetchable('ELM (variant')
        mseq = self.seq[0:self.resi - 1] + self.to_resn + self.seq[self.resi:]
        mseq_trim = mseq[max(self.resi - 30, 0):min(self.resi + 30, len(self.seq))]
        assert mseq_trim, 'Something went wrong in parsing {0}'.format(mseq)
        req = requests.get('http://elm.eu.org/start_search/{}'.format(mseq_trim))
        if req.status_code == 429:
            time.sleep(60)
            self.query_ELM()
            return
        elif req.status_code != 200:
            raise ConnectionError('Could not retrieve data: ' + req.text)
        response = req.text
        self.log('Retrieving ELM variant data from web')
        open(mfile, 'w').write(response)
    mdata = list(csv.DictReader(open(mfile, 'r'), delimiter='\t'))
    lbound = max(self.resi - 30, 0)
    ubound = min(self.resi + 30, len(self.seq))
    for entry in data:
        if int(entry['start']) < self.resi < int(entry['stop']):
            for m in range(len(mdata)):
                mentry = mdata[m]
                if entry['elm_identifier'] == mentry['elm_identifier'] and int(entry['start']) == int(
                        mentry['start']) + lbound:
                    del mdata[m]
                    self.mutational_effect.append(
                        'Possible linear motif <a target="_blank" href="http://elm.eu.org/elms/{0}">{0} <i class="fas fa-external-link-square"></i></a> spanning mutation: preserved'.format(
                            entry['elm_identifier']))
                    break
            else:
                self.mutational_effect.append(
                    'Possible linear motif <a target="_blank" href="http://elm.eu.org/elms/{0}">{0} <i class="fas fa-external-link-square"></i></a> spanning mutation: lost'.format(
                        entry['elm_identifier']))
    if mdata:
        for entry in mdata:
            if int(entry['start']) < 30 and int(entry['stop']) > 30:
                self.mutational_effect.append(
                    'Possible new linear motif <a target="_blank" href="http://elm.eu.org/elms/{0}">{0} <i class="fas fa-external-link-square"></i></a> spanning mutation'.format(
                        entry['elm_identifier']))






"""This chunk of code is for looking at the neighnbours
It will be implemented in to the code but not now.

"""

print('testing extended')
with open('data/human_prot_namedex.json') as f:
    namedex = json.load(f)

def get_friend(name):
    print(name)
    try:
        friend = Protein(uniprot=namedex[name], gene=name)
        friend.parse_uniprot()
        friend.parse_pLI()
        # friend.fetch_binders()
        return friend
    except Exception as err:
        print(err)
        return None

dock = get_friend('DOCK11')
dock.fetch_binders()
print(dock.partners)
with open('Dock11_test.csv', 'w', newline='') as fh:
    sheet = csv.DictWriter(fh, fieldnames='name uniprot uniprot_name group disease pLI pRec pNull'.split())
    sheet.writeheader()
    friends = set([f for l in dock.partners.values() for f in l])
    for friend in friends:
        groups = [g for g in dock.partners.keys() if friend in dock.partners[g]]
        fprot = get_friend(friend)
        if fprot:
            sheet.writerow({'name': friend,
                            'uniprot': fprot.uniprot,
                            'uniprot_name': fprot.uniprot_name,
                            'group': ' | '.join(groups),
                            'disease': ' | '.join([d['name'] for d in fprot.diseases]),
                            'pLI': fprot.pLI,
                            'pRec': fprot.pRec,
                            'pNull': fprot.pNull
                            })
        else:
            sheet.writerow({'name': friend})
