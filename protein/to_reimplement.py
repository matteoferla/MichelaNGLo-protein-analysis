
"""This chunk of code is for looking at the neighnbours
It will be implemented in to the code but not now.

"""


def junk():
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
