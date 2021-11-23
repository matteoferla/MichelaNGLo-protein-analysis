import re
from collections import namedtuple
from .mutation import Mutation # for the aa3to1

ethnicities = ['afr', 'amr', 'asj', 'eas', 'fin', 'mid', 'nfe', 'oth', 'sas']  # no Amish
defaults = tuple([None] * 11 + [{e: 0 for e in ethnicities}])
aa3to1 = Mutation.aa3to1

class Variant(namedtuple('Variant',
                         ['id',  # rs number
                          'x', 'y', 'description',  # position
                          'from_residue',
                          'residue_index',
                          'to_residue',
                          'impact',
                          'homozygous',
                          'frequency',
                          'N',  # count and index are reserved
                          'consequence',
                          'frequencies'
                          ],
                         defaults=defaults)):
    """
    This should be generated the once... the complication comes from updating legacy data.

    ``x`` and ``y`` are the same as ``residue_index`` but are basically for the feature viewer.

    Stores the gnomAD/Clinvar data for easy use by FeatureViewer and co.

    NB. Clinvar ``description`` has the phenotypes and ``count`` has the number of submitters.

    A note that is likely now invalid: "Can be converted to Mutation."
    """
    ethnicities = ethnicities

    @property
    def type(self) -> str:
        """
        self.consequence was added later.
        Plus there "nonsense" is a bunch of more accurate labels.

        :return:
        """
        # 3.1:
        if self.to_residue:
            return 'missense' if self.to_residue in 'ACDEFGHIKLMNPQRSTVWY' else 'nonsense'
        else: # pre-3.1
            rex = re.match(r'(\w)(\d+)([\w*]+)', self.description)
            if not rex:
                return 'other'
            elif rex.group(3) in 'ACDEFGHIKLMNPQRSTVWY':
                return 'missense'
            else:
                return 'nonsense'

    @property
    def mutation(self) -> str:
        return f'{self.from_residue}{self.residue_index}{self.to_residue}'

    def to_dict(self):
        return dict(**self._asdict(),
                    type=self.type)

    @classmethod
    def from_v31_dict(cls, data: dict) -> 'Variant':
        if 'description' in data:
            # false alarm. it is not a v3.1
            # return cls(**data)
            raise NotImplementedError('Impossible.')
        i = data['residue_index']
        description = f"{data['from_residue']}{i}{data['to_residue']}"
        if data['identifier']:
            description += f" ({data['identifier']})"
            identifier = f"gnomAD_{i}_{i}_{data['identifier']}"
        else:
            identifier = f"gnomAD_{i}_{i}_xxxx"
        frequencies = {key.split()[1]: data[key] for key in data if 'frequency ' in key}
        return cls(id=identifier,
                    x=i,
                    y=i,
                   frequencies=frequencies,
                   description=description,
                   N=data['count'],
                   **{key: data[key] for key in ('consequence',
                                                 'from_residue',
                                                 'residue_index',
                                                 'to_residue',
                                                 'frequency',
                                                 'impact',
                                                 'homozygous')}
                   )

    @classmethod
    def from_clinvar_dict(cls, data: dict) -> 'Variant':
        i = data['residue_index']
        if 'RS' in data:
            identifier = f"clinvar_{i}_{i}_rs{data['RS']}"
        else:
            # none...
            identifier = f"clinvar_{i}_{i}_{data['identifier']}"
        # Not used: 'name':          'NM_014855.3(AP5Z1):c.1936G>A (p.Val646Met)'
        return cls(id=identifier,
                   x=i,
                   y=i,
                   frequency=0,
                   frequencies=None,  # blank!
                   description='; '.join(data['phenotype']),   # coopted.
                   N=data['N_sub'],  # coopted.
                   from_residue=aa3to1(data['from_residue']),
                   to_residue=aa3to1(data['to_residue']),
                   **{key: data[key] for key in (
                                                 'residue_index',
                                                 'impact',
                                                 'consequence')}
                   )
