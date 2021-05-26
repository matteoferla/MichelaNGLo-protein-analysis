import re
from collections import namedtuple

class Variant(namedtuple('Variant',
                         ['id', 'x', 'y', 'impact', 'description', 'homozygous'],
                         defaults=(None, None, None, None, None, None))):
    """
    Stores the gnomAD data for easy use by FeatureViewer and co. Can be converted to Mutation.
    """
    @property
    def type(self) -> str:
        rex = re.match(r'(\w)(\d+)([\w*]+)', self.description)
        if not rex:
            return 'other'
        elif rex.group(3) in 'ACDEFGHIKLMNPQRSTVWY':
                return 'missense'
        else:
            return 'nonsense'

    def to_dict(self):
        return dict(id=self.id,
                    x=self.x, y=self.y,
                    impact=self.impact,
                    type=self.type,
                    description=self.description,
                    homozygous=self.homozygous)
