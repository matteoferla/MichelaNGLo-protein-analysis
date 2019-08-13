__description__ = """
The ET.Element gets expanded (monkey patched) with the following methods stored in ElementalExpansion before allocation.
The Element can be monkeypached by importing xmlschema, as opposed to using the __builtin__ workaround. Basically, I am piggybacking my monkeypatch on it, meaning that I don't need to copypaste from SO.

    ET.Element.ns_strip()  #returns the tag with the {namespace}
    ET.Element.is_tag(value) # boolean fx to check if tag == value
    ET.Element.describe # prints the content of the elemtn for debugging, similarly to dump but better.
    ET.Element.is_human() # boolean fx to check if human or dancer
    ET.Element.has_attr(key, opt_value) # boolean fx to check if it has key and optionally it key has the given value

To use: `from ET_monkeypatch import ET`
"""

import re

# Only the Python XML reader can be monkeypatched.
# needed to piggyback the monkeypatch. No schema validation used.
# copying method with xmlschema.
import sys, importlib
sys.modules.pop('xml.etree.ElementTree', None)
sys.modules['_elementtree'] = None
ET = importlib.import_module('xml.etree.ElementTree')
########

#### Expanding element tree element...

class NewElement(ET.Element):
    """
    This is a collection of methods that helps handle better the ET.Element instnaces. They are monkeypatched to the class object itself.
    """

    def ns_strip(self, ns='{http://uniprot.org/uniprot}'):
        return self.tag.replace(ns, '').replace('ns0', '').replace('{', '').replace('}', '')

    def is_tag(self, tag):
        if self.ns_strip() == tag:
            return True
        else:
            return False

    def describe(self):
        print('*' * 10)
        print('element', self)
        print('tag', self.tag)
        print('text', self.text)
        print('tail', self.tail)
        print('attr', self.attrib)
        print([child for child in self])
        print('*' * 10)

    def is_human(self):
        for elem in self:
            if elem.is_tag('organism'):
                for organism_el in list(elem):
                    if organism_el.text == 'Human':
                        return True
        else:
            return False

    def has_attr(self, key, value=None):
        if key in self.attrib:
            if not value:
                return True
            elif self.attrib[key] == value:
                return True
            elif isinstance(self.attrib[key], list) and value in self.attrib[key]:
                return True
        return False

    def get_attr(self, key):
        if key in self.attrib:
                return self.attrib[key]
        else:
            return ''

    def has_text(self):
        if not self.text:
            return False
        if re.match('\w', str(self.text)):
            return True
        else:
            return False

    def get_subtag(self, tag):
        """Gets only the first subtag."""
        for elem in self:
            if elem.is_tag(tag):
                return elem
        else:
            return None

    def get_sub_by_type(self, type_attr):
        """Gets only the first subtag."""
        for elem in self:
            if elem.has_attr('type', type_attr):
                return elem
        else:
            return None

setattr(ET,'Element', NewElement)


if __name__ == '__main__':
    xml = ET.fromstring('<hello type="test">World</hello>')
    print(type(ET.Element))
    print(type(xml))
    print(xml.is_tag('hello'))