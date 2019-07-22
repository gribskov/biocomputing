"""=================================================================================================
parsing xml with lxml


================================================================================================="""
#from lxml import etree
from io import StringIO

blast = etree.parse('hemoglobin_blast.xml')
print(blast)
#
# blast = etree.ElementTree(file='blast.xml')
# print(blast)

# blastxml = open('blast.xml', 'r')
# xml = blastxml.read()
# blast = etree.parse(StringIO(xml))
#
# print(etree.tostring(blast,pretty_print=True))
# print(blast)
# for keys in blast.attrib:
#     print('att:{}   val:{}'.format(key, blast.attrib[key]))

from lxml import etree
blast = etree.parse('hemoglobin_blast.xml')

hits = blast.xpath('//Hit')
print(hits)
for hit in hits:
    id = hit.xpath('Hit_id')
    print('\n{}'.format(id[0].text))

    hsps = hit.xpath('//Hsp_evalue')
    for hsp in hsps:
        print('    {}'.format(hsp.text))