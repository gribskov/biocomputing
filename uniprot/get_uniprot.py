'''
get_uniprot.py
retrieve family of sequences from uniprot based on id of the uniref90 representative
'''
import urllib.request, urllib.parse, urllib.error
import xml.etree.ElementTree as xml
import sqlite3 as db
from uniprot import *

id = 'K4A607_SETIT'
#id = 'K4A607_SETIT'
uniref_id = getUnirefIDByMember(90, id)
print('uniref:', uniref_id)

idlist = getUniprotIDByUniref(uniref_id)
print('idlist:', idlist)

# for id in idlist:
#     print( 'getting', id)
#     print(getUniprotByID(id,'txt'))
#
text = getUniProtByIDList(idlist,'xml')
# print(text)
uniprot = xml.fromstring(text)
# print('xml:', uniprot)
# print(uniprot.text)

# for entry in uniprot.iter('{http://uniprot.org/uniprot}entry'):
#     print('entryi:', entry.tag, entry.attrib)

ns = {'u':'http://uniprot.org/uniprot'}
entry = uniprot.findall('u:entry',ns)

for e in entry:
    # x = xml.dump(e)
    # x=xml.tostring(e,encoding='unicode',method='xml'
    x = xml.tostring(e, method='xml')
    print('entryf:', x.decode())