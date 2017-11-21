'''
get_uniprot.py
retrieve family of sequences from uniprot based on id of the uniref90 representative
'''
import urllib.request, urllib.parse, urllib.error
from uniprot import *

id = 'K4A607_SETIT'
uniref_id = getUnirefIDByMember(90, id)
print('uniref:', uniref_id)

idlist = getUniprotIDByUniref(uniref_id)
print('idlist:', idlist)

for id in idlist:
    print( 'getting', id)
    print(getUniprotByID(id,'txt'))
