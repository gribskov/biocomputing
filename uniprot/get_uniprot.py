'''
get_uniprot.py
retrieve family of sequences from uniprot based on id of the uniref90 representative
'''
import urllib.request, urllib.parse, urllib.error

uniprot = 'http://www.uniprot.org/uniprot/'
uniref = 'http://www.uniprot.org/uniref/'
id = 'K4A607_SETIT'
format = 'list'

# get the uniref90 group from the entry
# http://www.uniprot.org/uniref/?query=member:P99999+AND+identity:0.9

query = '{0}?query=member:{1}+AND+identity:0.9&format={2}'.format(uniref, id, format)
print(query)
fhand = urllib.request.urlopen(query)
uniref_id = fhand.readline().decode().strip()
print('uniref:', uniref_id)

# get the list of entries in the uniref group, store in idlist
#  http://www.uniprot.org/uniref/UniRef90_A0MWC0&format=list

query = '{0}{1}&format={2}'.format(uniref, uniref_id, format)
print(query)
fhand = urllib.request.urlopen(query)

idlist = []
for line in fhand:
    idlist.append(line.decode().strip())

print('idlist:', idlist)

# retrieve the whole uniprot entry for each member
format = 'txt'
for id in idlist:
    query = '{0}{1}&format={2}'.format(uniprot, id, format)
    print(query)
    fhand = urllib.request.urlopen(query)

    for line in fhand:
        print(line.decode().strip())
