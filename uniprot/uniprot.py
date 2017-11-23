'''
uniprot access using REST api
'''
import urllib.request, urllib.parse, urllib.error
import requests

uniprot = 'http://www.uniprot.org/uniprot/'
uniref = 'http://www.uniprot.org/uniref/'
batch = 'http://www.uniprot.org/uploadlists/'

def getUnirefIDByMember(level, id):
    '''
    get the uniref90 group from the entry
    sample query: http://www.uniprot.org/uniref/?query=member:P99999+AND+identity:0.9
    usage
        text = getUnirefByMember( 'K4A607_SETIT', 90)
    '''
    if level == 90:
        level = 0.9
    else:
        level = 0.5

    query = '{0}?query=member:{1}+AND+identity:{2}&format={3}'.format(uniref, id, level, 'list')
    # print(query)
    fhand = urllib.request.urlopen(query)
    return fhand.readline().decode().strip()


def getUniprotIDByUniref(uniref_id):
    '''
    get the list of entries in the uniref group, store in idlist
    Generated query: http://www.uniprot.org/uniref/UniRef90_A0MWC0&format=list
    usage
        idlist = getUniprotIDByUniref(uniref_id)
    '''
    format = 'list'
    query = '{0}{1}&format={2}'.format(uniref, uniref_id, format)
    # print(query)
    fhand = urllib.request.urlopen(query)

    idlist = []
    for line in fhand:
        idlist.append(line.decode().strip())

    return idlist


def getUniprotByID(id, format):
    '''
    retrieve the whole uniprot entry for each member
    Allowed formats: txt, xml, rdf, fasta, tab
    usage
        text = getUniprotByID(id,'xml')
    '''
    query = '{0}{1}&format={2}'.format(uniprot, id, format)
    # print(query)
    fhand = urllib.request.urlopen(query)

    text = ''
    for line in fhand:
        text += line.decode()

    return text


def getUniProtByIDList(idlist):
    """
    Retrieves a list sequences from UniProt. See http://www.uniprot.org/faq/28 for building queries.
    usage
        text = getUniProtByIDList(idlist)
    """
    # convert idlist to string
    idstr = ''
    for id in idlist:
        idstr += id + ' '

    params = {
        'from': 'ACC+ID',
        'to': 'ACC',
        'format': 'txt',
        'uploadQuery': idstr
    }
    r = requests.post(batch, data=params)
    return r.text


if __name__ == '__main__':
    id = 'K4A607_SETIT'
    uniref_id = getUnirefIDByMember(90, id)
    print('uniref:', uniref_id)

    idlist = getUniprotIDByUniref(uniref_id)
    print('idlist:', idlist)

    # for id in idlist:
    #     print(getUniprotByID(id,'txt'))

    text = getUniProtByIDList(idlist)
    print(text)
