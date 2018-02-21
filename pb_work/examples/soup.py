"""=================================================================================================
beautiful soup blast

    
================================================================================================="""
import requests
# import lxml
from bs4 import BeautifulSoup, Comment

blast = 'https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi'

program = 'blastp'
database = 'pdb'
query = '''>Aaa
AKELYRDYNHGVLKITICKSCQKPVDKYIEYDPVIILINAILCKAQAYRHILFNTQINIHGKLYLRWWQLQDSNQNTAPDDLIRYAKEWDF'''

command = 'Put&PROGRAM={}&DATABASE={}&QUERY={}'.format(program, database, query)
command = {'CMD': 'Put',
           'PROGRAM': program,
           'DATABASE': database,
           'QUERY': query,
           'EMAIL': 'gribskov@purdue.edu'
           }
response = requests.post(blast, command)
submit = BeautifulSoup(response.content, 'html.parser')
rid = submit.find('input', {'id': 'rid'})
rid = rid['value']
print('RID:', rid)

# this is the polling loop

import time

maxtries = 10
notready = 1
# command = 'CMD=Get&FORMAT_OBJECT=Searchinfo&RID={}'.format(rid)
command = {'CMD': 'Get',
           'FORMAT_OBJECT': 'Searchinfo',
           'RID': rid
           }
print(blast,':',command)
status = ''
while notready:
    print('    polling ... try {}'.format(notready))
    response = requests.post(blast, command)
    check = BeautifulSoup(response.content, 'lxml')

    # complicated way to get all html comments <!-- to -->
    comments = check.find_all(string=lambda text: isinstance(text, Comment))

    for c in comments:
        if 'Status=' in c:
            if 'READY' in c:
                notready = False
                status = 'READY'
                break

    notready += 1
    if status == 'READY' or notready >= maxtries:
        break

    # don't poll too often, ncbi requests no more than 1/min
    time.sleep(60)

# end of polling loop

if not status == 'READY':
    # polling reached limit
    print('unable to find result () in {} tries'.format(rid, notready))

# final result

command = 'CMD=Get&FORMAT_TYPE=XML&RID={}'.format(rid)
response = requests.post(blast, command)
# print(response.text)

blast = BeautifulSoup(response.content, 'lxml')
hits = blast.find_all('hit')
print()
for hit in hits:
    hit_id = hit.find('hit_id')
    hit_len = hit.find('hit_len')
    evalue = hit.find('hsp_evalue')

    print('{}\t{}\t{}'.format(evalue.get_text(), hit_len.get_text(), hit_id.get_text()))
