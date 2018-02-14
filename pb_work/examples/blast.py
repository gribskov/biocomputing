"""=================================================================================================
blast.py
Run blast search at ncbi vi blast API

    
================================================================================================="""
import requests

blast = 'https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi'

program = 'blastp'
database = 'pdb'
query = '''>AAG47671.1 ARV1 [Homo sapiens]
AMGNGGRSGCQYRCIECNQEAKELYRDYNHGVLKITICKSCQKPVDKYIEYDPVIILINAILCKAQAYRHILFNTQINIHGKLYLRWWQLQDSNQNTAPDDLIRYAKEWDF'''

# command = 'Put&PROGRAM={}&DATABASE={}&QUERY={}'.format(program, database, query)
# command = {'CMD': 'Put',
#            'PROGRAM': program,
#            'DATABASE': database,
#            'QUERY': query,
#            'EMAIL': 'gribskov@purdue.edu'
#            }
# print('command:', command)
# response = requests.post(blast, command)
# # print(response.url)
# # print('\n', response.text)
# # print('response:', response)
#
#
# info_key = 'QBlastInfoBegin\n    RID = '
# info_begin = response.text.find(info_key) + len(info_key)
# info_end = response.text.find(' ', info_begin)
# print(info_begin, info_end, )
# rid = response.text[info_begin:info_end]

rid = '7Y4XKJJS015'
command = 'CMD=Get&&FORMAT_OBJECT=Searchinfo&RID={}'.format(rid)

# import time
# maxtries = 10
# notready = 1
# while notready:
#     response = requests.post(blast, command)
#     # print('\n', response.text)
#
#     status_key = 'QBlastInfoBegin\n\tStatus='
#     status_begin = response.text.find(status_key) + len(status_key)
#     status_end = response.text.find('\n', status_begin)
#     # print(status_begin, status_end, )
#     status = response.text[status_begin:status_end]
#     print(status)
#     if status =='READY':
#         notready = 0
#         break
#     else:
#         notread += 1
#     if notready >= maxtries:
#         break
#
#     # don't poll too often, ncbi requests no more than 1/min
#     time.sleep(60)
#
# if notready > 0:
#     #polling reached lime
#     print('unable to find result () in {} tries'.format(rid, notready))

# get the final result

command = 'CMD=Get&FORMAT_TYPE=XML&RID={}'.format(rid)
response = requests.post(blast, command)
print(response.text)
