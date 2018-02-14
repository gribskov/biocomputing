"""=================================================================================================
beautiful soup blast

    
================================================================================================="""
import requests
#import lxml
from bs4 import BeautifulSoup

blast = 'https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi'

# program = 'blastp'
# database = 'pdb'
# # query = '''>AAG47671.1 ARV1 [Homo sapiens]
# # AMGNGGRSGCQYRCIECNQEAKELYRDYNHGVLKITICKSCQKPVDKYIEYDPVIILINAILCKAQAYRHILFNTQINIHGKLYLRWWQLQDSNQNTAPDDLIRYAKEWDF'''
# query = '''>At2
# NGGRSGCQYRCIECNQEAKELYRDYNHGVLKITICKSCQKPVDKYIEYDPVIILINAILCKAQAYRHILFNTQINIHGKLYLRWWQLQDSNQNTAPDDLIRYAKEWDF'''
#
# command = 'Put&PROGRAM={}&DATABASE={}&QUERY={}'.format(program, database, query)
# command = {'CMD': 'Put',
#            'PROGRAM': program,
#            'DATABASE': database,
#            'QUERY': query,
#            'EMAIL': 'gribskov@purdue.edu'
#            }
# print('command:', command)
# response = requests.post(blast, command)
# blast = BeautifulSoup(response.content, 'lxml')
# rid = blast.find('input',{'id':'rid'})
# rid = rid[0]
#print(rid)
print(rid['value'])

# print(blast.prettify())

# <input name="RID" type="hidden" value="88PNE34R014"/>
#      <input name="WWW_BLAST_TYPE" type="hidden" value=""/>
#       <!--QBlastInfoBegin
#     RID = 88PNE34R014
#     RTOE = 23
# QBlastInfoEnd
# -->

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

rid = '88SBTMYK014'
command = 'CMD=Get&&FORMAT_OBJECT=Searchinfo&RID={}'.format(rid)


