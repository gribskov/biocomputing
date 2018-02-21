"""=================================================================================================
blast with beautiful soup

    
================================================================================================="""
import requests
from bs4 import BeautifulSoup

blast = 'https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi'

program = 'blastp'
database = 'pdb'
query = '''>AAG47671.1 ARV1 [Homo sapiens]
AMGNGGRSGCQYRCIECNQEAKELYRDYNHGVLKITICKSCQKPVDKYIEYDPVIILINAILCKAQAYRHILFNTQINIHGKLYLRWWQLQDSNQNTAPDDLIRYAKEWDF'''

# command = 'Put&PROGRAM={}&DATABASE={}&QUERY={}'.format(program, database, query)
command = {'CMD': 'Put',
           'PROGRAM': program,
           'DATABASE': database,
           'QUERY': query,
           'EMAIL': 'gribskov@purdue.edu'
           }
print('command:', command)
response = requests.post(blast, command)

soup = BeautifulSoup(response.content, 'html.parser')
print(soup.prettify())
qblast = soup.find(string='RID = ')

print(qblast)
