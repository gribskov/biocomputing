"""=================================================================================================
eutil.py
retrieve information from NCBI  suing eutil REST interface
works 2 March 2021
================================================================================================="""
import requests

# https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=protein&term=asthma&usehistory=y

ncbi = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
esearch = 'esearch.fcgi?'

query = 'db=protein&term=asthma&usehistory=y'
response = requests.get(ncbi + esearch + query)
print(response.content)

print('\n', response.text)

params = {'db': 'protein',
          'term': 'asthma',
          'usehistory': 'y'}

esearch = 'esearch.fcgi'
response = requests.post(ncbi + esearch, params)
print('\n', response.text)
